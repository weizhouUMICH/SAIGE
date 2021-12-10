
// This includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
// std::this_thread::sleep_for (std::chrono::seconds(1));
#include <cstdio>         // std::remove

// Currently, omp does not work well, will check it later
// error: SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
// remove all Rcpp::List to check if it works
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]]

#include "Main.hpp"
#include "PLINK.hpp"
#include "BGEN.hpp"
#include "VCF.hpp"
#include "UTIL.hpp"
#include "DenseGRM.hpp"

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
static VCF::VcfClass* ptr_gVCFobj = NULL;

// global objects for different analysis methods
//single, SAIGE
//Region, SAIGE-GENE+


// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop", //drop is not allowed
double g_missingRate_cutoff;
unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_minMAC_cutoff;    // for Rare Variants (RVs) whose MAC < this value, we aggregate these variants like SAIGE-GENE+ 
double g_region_maxMAF_cutoff;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage

arma::uvec g_group;
bool g_ifOutGroup;
unsigned int g_nGroup;

// global variables for sparse GRM
arma::sp_mat g_SparseGRM;

// [[Rcpp::export]]
void setSparseGRMInCPP(Rcpp::List t_KinMatListR)
{
  arma::umat locations = t_KinMatListR["locations"];
  arma::vec values = t_KinMatListR["values"];
  int n = t_KinMatListR["nSubj"];
  // make a sparse matrix
  arma::sp_mat KinMat(locations, values, n, n);
  g_SparseGRM = KinMat;
}

// [[Rcpp::export]]
void setDenseGRMInCPP(double t_memoryChunk,
                      double t_minMafGRM,
                      double t_maxMissingGRM)
{
  ptr_gDenseGRMobj = new DenseGRM::DenseGRMClass(ptr_gPLINKobj, t_memoryChunk, t_minMafGRM, t_maxMissingGRM);
}

// [[Rcpp::export]]
arma::vec getDenseGRMInCPP(arma::vec t_bVec,
                           std::string t_excludeChr, 
                           int t_grainSize)
{
  arma::vec yVec = DenseGRM::getKinbVec(t_bVec, ptr_gDenseGRMobj, t_excludeChr, t_grainSize);
  return yVec;
}

// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               unsigned int t_omp_num_threads,
                               arma::uvec t_group,
                               bool t_ifOutGroup,
                               unsigned int t_nGroup)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_omp_num_threads = t_omp_num_threads;
  g_group = t_group;
  g_ifOutGroup = t_ifOutGroup;
  g_nGroup = t_nGroup;
}

// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_max_maf_region,
                               double t_min_mac_region,
                               unsigned int t_max_markers_region,
                               unsigned int t_omp_num_threads,
                               arma::uvec t_group,
                               bool t_ifOutGroup,
                               unsigned int t_nGroup)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_region_minMAC_cutoff = t_min_mac_region;
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_region_maxMarkers_cutoff = t_max_markers_region;
  g_omp_num_threads = t_omp_num_threads;
  g_group = t_group;
  g_ifOutGroup = t_ifOutGroup;
  g_nGroup = t_nGroup;
}

void updateGroupInfo(arma::vec t_GVec,
                     std::vector<uint32_t> t_indexForMissing,
                     arma::vec& nSamplesInGroupVec,
                     arma::vec& AltCountsInGroupVec,
                     arma::vec& AltFreqInGroupVec)
{
  unsigned int n1 = t_GVec.size();
  nSamplesInGroupVec.zeros();
  AltCountsInGroupVec.zeros();
  
  unsigned int i1 = 0;
  for(unsigned int i = 0; i < n1; i++){
    if(i == t_indexForMissing.at(i1)){
      if(i1 < t_indexForMissing.size() - 1)
        i1 ++;
    }else{
      unsigned int grp = g_group.at(i);
      
      // std::cout << "grp:\t" << grp << std::endl;
      // std::cout << "i:\t" << i << std::endl;
      // std::cout << "i1:\t" << i1 << std::endl;
      
      nSamplesInGroupVec.at(grp) += 1;
      AltCountsInGroupVec.at(grp) += t_GVec.at(i);
    }
  }
  
  AltFreqInGroupVec = AltCountsInGroupVec / nSamplesInGroupVec / 2;
}

//////// ---------- Main function for marker-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainMarkerInCPP(std::string t_method,       // "SAIGE"
                           std::string t_genoType,     // "PLINK", "BGEN", "VCF"
                           std::vector<uint32_t> t_genoIndex)  
{
  int q = t_genoIndex.size();  // number of markers
  
  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> missingRateVec(q);  // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);       
  std::vector<double> pvalVec(q, arma::datum::nan);
  std::vector<double> zScoreVec(q, arma::datum::nan);
  
  arma::mat nSamplesInGroup;
  arma::mat AltCountsInGroup;
  arma::mat AltFreqInGroup;

  if(g_ifOutGroup){
    nSamplesInGroup.resize(q, g_nGroup);
    AltCountsInGroup.resize(q, g_nGroup);
    AltFreqInGroup.resize(q, g_nGroup);
  }

  // std::cout << "Totally " << g_omp_num_threads << " thread(s) were used for parallel computation." << std::endl;
  
  // loop for all markers
//   omp_set_dynamic(0);     // Explicitly disable dynamic teams
//   omp_set_num_threads(g_omp_num_threads); // Use 4 threads for all consecutive parallel regions
//   
// #pragma omp parallel
// {
  for(int i = 0; i < q; i++){
    
    if(i % 1000 == 0){
      std::cout << "Completed " << i << "/" << q << " markers in the chunk." << std::endl;
    }
    
    // information of marker
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    
    uint32_t gIndex = t_genoIndex.at(i);
    
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    int n = GVec.size();
    
    // std::cout << "test1.1" << std::endl;
    
    if(g_ifOutGroup){
      arma::vec nSamplesInGroupVec(g_nGroup);
      arma::vec AltCountsInGroupVec(g_nGroup);
      arma::vec AltFreqInGroupVec(g_nGroup);
      
      // std::cout << "test1.2" << std::endl;
      // std::cout << "g_nGroup:\t" << g_nGroup << std::endl;
      
      updateGroupInfo(GVec, indexForMissing, nSamplesInGroupVec, AltCountsInGroupVec, AltFreqInGroupVec);
      
      // std::cout << "test1.3" << std::endl;
      
      nSamplesInGroup.row(i) = nSamplesInGroupVec.t();
      AltCountsInGroup.row(i) = AltCountsInGroupVec.t();
      AltFreqInGroup.row(i) = AltFreqInGroupVec.t();
    }
    
    // e.g. 21:1000234:A:T
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    
    // record basic information for the marker
    markerVec.at(i) = marker;               // marker IDs
    infoVec.at(i) = info;    // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    
    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * n * (1 - missingRate);
    
    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff))
      continue;
    
    // Check UTIL.cpp
    flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
    
    // analysis results for single-marker
    double Beta, seBeta, pval, zScore;
    
    Unified_getMarkerPval(t_method, GVec, 
                          false, // bool t_isOnlyOutputNonZero, 
                          indexForNonZero, Beta, seBeta, pval, zScore, altFreq);
    
    BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    seBetaVec.at(i) = seBeta;       
    pvalVec.at(i) = pval;
    zScoreVec.at(i) = zScore;
  }

  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("altCountsVec") = altCountsVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("pvalVec") = pvalVec,
                                          Rcpp::Named("beta") = BetaVec,
                                          Rcpp::Named("seBeta") = seBetaVec,
                                          Rcpp::Named("zScore") = zScoreVec,
                                          Rcpp::Named("nSamplesInGroup") = nSamplesInGroup,
                                          Rcpp::Named("AltCountsInGroup") = AltCountsInGroup,
                                          Rcpp::Named("AltFreqInGroup") = AltFreqInGroup);
  
  return OutList;  
}

//////// ---------- Main function for region-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainRegionURVInCPP(std::string t_method,       // "POLMM", "SPACox", "SAIGE" (to be continued)
                              std::string t_genoType,     // "PLINK", "BGEN"
                              std::vector<uint32_t> t_genoIndex,
                              unsigned int t_n)           // sample size
{
  unsigned int q = t_genoIndex.size();                 // number of Ultra-Rare Variants (URV) markers (after QC) in one region
  double Stat, Beta, seBeta, pval0, pval1;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  
  arma::vec GVecURV(t_n, arma::fill::zeros);    // aggregate ultra-rare variants (URV) whose MAC less than cutoff (g_region_minMAC_cutoff)
  
  for(unsigned int i = 0; i < q; i++)
  {
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;

    uint32_t gIndex = t_genoIndex.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
    
    if(altFreq < 0.5){
      GVecURV = arma::max(GVecURV, GVec);  // edited on 2021-08-20
    }else{
      GVecURV = arma::max(GVecURV, 2 - GVec); // edited on 2021-08-20
    }
  }
  
  Unified_getRegionPVec(t_method, GVecURV, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("Stat") = Stat,
                                          Rcpp::Named("pval1") = pval1);

  return OutList;  
}

// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(std::string t_method,       // "POLMM", "SPACox", "SAIGE" (to be continued)
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint32_t> t_genoIndex,
                           std::string t_outputFile,
                           unsigned int t_n,           // sample size  
                           arma::mat P1Mat,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
                           arma::mat P2Mat)
{
  unsigned int q = t_genoIndex.size();                 // number of markers (before QC) in one region
  
  // set up output (Ultra-Rare Variants, URV)  removed on 09-18-2021
  // std::vector<std::string> markerVec(q), markerURVVec(q);      // marker IDs
  // std::vector<std::string> infoVec(q), infoURVVec(q);          // marker information: CHR:POS:REF:ALT
  // std::vector<double> altFreqVec(q), altFreqURVVec(q);         // allele frequencies of the ALT allele, this is not always < 0.5.
  // std::vector<double> MACVec(q), MACURVVec(q);
  // std::vector<double> MAFVec(q), MAFURVVec(q);
  // std::vector<double> missingRateVec(q), missingRateURVVec(q); // missing rate
  
  // added on 09-18-2021
  arma::uvec indicatorVec(q, arma::fill::zeros);       // 0: does not pass QC, 1: non-URV, 2: URV
  Rcpp::StringVector markerVec(q);
  Rcpp::StringVector infoVec(q);
  arma::vec altFreqVec(q);         // allele frequencies of the ALT allele, this is not always < 0.5.
  arma::vec MACVec(q);
  arma::vec MAFVec(q);
  arma::vec missingRateVec(q);     // missing rate
  
  std::vector<double> BetaVec(q);            // beta value for ALT allele
  std::vector<double> seBetaVec(q);          // seBeta value
  std::vector<double> pval0Vec(q);           // p values from normal distribution approximation  // might be confused, is this needed?
  std::vector<double> pval1Vec(q);           // p values from more accurate methods including SPA and ER
  
  std::vector<double> StatVec(q);            // score statistics

  arma::mat nSamplesInGroup;
  arma::mat AltCountsInGroup;
  arma::mat AltFreqInGroup;
  
  if(g_ifOutGroup){
    nSamplesInGroup.resize(q, g_nGroup);
    AltCountsInGroup.resize(q, g_nGroup);
    AltFreqInGroup.resize(q, g_nGroup);
  }
  
  // example #1: (q = 999, m1 = 10) -> (nchunks = 100, m2 = 9)
  // example #2: (q = 1000, m1 = 10) -> (nchunks = 100, m2 = 10)
  // example #3: (q = 1001, m1 = 10) -> (nchunks = 101, m2 = 1)
  
  unsigned int m1 = g_region_maxMarkers_cutoff;     // number of markers in all chunks expect for the last chunk
  // P1Mat should be of dimension: m1 * t_n
  // P2Mat should be of dimension: t_n * m1
  
  // Suppose that 
  // n is the sample size in analysis 
  // m (<q) is the number of markers that pass the marker-level QC (e.g., g_missingRate_cutoff and g_region_maxMAF_cutoff)
  // VarMat (m x m) is the variance matrix of these m markers
  // VarMat = P1Mat %*% P2Mat, where P1Mat is of (m x n) and P2Mat is of (n x m)
  
  // Added on 09-17-2021: we collapse all ultra-rare variants (URV) to get one "fake" marker. 
  // That part has been moved to function mainRegionURVInCPP()
  
  std::vector<unsigned int> mPassCVVec;
    
  // conduct marker-level analysis
  double Stat, Beta, seBeta, pval0, pval1;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  
  // initiate chunk information
  unsigned int nchunks = 0;
  unsigned int ichunk = 0;
  unsigned int i1InChunk = 0;
  unsigned int i1 = 0;    // index of Markers (non-URV)
  unsigned int i2 = 0;    // index of Markers (Ultra-Rare Variants, URV)
  
  // cycle for q markers
  for(unsigned int i = 0; i < q; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    
    uint32_t gIndex = t_genoIndex.at(i);
    
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    
    flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
    
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
    
    // Quality Control (QC)
    if((missingRate > g_missingRate_cutoff) || (MAF > g_region_maxMAF_cutoff) || MAF == 0){
      continue;  // does not pass QC
    }

    if(g_ifOutGroup){
      arma::vec nSamplesInGroupVec(g_nGroup);
      arma::vec AltCountsInGroupVec(g_nGroup);
      arma::vec AltFreqInGroupVec(g_nGroup);
      
      // std::cout << "test1.2" << std::endl;
      // std::cout << "g_nGroup:\t" << g_nGroup << std::endl;
      
      updateGroupInfo(GVec, indexForMissing, nSamplesInGroupVec, AltCountsInGroupVec, AltFreqInGroupVec);
      
      // std::cout << "test1.3" << std::endl;
      
      nSamplesInGroup.row(i) = nSamplesInGroupVec.t();
      AltCountsInGroup.row(i) = AltCountsInGroupVec.t();
      AltFreqInGroup.row(i) = AltFreqInGroupVec.t();
    }
    
    markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    
    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants
      
      indicatorVec.at(i) = 1;
      
      if(i1InChunk == 0){
        std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
      }
      
      // markerVec.at(i1) = marker;             // marker IDs
      // infoVec.at(i1) = info;                 // marker information: CHR:POS:REF:ALT
      // altFreqVec.at(i1) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
      // missingRateVec.at(i1) = missingRate;
      // MACVec.at(i1) = MAC;
      // MAFVec.at(i1) = MAF;
      
      Unified_getRegionPVec(t_method, GVec, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
      
      // insert results to pre-setup vectors and matrix
      StatVec.at(i1) = Stat;        
      
      // BetaVec.at(i1) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true       
      // seBetaVec.at(i1) = seBeta;       
      // pval0Vec.at(i1) = pval0;
      BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true       
      seBetaVec.at(i) = seBeta;       
      pval0Vec.at(i) = pval0;
      
      pval1Vec.at(i1) = pval1;
      // adjPVec.at(i1) = pval1;
      
      P1Mat.row(i1InChunk) = P1Vec.t();
      P2Mat.col(i1InChunk) = P2Vec;
      
      i1 += 1;
      i1InChunk += 1;
      
    }else{   // Ultra-Rare Variants (URV)
      
      indicatorVec.at(i) = 2;
      
      // markerURVVec.at(i2) = marker;             // marker IDs
      // infoURVVec.at(i2) = info;                 // marker information: CHR:POS:REF:ALT
      // altFreqURVVec.at(i2) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
      // missingRateURVVec.at(i2) = missingRate;
      // MACURVVec.at(i2) = MAC;
      // MAFURVVec.at(i2) = MAF;
        
      i2 += 1;
    }
    
    if(i1InChunk == m1){
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      
      mPassCVVec.push_back(m1);
      ichunk += 1;
      i1InChunk = 0;
    }
    
    Rcpp::checkUserInterrupt();
  }
  
  if(i1 == 0){
    std::cout << "Only ultra-rare variants are found. This region will be skipped." << std::endl;
    Rcpp::List OutList = Rcpp::List::create();
    return OutList;
  }
  
  nchunks = ichunk + 1;
  arma::mat VarMat(i1, i1);
  
  // non Ultra Rare Variants
  // markerVec.resize(i1);
  // infoVec.resize(i1);              // marker information: CHR:POS:REF:ALT
  // altFreqVec.resize(i1);           // allele frequencies of ALT allele, this is not always < 0.5.
  // missingRateVec.resize(i1);
  // MACVec.resize(i1);
  // MAFVec.resize(i1);
  StatVec.resize(i1);        
  // BetaVec.resize(i1);              // Beta if flip = false, -1 * Beta is flip = true       
  // seBetaVec.resize(i1);       
  // pval0Vec.resize(i1);
  pval1Vec.resize(i1);
  // adjPVec.resize(i1);
  
  // Ultra Rare Variants
  // markerURVVec.resize(i2);          // marker IDs
  // infoURVVec.resize(i2);            // marker information: CHR:POS:REF:ALT
  // altFreqURVVec.resize(i2);         // allele frequencies of ALT allele, this is not always < 0.5.
  // missingRateURVVec.resize(i2);
  // MACURVVec.resize(i2);
  // MAFURVVec.resize(i2);
  
  mPassCVVec.push_back(i1InChunk);

  if(i1InChunk != 0){
    P1Mat = P1Mat.rows(0, i1InChunk - 1);
    P2Mat = P2Mat.cols(0, i1InChunk - 1);
    if(nchunks != 1){
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    }
  }
  
  // not so many markers in the region, so all matrix is in memory
  if(nchunks == 1){
    VarMat = P1Mat * P2Mat;
  }

  // the region includes more markers than limitation, so P1Mat and P2Mat have been put in hard drive
  if(nchunks > 1)
  {
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;
    
    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      last_row = first_row + mPassCVVec.at(index1) - 1;
      
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      
      P1Mat.load(P1MatFile);
      
      if(P1Mat.n_cols == 0) continue;
      
      // off-diagonal sub-matrix
      for(unsigned int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;
        
        P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + ".bin");
        
        if(P2Mat.n_cols == 0) continue;
        
        arma::mat offVarMat = P1Mat * P2Mat;
        
        // last_col = first_col + mPassQCVec.at(index2) - 1;
        last_col = first_col + mPassCVVec.at(index2) - 1;
        
        VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
        VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
        
        first_col = last_col + 1;
      }
      
      // diagonal sub-matrix
      last_col = first_col + mPassCVVec.at(index1) - 1;
      std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
      P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin");
      
      arma::mat diagVarMat = P1Mat * P2Mat;
      
      VarMat.submat(first_row, first_col, last_row, last_col) = diagVarMat;
      
      first_row = last_row + 1;
      first_col = 0;
      Rcpp::checkUserInterrupt();
    }
    
    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      std::string P2MatFile = t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin";
      const char* File1 = P1MatFile.c_str();
      const char* File2 = P2MatFile.c_str();
      std::remove(File1);
      std::remove(File2);
    }
  }
  
  // calculate p-values for the burden test
  
  // To be added later (2021-09-17)
  // double pval0Burden, pval1Burden;
  // Unified_getRegionPVec(t_method, GVec, Beta, seBeta, pval, P1Vec, P2Vec);
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("VarMat") = VarMat,
                                          Rcpp::Named("indicatorVec") = indicatorVec,
                                          Rcpp::Named("markerVec") = markerVec,
                                          // Rcpp::Named("markerURVVec") = markerURVVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          // Rcpp::Named("infoURVVec") = infoURVVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          // Rcpp::Named("altFreqURVVec") = altFreqURVVec,
                                          Rcpp::Named("MAFVec") = MAFVec,
                                          // Rcpp::Named("MAFURVVec") = MAFURVVec,
                                          Rcpp::Named("MACVec") = MACVec,
                                          // Rcpp::Named("MACURVVec") = MACURVVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          // Rcpp::Named("missingRateURVVec") = missingRateURVVec,
                                          Rcpp::Named("StatVec") = StatVec,
                                          Rcpp::Named("beta") = BetaVec,
                                          Rcpp::Named("seBeta") = seBetaVec,
                                          Rcpp::Named("PvalueNorm") = pval0Vec,  // If this line is uncommented, then error comes up, maybe reach a number limit?
                                          Rcpp::Named("pval1Vec") = pval1Vec,
                                          Rcpp::Named("nSamplesInGroup") = nSamplesInGroup,
                                          Rcpp::Named("AltCountsInGroup") = AltCountsInGroup,
                                          Rcpp::Named("AltFreqInGroup") = AltFreqInGroup);
  
  return OutList;
}


//////// ---------- Main function for genotype extraction --------- ////////////

// [[Rcpp::export]]
arma::mat getGenoInCPP(std::string t_genoType,
                       Rcpp::DataFrame t_markerInfo,
                       int n,
                       std::string t_imputeMethod)
{
  int q = t_markerInfo.nrow();         // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  arma::mat GMat(n, q);
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType,          // "PLINK", "BGEN"
                                          gIndex,              // different meanings for different genoType
                                          ref,                 // REF allele
                                          alt,                 // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,              // marker ID extracted from genotype file
                                          pd,                  // base position
                                          chr,                 // chromosome
                                          altFreq,             // frequency of ALT allele
                                          altCounts,           // counts of ALT allele
                                          missingRate,         // missing rate
                                          imputeInfo,          // imputation information score, i.e., R2 (all 1 for PLINK)
                                          true,                // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  // check UTIL.cpp
    GMat.col(i) = GVec;
  }
  
  return GMat;
}

// [[Rcpp::export]]
arma::sp_mat getSpGenoInCPP(std::string t_genoType,
                            Rcpp::DataFrame t_markerInfo,
                            int n,
                            std::string t_imputeMethod)
{
  int q = t_markerInfo.nrow();         // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  arma::sp_mat GMat(n, q);             // change #1 compared to getGenoInCPP()
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType,    // "PLINK", "BGEN"
                                          gIndex,        // different meanings for different genoType
                                          ref,           // REF allele
                                          alt,           // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,        // marker ID extracted from genotype file
                                          pd,            // base position
                                          chr,           // chromosome
                                          altFreq,       // frequency of ALT allele
                                          altCounts,     // counts of ALT allele
                                          missingRate,   // missing rate
                                          imputeInfo,    // imputation information score, i.e., R2 (all 1 for PLINK)
                                          true,         // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  // check UTIL.hpp
    GMat.col(i) = arma::sp_mat(GVec); // change #2 compared to getGenoInCPP()
  }
  
  return GMat;
}

// a unified function to get single marker from genotype file
arma::vec Unified_getOneMarker(std::string t_genoType,   // "PLINK", "BGEN"
                               uint64_t t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint32_t>& t_indexForMissing,     // index of missing genotype data
                               bool t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint32_t>& t_indexForNonZero)     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
  arma::vec GVec;
  if(t_genoType == "PLINK"){
    GVec = ptr_gPLINKobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                       t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                       true);   // t_isTrueGenotype, only used for PLINK format.
  }
  
  if(t_genoType == "BGEN"){
    bool isBoolRead;
    GVec = ptr_gBGENobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead);
  }
  
  return GVec;
}

// a unified function to get marker-level p-value
void Unified_getMarkerPval(std::string t_method,   // "POLMM", "SPACox", "SAIGE"
                           arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq)
{
  if(t_method == "SAIGE"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAIGE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    
    ptr_gSAIGEobj->getMarkerPval(t_GVec, t_Beta, t_seBeta, t_pval, t_altFreq, t_zScore);
  }
  
}

// a unified function to get marker-level information for region-level analysis

// Unified_getRegionPVec(t_method, GVec, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
void Unified_getRegionPVec(std::string t_method, 
                           arma::vec t_GVec, 
                           double& t_Stat,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval0, 
                           double& t_pval1,
                           arma::vec& t_P1Vec, 
                           arma::vec& t_P2Vec)
{
  // Check src/POLMM.cpp and src/POLMM.hpp
  if(t_method == "POLMM")
  {
    ptr_gPOLMMobj->getRegionPVec(t_GVec, t_Stat, t_Beta, t_seBeta, t_pval0, t_pval1, t_P1Vec, t_P2Vec);
  }
  
  // Check src/SPACox.cpp and src/SPACox.hpp
  if(t_method == "SPACox")
  {
    ptr_gSPACoxobj->getRegionPVec(t_GVec, t_Stat, t_pval0, t_pval1, t_P1Vec, t_P2Vec);
  }
  
}

//////// ---------- Main functions to set objects for different genotype format --------- ////////////

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> t_SampleInModel,
                      std::string t_AlleleOrder)
{
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                        t_famFile,
                                        t_bedFile,
                                        t_SampleInModel,
                                        t_AlleleOrder);
  
  int n = ptr_gPLINKobj->getN();
  std::cout << "Number of samples:\t" << n << std::endl;
  
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> t_SampleInBgen,
                     std::vector<std::string> t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropmissingdosagesInBgen,
                     std::string t_AlleleOrder)
{
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                     t_bgenFileIndex,
                                     t_SampleInBgen,
                                     t_SampleInModel,
                                     t_isSparseDosageInBgen,
                                     t_isDropmissingdosagesInBgen,
                                     t_AlleleOrder);
  int n = ptr_gBGENobj->getN();
  std::cout << "Number of samples:\t" << n << std::endl;
}


//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

// [[Rcpp::export]]
void setPOLMMobjInCPP(arma::mat t_muMat,
                      arma::mat t_iRMat,
                      arma::mat t_Cova,
                      arma::uvec t_yVec,
                      double t_tau,
                      bool t_printPCGInfo,
                      double t_tolPCG,
                      int t_maxiterPCG,
                      double t_varRatio, 
                      double t_SPA_cutoff,
                      bool t_flagSparseGRM)
{
  // check POLMM.cpp
  ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
                                        t_iRMat,
                                        t_Cova,
                                        t_yVec,
                                        g_SparseGRM,
                                        t_tau,
                                        t_printPCGInfo,
                                        t_tolPCG,
                                        t_maxiterPCG,
                                        t_varRatio, 
                                        t_SPA_cutoff,
                                        t_flagSparseGRM);
}

// [[Rcpp::export]]
Rcpp::List setPOLMMobjInCPP_NULL(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                                 arma::mat t_Cova,
                                 arma::uvec t_yVec,     // should be from 0 to J-1
                                 arma::vec t_beta,
                                 arma::vec t_bVec,
                                 arma::vec t_eps,           // 
                                 double t_tau,
                                 Rcpp::List t_SPmatR,    // output of makeSPmatR()
                                 Rcpp::List t_controlList)
{
  // arma::umat locations = t_SPmatR["locations"];
  // arma::vec values = t_SPmatR["values"];
  // std::cout << "Setting Sparse GRM...." << std::endl;
  // arma::sp_mat SparseGRM = arma::sp_mat(locations, values);
  
  // std::cout << "test, t_flagSparseGRM" << t_flagSparseGRM << std::endl;
  
  // The following function is in POLMM.cpp
  ptr_gPOLMMobj = new POLMM::POLMMClass(t_flagSparseGRM,       // if 1, then use Sparse GRM, otherwise, use Dense GRM
                                        ptr_gDenseGRMobj,
                                        ptr_gPLINKobj,
                                        t_Cova,
                                        t_yVec,     // should be from 0 to J-1
                                        t_beta,
                                        t_bVec,
                                        t_eps,           // 
                                        t_tau,
                                        g_SparseGRM,    // results of function setSparseGRMInCPP()
                                        t_controlList);
  
  ptr_gPOLMMobj->fitPOLMM();
  
  Rcpp::List outList = ptr_gPOLMMobj->getPOLMM();
  
  // ptr_gDenseGRMobj->closeDenseGRMObj();
  return outList;
}


// [[Rcpp::export]]
void setSPACoxobjInCPP(arma::mat t_cumul,
                       arma::vec t_mresid,
                       arma::mat t_XinvXX,
                       arma::mat t_tX,
                       int t_N,
                       double t_pVal_covaAdj_Cutoff,
                       double t_SPA_Cutoff)
{
  ptr_gSPACoxobj = new SPACox::SPACoxClass(t_cumul,
                                           t_mresid,
                                           t_XinvXX,
                                           t_tX,
                                           t_N,
                                           t_pVal_covaAdj_Cutoff,
                                           t_SPA_Cutoff);
}

