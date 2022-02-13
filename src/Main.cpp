// This includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <vector>
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
#include "SAIGE_test.hpp"
#include "UTIL.hpp"

#include <Rcpp.h>
#include "getMem.hpp"

#include <boost/math/distributions/beta.hpp>

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
static VCF::VcfClass* ptr_gVCFobj = NULL;
// global objects for different analysis methods
static SAIGE::SAIGEClass* ptr_gSAIGEobj = NULL;
//single, SAIGE
//Region, SAIGE-GENE+


// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop", //drop is not allowed
double g_missingRate_cutoff;
unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_minMAC_cutoff;    // for Rare Variants (RVs) whose MAC < this value, we aggregate these variants like SAIGE-GENE+ 
double g_marker_minINFO_cutoff;
arma::vec g_region_maxMAF_cutoff;
double g_maxMAFLimit;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage
bool g_isOutputMoreDetails;
int g_marker_chunksize;

std::string g_method_to_CollapseUltraRare;
double g_DosageCutoff_for_UltraRarePresence;

double g_dosage_zerod_MAC_cutoff;
double g_dosage_zerod_cutoff;
bool g_markerTestEnd = false;
arma::vec g_weights_beta(2);

bool  g_is_Firth_beta;
double g_pCutoffforFirth;

// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
			       double t_min_info_marker,
                               unsigned int t_omp_num_threads,
			       bool t_isOutputMoreDetails,
			       int t_marker_chunksize,
			       double t_dosage_zerod_cutoff,
			       double t_dosage_zerod_MAC_cutoff,
			       arma::vec & t_weights_beta
			       )

{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_marker_minINFO_cutoff = t_min_info_marker;
  g_omp_num_threads = t_omp_num_threads;
  g_isOutputMoreDetails = t_isOutputMoreDetails;
  g_marker_chunksize = t_marker_chunksize;
  g_dosage_zerod_cutoff = t_dosage_zerod_cutoff;
  g_dosage_zerod_MAC_cutoff = t_dosage_zerod_MAC_cutoff;
  g_weights_beta = t_weights_beta;
}


//double t_DosageCutoff_for_UltraRarePresence,
			       //std::string t_method_to_CollapseUltraRare,

  //g_method_to_CollapseUltraRare = t_method_to_CollapseUltraRare;
  //g_DosageCutoff_for_UltraRarePresence = t_DosageCutoff_for_UltraRarePresence;
// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               arma::vec t_max_maf_region,
                               unsigned int t_max_markers_region,
                               unsigned int t_omp_num_threads,
			       double t_MACCutoff_to_CollapseUltraRare,
			       double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
			       arma::vec & t_weights_beta)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_maxMAFLimit = g_region_maxMAF_cutoff.max();
  g_region_maxMarkers_cutoff = t_max_markers_region;
  g_omp_num_threads = t_omp_num_threads;
  g_region_minMAC_cutoff = t_MACCutoff_to_CollapseUltraRare;
  g_dosage_zerod_cutoff = t_dosage_zerod_cutoff;
  g_dosage_zerod_MAC_cutoff = t_dosage_zerod_MAC_cutoff;
  g_weights_beta = t_weights_beta;
}



//////// ---------- Main function for marker-level analysis --------- ////////////
                           //std::vector<uint32_t> & t_genoIndex,

// [[Rcpp::export]]
Rcpp::DataFrame mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
			   std::string & t_traitType,
			   std::vector<std::string> & t_genoIndex,
			   bool & t_isMoreOutput,
			   bool & t_isImputation)
{

  int q = t_genoIndex.size();  // number of markers
  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs

  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q);  
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);       
  std::vector<double> pvalVec(q, arma::datum::nan);
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  std::vector<double> pvalNAVec(q, arma::datum::nan);

  bool isCondition = ptr_gSAIGEobj->m_isCondition;
  //if(isCondition){
  std::vector<double> Beta_cVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBeta_cVec(q, arma::datum::nan);
  std::vector<double> pval_cVec(q, arma::datum::nan);
  std::vector<double> Tstat_cVec(q, arma::datum::nan);
  std::vector<double> varT_cVec(q, arma::datum::nan);
  std::vector<double> pvalNA_cVec(q, arma::datum::nan);
  //}
  arma::rowvec G1tilde_P_G2tilde_Vec(ptr_gSAIGEobj->m_numMarker_cond);

  std::vector<bool>  isSPAConvergeVec(q);
  std::vector<double>  AF_caseVec(q);
  std::vector<double>  AF_ctrlVec(q);
  std::vector<uint32_t>  N_caseVec(q);
  std::vector<uint32_t>  N_ctrlVec(q);
    //if(t_isMoreOutput){
  std::vector<double>  N_case_homVec(q);
  std::vector<double>  N_ctrl_hetVec(q);
  std::vector<double>  N_case_hetVec(q);
  std::vector<double>  N_ctrl_homVec(q);
  std::vector<uint32_t>  N_Vec(q);
  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;
  std::vector<uint> indexForMissing;

  int n = ptr_gSAIGEobj->m_n;
  //std::vector<double> t_GVec0;
  //std::vector<double> t_GVec0(n);
  arma::vec t_GVec(n);
  arma::vec gtildeVec(n);
  arma::vec t_P2Vec;
//  }	  
  //ptr_gSAIGEobj->assignSingleVarianceRatio();
  bool hasVarRatio = true;;
  bool isSingleVarianceRatio = true;
  if((ptr_gSAIGEobj->m_varRatio).n_elem == 1){
        ptr_gSAIGEobj->assignSingleVarianceRatio();
  }else{		
	isSingleVarianceRatio = false;
  }

  for(int i = 0; i < q; i++){
    if((i+1) % g_marker_chunksize == 0){
      std::cout << "Completed " << (i+1) << "/" << q << " markers in the chunk." << std::endl;
    }
    // information of marker
    double altFreq, altCounts, missingRate, imputeInfo, AF_case, AF_ctrl, N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom; 
    std::string chr, ref, alt, marker;
    uint32_t pd, N_case, N_ctrl, N;
    bool flip = false;
    std::string t_genoIndex_str = t_genoIndex.at(i);

    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );
    //free(end);
    std::remove(end);

    //Main.cpp
    //PLINK or BGEN 
    //uint32_t gIndex_temp = gIndex; 
    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false; 
   
   //clear vectors
   indexZeroVec.clear();
   indexNonZeroVec.clear();
   indexForMissing.clear();
   //t_GVec0.clear();
   //t_GVec.clear();
   //arma::vec timeoutput1 = getTime();
   bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, t_GVec, t_isImputation);
   //arma::vec timeoutput2 = getTime();   
//printTime(timeoutput1, timeoutput2, "Unified_getOneMarker"); 
//
   if(!isReadMarker){
      //std::cout << "isReadMarker " << isReadMarker << std::endl;
      g_markerTestEnd = true;
      bool isEndFile = check_Vcf_end();
      break;
    }


   //std::cout << "t_GVec0.size()) " << t_GVec0.size() << std::endl;
   //arma::vec t_GVec(t_GVec0.size());
   //arma::vec t_GVec = arma::conv_to< arma::colvec >::from(t_GVec0);

   //arma::vec t_GVec(t_GVec0);
   //t_GVec0.clear(); 

   //for(uint j = 0; j < n; j++){
   //	t_GVec(j) = t_GVec0.at(j);	
   //}

    //for(int indi = 0; indi < indexForNonZero.size(); indi++){
    //  std::cout << indexForNonZero[indi] << std::endl;
    //}
//   std::cout << "marker " << marker << std::endl;
//   std::cout << "indexForMissing.size() " << indexForMissing.size() << std::endl;
//   std::cout << "indexNonZeroVec.size() " << indexNonZeroVec.size() << std::endl;
    //int n = t_GVec.size();
    //arma::vec gtildeVec(n);

  


   std::string pds = std::to_string(pd); 
    std::string info = chr+":"+pds+":"+ref+":"+alt;

    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    refVec.at(i) = ref;
    altVec.at(i) = alt; 
    // record basic information for the marker
    markerVec.at(i) = marker;               // marker IDs
    infoVec.at(i) = info;    // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    //altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    imputationInfoVec.at(i) = imputeInfo;



    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * n * (1 - missingRate) *2;

    
    
   /* 
    
    std::cout << "missingRate " << missingRate << std::endl;
   std::cout << "MAF " << MAF << std::endl;
   std::cout << "MAC " << MAC << std::endl;
   std::cout << "altFreq " << altFreq << std::endl;
   std::cout << "n " << n << std::endl;
   */ 


    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
      continue;
    }else{
    // Check UTIL.cpp
    //
    //
    arma::vec timeoutput3 = getTime();
    indexZeroVec.clear();
    indexNonZeroVec.clear();


    flip = imputeGenoAndFlip(t_GVec, altFreq, altCounts,indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
   

//arma::vec timeoutput4 = getTime();
//printTime(timeoutput3, timeoutput4, "imputeGenoAndFlip");


    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.

   //std::cout << "MAC " << MAC << std::endl; 
   //std::cout << "info " << info << std::endl; 
    // analysis results for single-marker
    double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy;
    double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;

    bool isSPAConverge, is_gtilde;
    //arma::vec t_P2Vec;
    //arma::vec t_P2Vec;

    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
    indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
    indexZeroVec.clear();
    indexNonZeroVec.clear();
    t_P2Vec.clear();
    G1tilde_P_G2tilde_Vec.clear();    
   arma::vec timeoutput5 = getTime(); 
 
 
   if(!isSingleVarianceRatio){ 
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC);
        if(!hasVarRatio){
                //std::cout << "WARNING! Marker " << info << " has MAC " << MAC << " and does not have variance ratio estimated, so the first variance ratio in the variance ratio file is used." << std::endl;
		//std::cout << ptr_gSAIGEobj->m_varRatio.front() << std::endl;
                ptr_gSAIGEobj->assignSingleVarianceRatio_withinput(ptr_gSAIGEobj->m_varRatio.front());
        //        exit(EXIT_FAILURE);
        }
   } 
   
    //check 'Main.cpp'
    bool is_region = false; 
    Unified_getMarkerPval( 
		    t_GVec, 
                          false, // bool t_isOnlyOutputNonZero, 
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,   
			  altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec);

//arma::vec timeoutput6 = getTime();
//printTime(timeoutput5, timeoutput6, "Unified_getMarkerPval");

   indexNonZeroVec_arma.clear();
   indexZeroVec_arma.clear();
   //std::cout << "isSPAConverge " << isSPAConverge << std::endl;
    BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    seBetaVec.at(i) = seBeta;       
    pvalVec.at(i) = pval;
    pvalNAVec.at(i) = pval_noSPA;
    TstatVec.at(i) = Tstat * (1 - 2*flip);
    varTVec.at(i) = varT;

    if(isCondition){
    	Beta_cVec.at(i) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true
    	seBeta_cVec.at(i) = seBeta_c;
    	pval_cVec.at(i) = pval_c;
    	pvalNA_cVec.at(i) = pval_noSPA_c;
    	Tstat_cVec.at(i) = Tstat_c * (1 - 2*flip);
    	varT_cVec.at(i) = varT_c;
    }
	
    if(t_traitType == "binary"){ 
	    arma::vec dosage_case = t_GVec.elem(ptr_gSAIGEobj->m_case_indices);
	    arma::vec dosage_ctrl = t_GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
      AF_case = arma::mean(dosage_case) /2;
      AF_ctrl = arma::mean(dosage_ctrl) /2;
      N_case = dosage_case.n_elem;
      N_ctrl = dosage_ctrl.n_elem;
      if(flip){
         AF_case = 1-AF_case;
         AF_ctrl = 1-AF_ctrl;
      }
      isSPAConvergeVec.at(i) = isSPAConverge;
      AF_caseVec.at(i) = AF_case;
      AF_ctrlVec.at(i) = AF_ctrl;

      N_caseVec.at(i) = N_case;
      N_ctrlVec.at(i) = N_ctrl;

      arma::uvec N_case_ctrl_het_hom0;
      if(t_isMoreOutput){	
   	N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5); 
    	N_case_homVec.at(i)  = N_case_ctrl_het_hom0.n_elem;
    	N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
    	N_case_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
    	N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
    	N_ctrl_homVec.at(i) = N_case_ctrl_het_hom0.n_elem;
    	N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
    	N_ctrl_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
	if(flip){
		N_case_homVec.at(i) = N_case - N_case_hetVec.at(i) -  N_case_homVec.at(i);
		N_ctrl_homVec.at(i) = N_ctrl - N_ctrl_hetVec.at(i) - N_ctrl_homVec.at(i);
	}		
      }	
    }else if(t_traitType == "quantitative"){
      N_Vec.at(i) = n;

    }

    
   } //    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
 
  
  
    //t_GVec.clear();
  }

  //Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
  Rcpp::DataFrame OUT_DF = Rcpp::DataFrame::create(
  //Rcpp::List OutList = Rcpp::List::create(
  	  Rcpp::Named("CHR") = chrVec,
	  Rcpp::Named("POS") = posVec,
	  Rcpp::Named("MarkerID") = markerVec, 
	  Rcpp::Named("Allele1") = refVec, 
	  Rcpp::Named("Allele2") = altVec,
	  Rcpp::Named("AC_Allele2") = altCountsVec,
	  Rcpp::Named("AF_Allele2") = altFreqVec);
	
	 if(t_isImputation){
		OUT_DF["imputationInfo"] = imputationInfoVec;
	 }else{
		OUT_DF["MissingRate"] = missingRateVec;
	 }

	OUT_DF["BETA"] = BetaVec;
	OUT_DF["SE"] = seBetaVec;
	OUT_DF["Tstat"] = TstatVec;
	OUT_DF["var"] = varTVec;
	OUT_DF["p.value"] = pvalVec;

	if(t_traitType == "binary"){
		OUT_DF["p.value.NA"] = pvalNAVec;
		OUT_DF["Is.SPA.converge"] = isSPAConvergeVec;
	    if(isCondition){
		OUT_DF["BETA_c"] = Beta_cVec;
		OUT_DF["SE_c"] = seBeta_cVec;
		OUT_DF["Tstat_c"] = Tstat_cVec;
		OUT_DF["var_c"] = varT_cVec;
		OUT_DF["p.value_c"] = pval_cVec;
		OUT_DF["p.value.NA_c"] = pvalNA_cVec;
	     }
	     OUT_DF["AF_case"] = AF_caseVec;
	     OUT_DF["AF_ctrl"] = AF_ctrlVec;
	     OUT_DF["N_case"] = N_caseVec;
	     OUT_DF["N_ctrl"] = N_ctrlVec;
	     
	     if(t_isMoreOutput){
		OUT_DF["N_case_hom"] = N_case_homVec;
		OUT_DF["N_case_het"] = N_case_hetVec;
		OUT_DF["N_ctrl_hom"] = N_ctrl_homVec;
		OUT_DF["N_ctrl_het"] = N_ctrl_hetVec;
	     }
	}else if(t_traitType == "quantitative"){	
	    if(isCondition){
		OUT_DF["BETA_c"] = Beta_cVec;
		OUT_DF["SE_c"] = seBeta_cVec;
		OUT_DF["Tstat_c"] = Tstat_cVec;
		OUT_DF["var_c"] = varT_cVec;
		OUT_DF["p.value_c"] = pval_cVec;
	    }
	     OUT_DF["N"] = N_Vec;
	}	
	return(OUT_DF);
}




// a unified function to get single marker from genotype file
bool Unified_getOneMarker(std::string & t_genoType,   // "PLINK", "BGEN", "Vcf"
                               uint64_t & t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool & t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint>& t_indexForMissing,     // index of missing genotype data
                               bool & t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint>& t_indexForNonZero, //
			       arma::vec & t_GVec,
			       bool t_isImputation 
			       )     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
  //arma::vec GVec(ptr_gSAIGEobj->m_n);
  bool isBoolRead = true;
  if(t_genoType == "plink"){
   bool isTrueGenotype = true;
   ptr_gPLINKobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                       t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                       isTrueGenotype, t_GVec);   // t_isTrueGenotype, only used for PLINK format.
  }
  
  if(t_genoType == "bgen"){
    //bool isBoolRead = true;
    ptr_gBGENobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead, t_GVec, t_isImputation);
  }

  if(t_genoType == "vcf"){
    ptr_gVCFobj->getOneMarker(t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero, isBoolRead, t_GVec, t_isImputation);
    ptr_gVCFobj->move_forward_iterator(1);
  }	  
  
  return isBoolRead;
}

// a unified function to get marker-level p-value
void Unified_getMarkerPval(
                           arma::vec & t_GVec,
                           bool t_isOnlyOutputNonZero,
			   arma::uvec & t_indexForNonZero_vec,
			   arma::uvec & t_indexForZero_vec,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
			   double& t_pval_noSPA,
                           double& t_Tstat,
			   double& t_gy,
                           double& t_varT,
                           double t_altFreq, 
			   bool & t_isSPAConverge,
			   arma::vec & t_gtilde,
			   bool & is_gtilde,
			   bool  is_region,
                           arma::vec & t_P2Vec,
			   bool t_isCondition,
                           double& t_Beta_c, 
                           double& t_seBeta_c, 
                           double& t_pval_c,
			   double& t_pval_noSPA_c,
                           double& t_Tstat_c,
                           double& t_varT_c,
			   arma::rowvec & t_G1tilde_P_G2tilde_Vec)
{
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAIGE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' should be false.");   

    ptr_gSAIGEobj->getMarkerPval(t_GVec, t_indexForNonZero_vec, t_indexForZero_vec, t_Beta, t_seBeta, t_pval, t_pval_noSPA, t_altFreq, t_Tstat, t_gy, t_varT, t_isSPAConverge, t_gtilde, is_gtilde, is_region, t_P2Vec, t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c, t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec); //SAIGE_new.cpp
    
    //t_indexForNonZero_vec.clear();
  
}

//////// ---------- Main functions to set objects for different genotype format --------- ////////////

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> & t_SampleInModel,
                      std::string t_AlleleOrder)
{
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                        t_famFile,
                                        t_bedFile,
                                        t_SampleInModel,
                                        t_AlleleOrder);
  
  
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> & t_SampleInBgen,
                     std::vector<std::string> & t_SampleInModel,
                     std::string t_AlleleOrder)
{
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                     t_bgenFileIndex,
                                     t_SampleInBgen,
                                     t_SampleInModel,
				     false,
				     false,
                                     t_AlleleOrder);
  //int n = ptr_gBGENobj->getN();
}


// [[Rcpp::export]]
void setVCFobjInCPP(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::vector<std::string> & t_SampleInModel)
{
  ptr_gVCFobj = new VCF::VcfClass(t_vcfFileName,
		  		t_vcfFileIndex,
				t_vcfField,
				false,
				t_SampleInModel);
		  
}



//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

// [[Rcpp::export]]
void setSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
        arma::mat & t_X,
        arma::vec &  t_S_a,
        arma::vec & t_res,
        arma::vec & t_mu2,
        arma::vec & t_mu,
        arma::vec & t_varRatio,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::vec & t_tauvec,
        std::string t_traitType,
        arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
        arma::umat & t_locationMat,
        arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
	std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
	double t_pCutoffforFirth,
	arma::vec & t_offset)
{
  // check SAIGE.cpp
  ptr_gSAIGEobj = new SAIGE::SAIGEClass(
	t_XVX,
        t_XXVX_inv,
        t_XV,
        t_XVX_inv_XV,
        t_X,
        t_S_a,
        t_res,
        t_mu2,
        t_mu,
        t_varRatio,
	t_cateVarRatioMinMACVecExclude,
	t_cateVarRatioMaxMACVecInclude,
        t_SPA_Cutoff,
        t_tauvec,
        t_traitType,
        t_y,
	t_impute_method,
	t_flagSparseGRM,
	t_locationMat,
	t_valueVec,
	t_dimNum, 
	t_isCondition,
	t_condition_genoIndex,
	t_is_Firth_beta,
        t_pCutoffforFirth,
	t_offset);
  //ptr_gSAIGEobj->m_flagSparseGRM = false;
		 	  
}



// [[Rcpp::export]]
void setSparseSigmaInCPP(int r, arma::umat & t_locationMatinR, arma::vec & t_valueVecinR)
{
  ptr_gSAIGEobj->setupSparseMat(r, t_locationMatinR, t_valueVecinR);
  ptr_gSAIGEobj->m_flagSparseGRM = true;
}


// [[Rcpp::export]]
Rcpp::List RegionSetUpConditional_binary_InCPP(arma::vec & t_weight_cond){

	unsigned int q_cond = (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows;
  	boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);
  	arma::vec w0G2Vec_cond(q_cond);
  	double w0G2_cond, MAFG2_cond;
        for(unsigned int ci = 0; ci < q_cond; ci++){
		if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		}else{
                	MAFG2_cond = (ptr_gSAIGEobj->m_MAF_cond)[ci];
                	w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		}
                w0G2Vec_cond.at(ci) = w0G2_cond;
        }
	arma::mat m_VarMat_weighted_cond = (w0G2Vec_cond * (w0G2Vec_cond.t())) % (ptr_gSAIGEobj->m_VarMat_cond);

	 Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("VarMat_G2_cond") = m_VarMat_weighted_cond,
                                          Rcpp::Named("Score_G2_cond") = ptr_gSAIGEobj->m_Tstat_cond,
                                          Rcpp::Named("pval_G2_cond") = ptr_gSAIGEobj->m_p_cond,
                                          Rcpp::Named("gsum_G2_cond") = ptr_gSAIGEobj->m_gsum_cond,
                                          Rcpp::Named("qsum_G2_cond") = ptr_gSAIGEobj->m_qsum_cond
					  );
	 return(OutList);

}


//////// ---------- Main function for region-level analysis --------- ////////////
// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_genoIndex,
			   arma::mat & annoIndicatorMat,
			   arma::vec & maxMAFVec, 
                           std::string t_outputFile,
			   std::string t_traitType,
                           unsigned int t_n,           // sample size  
                           arma::mat P1Mat,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
                           arma::mat P2Mat, 
			   std::string t_regionTestType, 
			   bool t_isImputation,
			   arma::vec & t_weight,
			   arma::vec & t_weight_cond)
{

  bool isWeightCustomized = false;
	

  //std::cout << "okk1" << std::endl;
  unsigned int q0 = t_genoIndex.size();                 // number of markers (before QC) in one region
  if(!(t_weight.is_zero()) && t_weight.n_elem == q0){
     isWeightCustomized = true;	
  } 	  
  unsigned int q_anno = annoIndicatorMat.n_cols;
  unsigned int q_maf = maxMAFVec.n_elem;
  unsigned int q_anno_maf = q_anno*q_maf;
  arma::mat genoURMat(t_n, q_anno_maf, arma::fill::zeros);

  arma::vec weightURVec(q_anno_maf, arma::fill::zeros);


  unsigned int q = q0 + q_anno_maf;
  arma::imat annoMAFIndicatorMat(q, q_anno_maf, arma::fill::zeros);
  //std::cout << "okk2" << std::endl;
  arma::ivec annoMAFIndicatorVec(q_anno_maf);
  annoMAFIndicatorVec.zeros();
	
  unsigned int q_cond = (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows;	
  arma::rowvec G1tilde_P_G2tilde_Vec(q_cond);
  //std::cout << "okk3" << std::endl;
  boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);

  bool isCondition = ptr_gSAIGEobj->m_isCondition;
  arma::vec w0G2Vec_cond(q_cond);
  double w0G2_cond, MAFG2_cond;
  if(isCondition){
	for(unsigned int ci = 0; ci < q_cond; ci++){
		if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		}else{	
			MAFG2_cond = (ptr_gSAIGEobj->m_MAF_cond)[ci];	
  			w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		}	
		w0G2Vec_cond.at(ci) = w0G2_cond;
	}
  }
  arma::mat w0G2Mat_cond(q_cond, q_cond);
  w0G2Mat_cond = w0G2Vec_cond * (w0G2Vec_cond.t());
  arma::mat genoSumMat(t_n, q_anno_maf, arma::fill::zeros); //for Phi_cc for binary traits

  //if(isCondition){
	std::vector<double> Beta_cVec(q, arma::datum::nan);         // beta value for ALT allele
        std::vector<double> seBeta_cVec(q, arma::datum::nan);
        std::vector<double> pval_cVec(q, arma::datum::nan);
        std::vector<double> Tstat_cVec(q, arma::datum::nan);
        std::vector<double> varT_cVec(q, arma::datum::nan);
        std::vector<double> pvalNA_cVec(q, arma::datum::nan);
  	arma::mat G1tilde_P_G2tilde_Weighted_Mat(q, q_cond);
  //}
  //group test output
  //arma::vec MAC_GroupVec = arma::zeros<vec>(q_anno_maf);
  arma::vec MAC_GroupVec(q_anno_maf);
  MAC_GroupVec.zeros();
  arma::vec MACCase_GroupVec(q_anno_maf);
  MACCase_GroupVec.zeros();
  arma::vec MACControl_GroupVec(q_anno_maf);
  MACControl_GroupVec.zeros();
  arma::vec NumRare_GroupVec(q_anno_maf);
  NumRare_GroupVec.zeros();
  arma::vec NumUltraRare_GroupVec(q_anno_maf);
  NumUltraRare_GroupVec.zeros();
  //arma::uvec g_case_indices;
  //arma::uvec g_ctrl_indices;
  arma::vec gtildeVec;
  double MACgroup, MACcasegroup, MACcontrolgroup, AF_case, AF_ctrl;
  //if(t_traitType == "binary"){
  //  ptr_gSAIGEobj->getindices(g_case_indices, g_ctrl_indices);
  //}

  // added on 09-18-2021
  arma::uvec indicatorVec(q, arma::fill::zeros);       // 0: does not pass QC, 1: non-URV, 2: URV
  std::vector<std::string> markerVec(q);

  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs

  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MACVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MAFVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> N_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> N_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> N_Vec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q, arma::datum::nan);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q, arma::datum::nan);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q, arma::datum::nan);
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);
  std::vector<double> pvalVec(q, arma::datum::nan);
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> TstatVec_flip(q, arma::datum::nan);
  std::vector<double> gyVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  std::vector<double> pvalNAVec(q, arma::datum::nan);  
  std::vector<bool>  isSPAConvergeVec(q);


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
  double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy;
  double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;
  bool isSPAConverge, is_gtilde;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  //std::vector<double> GVec0(t_n);
  arma::vec GVec(t_n);
  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;


  bool hasVarRatio = true;;
  bool isSingleVarianceRatio = true;
  if((ptr_gSAIGEobj->m_varRatio).n_elem == 1){
        ptr_gSAIGEobj->assignSingleVarianceRatio();
  }else{
        isSingleVarianceRatio = false;
  }



  // initiate chunk information
  unsigned int nchunks = 0; //number of chunks
  unsigned int ichunk = 0; //ith chunk
  unsigned int i1InChunk = 0; //i1th marker in ith chunk
  unsigned int i1 = 0;    // index of Markers (non-URV)
  unsigned int i2 = 0;    // index of Markers (Ultra-Rare Variants, URV)

  unsigned int jm; 
  //std::cout << "q0 " << q0 << std::endl; 
  // cycle for q markers
  for(unsigned int i = 0; i < q0; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
   
    std::string t_genoIndex_str = t_genoIndex.at(i);

    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );
    std::remove(end);

    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false;
    //  std::vector<uint> indexZeroVec;
  //std::vector<uint> indexNonZeroVec;
    bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec,
					  GVec,
					  t_isImputation);
    if(!isReadMarker){
      break;
    }	    
   //arma::vec GVec(GVec0);
   //GVec0.clear();
    std::string pds = std::to_string(pd);
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;


    double MAF = std::min(altFreq, 1 - altFreq);
    double w0;
    if(isWeightCustomized){
	w0 = t_weight(i);
    }else{
	w0 = boost::math::pdf(beta_dist, MAF);
    }
    //double w0 = boost::math::pdf(beta_dist, MAF);
    double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
    flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);

       arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
       indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
       indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);


    MAF = std::min(altFreq, 1 - altFreq);
    MAC = std::min(altCounts, t_n *2 - altCounts);
    // Quality Control (QC)
    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    refVec.at(i) = ref;
    altVec.at(i) = alt;


    markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    altCountsVec.at(i) = altCounts;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    imputationInfoVec.at(i) = imputeInfo;
    

    if((missingRate > g_missingRate_cutoff) || (MAF > g_maxMAFLimit) || MAF == 0 || (imputeInfo < g_marker_minINFO_cutoff)){
      continue;  // does not pass QC
    }




    /*
    markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    altCountsVec.at(i) = altCounts;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    imputationInfoVec.at(i) = imputeInfo;
    */
    //ptr_gSAIGEobj->assignVarianceRatio(MAC);




    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants
    
      if(!isSingleVarianceRatio){	    
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC);
        if(!hasVarRatio){
                std::cout << "Error! Marker " << info << " has MAC " << MAC << " and does not have variance ratio estimated." << std::endl;
                exit(EXIT_FAILURE);
        }
      }
      
      indicatorVec.at(i) = 1;
      
      if(i1InChunk == 0){
        std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
      }
     
       
      Unified_getMarkerPval(
                    GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec);
      BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true       
      seBetaVec.at(i) = seBeta;       
      pvalVec.at(i) = pval;
      pvalNAVec.at(i) = pval_noSPA;
      TstatVec.at(i) = Tstat * (1 - 2*flip);
      //gyVec.at(i) = gy * (1-2*flip);
      TstatVec_flip.at(i) = Tstat;
      gyVec.at(i) = gy;
      varTVec.at(i) = varT;
      isSPAConvergeVec.at(i) = isSPAConverge;

      if(isCondition){ 	
      	Beta_cVec.at(i) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
      	seBeta_cVec.at(i) = seBeta_c;
      	pval_cVec.at(i) = pval_c;
      	pvalNA_cVec.at(i) = pval_noSPA_c;
      	Tstat_cVec.at(i) = Tstat_c * (1 - 2*flip);
      	varT_cVec.at(i) = varT_c;
	G1tilde_P_G2tilde_Weighted_Mat.row(i) = G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * w0;	
	//G1tilde_P_G2tilde_Mat.row(i) = G1tilde_P_G2tilde_Vec;	
      }

      int n = GVec.size();
      if(t_regionTestType != "BURDEN"){ 
        P1Mat.row(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
        P2Mat.col(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
      }
      //P1Mat.row(i1InChunk) = gtildeVec.t();
      //P2Mat.col(i1InChunk) = P2Vec;
    
      i1 += 1;
      i1InChunk += 1;

      //std::cout << "i1InChunk is " << i1InChunk << std::endl;
      //std::cout << "ichunk is " << ichunk << std::endl;
     //std::cout << "i is " << i << std::endl;



      arma::vec dosage_case, dosage_ctrl;
      if(t_traitType == "binary"){
                        dosage_case = GVec.elem(ptr_gSAIGEobj->m_case_indices);
                        dosage_ctrl = GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);	
      }

      arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
      MAFIndicatorVec.zeros();
      MAFIndicatorVec.elem( find(maxMAFVec >= MAF) ).ones();	

      annoMAFIndicatorVec.zeros();
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
		for(unsigned int m = 0; m < q_maf; m++){
			if(MAFIndicatorVec(m) == 1){
				jm = j*q_maf + m;	
				annoMAFIndicatorVec(jm) = 1;
				MAC_GroupVec(jm) = MAC_GroupVec(jm) + MAC;
				if(t_traitType == "binary"){
					MACCase_GroupVec(jm) = MACCase_GroupVec(jm) + MACcasegroup;
					MACControl_GroupVec(jm) = MACControl_GroupVec(jm) + MACcontrolgroup;
  //					double w0 = boost::math::pdf(beta_dist, MAF);
					genoSumMat.col(jm) = genoSumMat.col(jm) + w0*GVec;


				}
				NumRare_GroupVec(jm) = NumRare_GroupVec(jm) + 1;
			}

		}	
        }
      }
     annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
      if(t_traitType == "binary"){
      AF_case = arma::mean(dosage_case) /2;
      AF_ctrl = arma::mean(dosage_ctrl) /2;
      if(flip){
         AF_case = 1-AF_case;
         AF_ctrl = 1-AF_ctrl;
      }
      AF_caseVec.at(i) = AF_case;
      AF_ctrlVec.at(i) = AF_ctrl;
      N_caseVec.at(i) = dosage_case.n_elem;
      N_ctrlVec.at(i) = dosage_ctrl.n_elem;
     }else if(t_traitType == "quantitative"){
      N_Vec.at(i) = n;
    }      

      
    }else{   // Ultra-Rare Variants (URV)
      //std::cout << "MAC: " << MAC << std::endl;	    
      indicatorVec.at(i) = 2;
      arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
      MAFIndicatorVec.zeros();
      MAFIndicatorVec.elem( find(maxMAFVec >= MAF) ).ones();
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
		for(unsigned int m = 0; m < q_maf; m++){
                        if(MAFIndicatorVec(m) == 1){
                        	jm = j*q_maf + m;
				genoURMat.col(jm) = arma::max(genoURMat.col(jm), GVec);

				if(isWeightCustomized){
					weightURVec(jm) = std::max(weightURVec(jm), t_weight(i));  
				}	
				//arma::vec genoURMatcol_jm = genoURMat.col(jm);
				//arma::uvec genoURMatcol_jm_nonzero = arma::find(genoURMatcol_jm != 0);
				//arma::vec genoURMatcol_jm_sub = genoURMatcol_jm(genoURMatcol_jm_nonzero);
				//genoURMatcol_jm_sub.print();
				NumUltraRare_GroupVec(jm) = NumUltraRare_GroupVec(jm) + 1;
			}
		}
	}
      }	
      i2 += 1;
    }
  

    
    if(i1InChunk == m1 ){
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      if(t_regionTestType != "BURDEN"){
        P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      }

      mPassCVVec.push_back(m1);
      ichunk += 1;
      i1InChunk = 0;
      nchunks = nchunks + 1;
    }
/*
    std::cout << "i1InChunk is " << i1InChunk << std::endl;
    std::cout << "ichunk is " << ichunk << std::endl;
    std::cout << "m1 is " << m1 << std::endl;
    std::cout << "P1Mat.n_rows is " << P1Mat.n_rows << std::endl;
*/
    Rcpp::checkUserInterrupt();
  }
  /*
  for(unsigned int maxi = 0; maxi < MACVec.size(); maxi++){  
    std::cout << "MACVec.at(mi) " << MACVec.at(maxi) << std::endl;
  }
*/
/*
  if(i1 == 0){
    std::cout << "Only ultra-rare variants are found. This region will be skipped." << std::endl;
    Rcpp::List OutList = Rcpp::List::create();
    return OutList;
  }
*/

/*
  std::cout << "ok1" << std::endl;
  std::cout << "i1InChunk is " << i1InChunk << std::endl;
  std::cout << "ichunk is " << ichunk << std::endl;
  std::cout << "m1 is " << m1 << std::endl;
  std::cout << "ok2" << std::endl;
*/
//the second last chunk
  if(i1InChunk != 0){

    std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
   if(t_regionTestType != "BURDEN"){	  
    P1Mat = P1Mat.rows(0, i1InChunk - 1);
    P2Mat = P2Mat.cols(0, i1InChunk - 1);
    //std::cout << "nchunks " << nchunks << std::endl;
    //if(nchunks != 1){
    P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
   }
    ichunk = ichunk + 1;
    //}
    mPassCVVec.push_back(i1InChunk);
    //nchunks = nchunks + 1;
/*
    std::cout << "i1InChunk 2 is " << i1InChunk << std::endl;
    std::cout << "ichunk 2 is " << ichunk << std::endl;
    std::cout << "m1 2 is " << m1 << std::endl;
    std::cout << "P1Mat.n_rows " << P1Mat.n_rows << std::endl;
    std::cout << "P1Mat.n_cols " << P1Mat.n_cols << std::endl;
*/
    nchunks = nchunks + 1; 
    i1InChunk = 0;
  }
  //else{
  //  ichunk = ichunk - 1;
  //}	  


  //std::vector<unsigned int> URindVec;
//if(i2 != 0){  //if there are ultra-rare variants
//
  //arma::mat P1Mat0 = P1Mat;
  //arma::mat P2Mat0 = P1Mat;

  if(i2 > 0){
  int m1new = std::max(m1, q_anno_maf);
 if(t_regionTestType != "BURDEN"){
  P1Mat.resize(m1new, P1Mat.n_cols);
  P2Mat.resize(P2Mat.n_rows,m1new);
 }  
  arma::mat XV, XXVX_inv;
  ptr_gSAIGEobj->extract_XV_XXVX_inv(XV, XXVX_inv);
  //the last chunk for UR
  unsigned int i;
  for(unsigned int j = 0; j < q_anno; j++){
     for(unsigned int m = 0; m < q_maf; m++){
	jm = j*q_maf+m;
	arma::vec genoURVec = genoURMat.col(jm);
	arma::uvec indexForNonZero = arma::find(genoURVec != 0);
//arma::uvec indexForNonZero = arma::find(genoURVec > 0.2);
//std::cout << "indexForNonZero " << std::endl;
//std::cout << "genoURVecNonZero " << std::endl;
//arma::uvec indexForNonZero_sort = arma::sort(indexForNonZero);
//indexForNonZero_sort.print();
//arma::vec genoURVecNonZero = genoURVec(indexForNonZero);
//genoURVecNonZero.print();
//std::cout << genoURVec[indexForNonZero_sort[0]] << std::endl;
//	std::cout << "jm " << jm << std::endl;
//	std::cout << "indexForNonZero.n_elem" << indexForNonZero.n_elem << std::endl;
	i = q0 + jm;
	markerVec.at(i) = "UR";             // marker IDs
	if(indexForNonZero.n_elem > 0){
	//URindVec.push_back(jm+1);
	double altFreq = arma::mean(genoURVec)/2;
	double altCounts = arma::accu(genoURVec);
	double missingRate = 0;
	double imputeInfo = 1;
    	std::string chr, ref, alt, marker;
    	//uint32_t pd;
    	bool flip = false;
	std::string info = "UR";	
    	double MAF = std::min(altFreq, 1 - altFreq);
	double w0;
	if(isWeightCustomized){
            w0 = weightURVec(jm);		    
    	}else{
            w0 = boost::math::pdf(beta_dist, MAF);
    	}
	//double w0 = boost::math::pdf(beta_dist, MAF);	
	genoSumMat.col(jm) = genoSumMat.col(jm) + genoURVec * w0;
	arma::vec genoSumMatvec1 = genoSumMat.col(jm);
	arma::vec genoSumMatvec2 = XV * genoSumMatvec1;
	arma::vec genoSumMatvec3 = genoSumMatvec1 - XXVX_inv * genoSumMatvec2;
	genoSumMat.col(jm) = genoSumMatvec3;

	double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021

	if(!isSingleVarianceRatio){	
        	hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC);
        	if(!hasVarRatio){
                	//std::cout << "Error! Collapsed ultra rare marker with the annotation " << j+1 << " and MAF <= " << maxMAFVec(m) << " has MAC " << MAC << " and does not have variance ratio estimated." << std::endl;
			//std::cout << "Please fit the null model in Step 1 using a sparse GRM without variance ratio estimation or fit the null model in Step 1 using a full GRM with variance ratio estimated for lower MAC categories." << std::endl;"
			//exit(EXIT_FAILURE);
			hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(g_region_minMAC_cutoff);
		}

		if(!hasVarRatio){
			ptr_gSAIGEobj->assignSingleVarianceRatio_withinput(1.0);	

		}	
  	}



	annoMAFIndicatorVec.zeros();
	annoMAFIndicatorVec(jm) = 1;
	annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
    	infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    	altFreqVec.at(i) = altFreq;	// allele frequencies of ALT allele, this is not always < 0.5.
	altCountsVec.at(i) = altCounts;
	missingRateVec.at(i) = missingRate;
    	MACVec.at(i) = MAC;
    	MAFVec.at(i) = MAF;
	std::vector<uint32_t> indexForMissing;

    	flip = imputeGenoAndFlip(genoURVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);


     arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
       indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
       indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);



	ptr_gSAIGEobj->getMarkerPval(genoURVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec);
	int n = genoURVec.size();

	//Unified_getMarkerPval(
        //            genoURVec,
        //                  false, // bool t_isOnlyOutputNonZero,
        //                  indexForNonZero, Beta, seBeta, pval, pval_noSPA, Tstat, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec);

      BetaVec.at(i) = Beta* (1 - 2*flip);
      seBetaVec.at(i) = seBeta;
      pvalVec.at(i) = pval;
      pvalNAVec.at(i) = pval_noSPA;
      TstatVec.at(i) = Tstat * (1 - 2*flip);
      TstatVec_flip.at(i) = Tstat;
      gyVec.at(i) = gy;
      //gyVec.at(i) = gy * (1 - 2*flip);
      varTVec.at(i) = varT;
      isSPAConvergeVec.at(i) = isSPAConverge;
      if(isCondition){
        Beta_cVec.at(i) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
        seBeta_cVec.at(i) = seBeta_c;
        pval_cVec.at(i) = pval_c;
        pvalNA_cVec.at(i) = pval_noSPA_c;
        Tstat_cVec.at(i) = Tstat_c * (1 - 2*flip);
        varT_cVec.at(i) = varT_c;
        G1tilde_P_G2tilde_Weighted_Mat.row(i) = G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * w0;
        //G1tilde_P_G2tilde_Mat.row(i) = G1tilde_P_G2tilde_Vec;
      }

      arma::vec dosage_case, dosage_ctrl;
      MAC_GroupVec(jm) = MAC_GroupVec(jm) + MAC;
                if(t_traitType == "binary"){
                        dosage_case = genoURVec.elem(ptr_gSAIGEobj->m_case_indices);
                        dosage_ctrl = genoURVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);
                        MACCase_GroupVec(jm) = MACCase_GroupVec(jm) + MACcasegroup;
                        MACControl_GroupVec(jm) = MACControl_GroupVec(jm) + MACcontrolgroup;
                }

      if(t_traitType == "binary"){
      AF_case = arma::mean(dosage_case) /2;
      AF_ctrl = arma::mean(dosage_ctrl) /2;
      if(flip){
         AF_case = 1-AF_case;
         AF_ctrl = 1-AF_ctrl;
      }
      AF_caseVec.at(i) = AF_case;
      AF_ctrlVec.at(i) = AF_ctrl;
      N_caseVec.at(i) = dosage_case.n_elem;
      N_ctrlVec.at(i) = dosage_ctrl.n_elem;
     }else if(t_traitType == "quantitative"){
      N_Vec.at(i) = n;
    }

      // adjPVec.at(i1) = pval1;
    /*std::cout << "here9" << std::endl;
    std::cout << "i1InChunk " << i1InChunk  << std::endl;
    std::cout << P1Mat.n_rows << std::endl;
    std::cout << P1Mat.n_cols << std::endl;
    std::cout << gtildeVec.n_elem << std::endl;
    */
    if(t_regionTestType != "BURDEN"){  
      P1Mat.row(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
      P2Mat.col(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
    }
    //P1Mat.row(i1InChunk) = gtildeVec.t();
    //P2Mat.col(i1InChunk) = P2Vec;
    //std::cout << "here11" << std::endl;
    i1InChunk = i1InChunk + 1;
    i1 = i1 + 1;
    }
	//else{//if(indexForNonZero.n_elem > 0){
	//P1Mat.row(i1InChunk) =		    
    //}	    
   }  
  }


  //std::cout << "i1 = i1 + q_anno;" << std::endl;
  //i1 = i1 + q_anno_maf;
  //if(i1InChunk > 0){
   // }
//} //if(i2 != 0){if there are ultra-rare variants
  // std::cout << "P1Mat.n_rows ok1 " << P1Mat.n_rows << std::endl; 
  if(i1InChunk != 0){
    nchunks = nchunks + 1;
    if(t_regionTestType != "BURDEN"){
      P1Mat = P1Mat.rows(0, i1InChunk - 1);
      P2Mat = P2Mat.cols(0, i1InChunk - 1);
//if(nchunks != 1){
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    }
      ichunk = ichunk + 1;
//    }
    mPassCVVec.push_back(i1InChunk);
  }
   //std::cout << "P1Mat.n_rows ok2 " << P1Mat.n_rows << std::endl; 

  }// if(i2 > 0)    
  int mPassCVVecsize = mPassCVVec.size();
  nchunks = mPassCVVecsize;
  // not so many markers in the region, so all matrix is in memory
  //nchunks = ichunk + 1;
  //


arma::mat VarMat(i1, i1);

if(t_regionTestType != "BURDEN"){
  if(nchunks == 1){
    VarMat = P1Mat * P2Mat;
  //std::cout << "VarMat.n_rows " << VarMat.n_rows << std::endl;
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

	//std::cout << "P1Mat.n_cols " << P1Mat.n_cols << std::endl; 
        //std::cout << "P2Mat.n_cols " << P2Mat.n_cols << std::endl; 
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
      //std::cout << "P2Mat.n_cols " << P2Mat.n_cols << std::endl;
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
}

   Rcpp::DataFrame OUT_DF = Rcpp::DataFrame::create(
  //Rcpp::List OutList = Rcpp::List::create(
          Rcpp::Named("CHR") = chrVec,
          Rcpp::Named("POS") = posVec,
          Rcpp::Named("MarkerID") = markerVec,
          Rcpp::Named("Allele1") = refVec,
          Rcpp::Named("Allele2") = altVec,
          Rcpp::Named("AC_Allele2") = altCountsVec,
          Rcpp::Named("AF_Allele2") = altFreqVec);

    if(t_isImputation){
                OUT_DF["imputationInfo"] = imputationInfoVec;
         }else{
                OUT_DF["MissingRate"] = missingRateVec;
         }

        OUT_DF["BETA"] = BetaVec;
        OUT_DF["SE"] = seBetaVec;
        OUT_DF["Tstat"] = TstatVec;
        OUT_DF["var"] = varTVec;
        OUT_DF["p.value"] = pvalVec;

        if(t_traitType == "binary"){
                OUT_DF["p.value.NA"] = pvalNAVec;
                OUT_DF["Is.SPA.converge"] = isSPAConvergeVec;
            if(isCondition){
                OUT_DF["BETA_c"] = Beta_cVec;
                OUT_DF["SE_c"] = seBeta_cVec;
                OUT_DF["Tstat_c"] = Tstat_cVec;
                OUT_DF["var_c"] = varT_cVec;
                OUT_DF["p.value_c"] = pval_cVec;
                OUT_DF["p.value.NA_c"] = pvalNA_cVec;
             }
             OUT_DF["AF_caseVec"] = AF_caseVec;
             OUT_DF["AF_ctrlVec"] = AF_ctrlVec;
             OUT_DF["N_caseVec"] = N_caseVec;
             OUT_DF["N_ctrlVec"] = N_ctrlVec;

	     /*
             if(t_isMoreOutput){
                OUT_DF["N_case_hom"] = N_case_homVec;
                OUT_DF["N_case_het"] = N_case_hetVec;
                OUT_DF["N_ctrl_hom"] = N_ctrl_homVec;
                OUT_DF["N_ctrl_het"] = N_ctrl_hetVec;
             }
	     */
        }else if(t_traitType == "quantitative"){
            if(isCondition){
                OUT_DF["BETA_c"] = Beta_cVec;
                OUT_DF["SE_c"] = seBeta_cVec;
                OUT_DF["Tstat_c"] = Tstat_cVec;
                OUT_DF["var_c"] = varT_cVec;
                OUT_DF["p.value_c"] = pval_cVec;
            }
             OUT_DF["N"] = N_Vec;
        }
  

  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("OUT_DF") = OUT_DF,
		  			  Rcpp::Named("VarMat") = VarMat,
		  			  Rcpp::Named("annoMAFIndicatorMat") = annoMAFIndicatorMat,
					  Rcpp::Named("MAC_GroupVec") = MAC_GroupVec,
					  Rcpp::Named("MAFVec") = MAFVec,
                                          Rcpp::Named("TstatVec_flip") = TstatVec_flip,
					  Rcpp::Named("NumRare_GroupVec") =NumRare_GroupVec,
					  Rcpp::Named("NumUltraRare_GroupVec") = NumUltraRare_GroupVec
/*
					  Rcpp::Named("markerVec") = markerVec,
              				  Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("altCountsVec") = altCountsVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
					  Rcpp::Named("imputationInfoVec") =imputationInfoVec,
					  Rcpp::Named("isSPAConvergeVec") =isSPAConvergeVec,
                                          Rcpp::Named("pvalVec") = pvalVec,
                                          Rcpp::Named("BetaVec") = BetaVec,
                                          Rcpp::Named("seBetaVec") = seBetaVec,
                                          Rcpp::Named("TstatVec") = TstatVec,
                                          Rcpp::Named("varTVec") = varTVec
					  */
					  );

  if(t_traitType == "binary"){
    OutList.push_back(MACCase_GroupVec, "MACCase_GroupVec");
    OutList.push_back(MACControl_GroupVec, "MACCtrl_GroupVec");
    /*
    OutList.push_back(pvalNAVec, "pvalNAVec");
    OutList.push_back(AF_caseVec, "AF_caseVec");
    OutList.push_back(AF_ctrlVec, "AF_ctrlVec");
    OutList.push_back(N_caseVec, "N_caseVec");
    OutList.push_back(N_ctrlVec, "N_ctrlVec");
    */
    OutList.push_back(genoSumMat, "genoSumMat");
    OutList.push_back(gyVec, "gyVec");
  }
  //else if(t_traitType == "quantitative"){
  //  OutList.push_back(N_Vec, "N_Vec");
  //}

  //arma::mat scaled_m_VarInvMat_cond;
  if(isCondition){
  //std::cout << "okk5" << std::endl;
    arma::mat AdjCondMat = G1tilde_P_G2tilde_Weighted_Mat * (ptr_gSAIGEobj->m_VarInvMat_cond / (w0G2Mat_cond));
    arma::mat VarMatAdjCond = AdjCondMat * (G1tilde_P_G2tilde_Weighted_Mat.t());
    arma::vec TstatAdjCond = AdjCondMat * (ptr_gSAIGEobj->m_Tstat_cond % w0G2Vec_cond ); 
    OutList.push_back(G1tilde_P_G2tilde_Weighted_Mat, "G1tilde_P_G2tilde_Weighted_Mat"); 
    OutList.push_back(ptr_gSAIGEobj->m_scalefactor_G2_cond, "scalefactor_G2_cond");
    OutList.push_back(ptr_gSAIGEobj->m_VarInvMat_cond_scaled_weighted, "VarInvMat_G2_cond_scaled"); 
    OutList.push_back(ptr_gSAIGEobj->m_Tstat_cond, "Tstat_G2_cond"); //m_Tstat_cond is weighted
    OutList.push_back(ptr_gSAIGEobj->m_G2_Weight_cond, "G2_Weight_cond");
    OutList.push_back(TstatAdjCond, "TstatAdjCond");
    OutList.push_back(VarMatAdjCond, "VarMatAdjCond"); 
    //OutList.push_back(Beta_cVec, "Beta_cVec");
    //OutList.push_back(seBeta_cVec, "seBeta_cVec");
    //OutList.push_back(pval_cVec, "pval_cVec");
    //OutList.push_back(Tstat_cVec, "Tstat_cVec");
    //OutList.push_back(varT_cVec, "varT_cVec");
  }  
    //if(t_traitType == "binary"){
    //	OutList.push_back(pvalNA_cVec, "pvalNA_cVec");
    // }
//}
  return OutList;
}



// [[Rcpp::export]]
void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "plink", "bgen", "vcf"
                           std::vector<std::string> & t_genoIndex,
                           unsigned int t_n, 
			   arma::vec & t_weight_cond
			   )           // sample size
{
  bool isImpute = false;	
  unsigned int q = t_genoIndex.size();
  arma::mat P1Mat(q, t_n);
  arma::mat P2Mat(t_n, q);
  arma::mat VarInvMat(q, q);
  arma::vec TstatVec(q);
  arma::vec pVec(q);
  arma::vec MAFVec(q);
  arma::vec gyVec(q);
  arma::vec w0G2_cond_Vec(q);
  arma::vec gsumVec(t_n, arma::fill::zeros);
  //double beta1 = g_weights_beta[0];
  //double beta2 = g_weights_beta[1];  
  boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);
  //boost::math::beta_distribution<> beta_dist(beta1, beta2);
  //g_weights_beta.print();
  //boost::math::beta_distribution<> beta_dist(1, 25);
  //std::vector<double> GVec0(t_n);
  arma::vec GVec(t_n);
  double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy, w0G2_cond;
  bool isSPAConverge, is_gtilde;
  arma::vec P2Vec(t_n);

  //std::vector<uint> indexZeroVec;
  //std::vector<uint> indexNonZeroVec;
  double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;
  arma::rowvec G1tilde_P_G2tilde_Vec;
  bool isCondition = false;
  for(unsigned int i = 0; i < q; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::vector<uint> indexZeroVec;
    std::vector<uint> indexNonZeroVec;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;

    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false; 
    std::string t_genoIndex_str = t_genoIndex.at(i);

    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );
    std::remove(end);
    bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, GVec, isImpute);
     //arma::vec GVec(GVec0);
     //GVec0.clear();
    if(!isReadMarker){
      break;
    }

    std::string info = chr+":"+std::to_string(pd)+"_"+ref+"/"+alt;

  double MAF = std::min(altFreq, 1 - altFreq);
  double MAC = MAF * 2 * t_n * (1 - missingRate);

  bool hasVarRatio;
  if((ptr_gSAIGEobj->m_varRatio).n_elem == 1){
	ptr_gSAIGEobj->assignSingleVarianceRatio();	
  }else{
	hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC);
	if(!hasVarRatio){
		std::cout << "Error! Conditioning marker " << info << " has MAC " << MAC << " and does not have variance ratio estimated." << std::endl;
		exit(EXIT_FAILURE);
	}	
  }	  
  
  flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);


 arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
       indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
       indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);



  MAF = std::min(altFreq, 1 - altFreq);


  arma::vec gtildeVec;
   Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
                    indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec);
      P1Mat.row(i) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
      P2Mat.col(i) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
      //P1Mat.row(i) = gtildeVec.t();
      //P2Mat.col(i) = P2Vec;
     MAFVec(i) = MAF;
     //w0G2_cond = boost::math::pdf(beta_dist, MAF);

     if(!t_weight_cond.is_zero()){
	 w0G2_cond = t_weight_cond(i);
    }else{
	 w0G2_cond = boost::math::pdf(beta_dist, MAF);
    }
     w0G2_cond_Vec(i) = w0G2_cond;
     gyVec(i) = gy * w0G2_cond;
     gsumVec = gsumVec + GVec * w0G2_cond;
     TstatVec(i) = Tstat;
     pVec(i) = pval;
  }
  arma::mat VarMat = P1Mat * P2Mat;

  VarInvMat = VarMat.i();   
  double qsum = arma::accu(gyVec);
  arma::vec gsumtildeVec; 
  ptr_gSAIGEobj->getadjG(gsumVec, gsumtildeVec);
  ptr_gSAIGEobj->assignConditionFactors(
		   			P2Mat,
					VarInvMat,
					VarMat,
					TstatVec,
				        w0G2_cond_Vec,	
					MAFVec,
					qsum,
					gsumtildeVec,
					pVec);

}

// [[Rcpp::export]]
void assign_conditionMarkers_factors_binary_region(
			   arma::vec & scalefactor_G2_cond){
	//std::cout << "assign_conditionMarkers_factors_binary_region" << std::endl;
	ptr_gSAIGEobj->assignConditionFactors_scalefactor(scalefactor_G2_cond);
}

// [[Rcpp::export]]
void set_iterator_inVcf(std::string & variantList, std::string & chrom, int & beg_pd, int & end_pd){
   if(!variantList.empty()){
	ptr_gVCFobj->set_iterator(variantList);	
   }else{
	ptr_gVCFobj->set_iterator(chrom, beg_pd, end_pd);
   }	   
}	

// [[Rcpp::export]]
bool check_Vcf_end(){
	bool isEnd = false;
	isEnd = ptr_gVCFobj->check_iterator_end();
	return(isEnd);
}


// [[Rcpp::export]]
void move_forward_iterator_Vcf(int i){
	ptr_gVCFobj->move_forward_iterator(i);
}



// [[Rcpp::export]]
arma::vec fast_logistf_fit(arma::mat & x,
		arma::vec & y,
		arma::vec & weight,
		arma::vec & offset,
		bool firth,
		arma::uvec & col_fit,
    	arma::vec init, 
	int maxit, 
	int maxstep, 
	int maxhs, 
	double lconv, 
	double gconv, 
	double xconv){
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  arma::vec pi_0 = -x * beta - offset; 
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  while(iter <= maxit){
	beta_old = beta;
	arma::vec wpi = weight % pi % (1 - pi);
	arma::vec wpi_sqrt = arma::sqrt(wpi);
	arma::vec W2 = weight % wpi_sqrt;
	arma::mat XW2(n, k, arma::fill::zeros);
	for(int j = 0; j < k; j++){
                XW2.col(j) = x.col(j) % W2;
        }       

	arma::mat Q;
	arma::mat R;
	arma::qr_econ(Q, R, XW2);
	arma::vec h = Q % Q * oneVec;
	arma::vec U_star(2, arma::fill::zeros);
	arma::vec ypih;
	if(firth){
		ypih = (weight % (y - pi)) + (h % (0.5 - pi));
	}else{
		ypih = (weight % (y - pi));
	}
	//ypih.print();
	arma::vec xcol(n, arma::fill::zeros);
	U_star = x.t() * ypih;
	
	arma::mat XX_XW2(n, k, arma::fill::zeros);
	for(int j = 0; j < k; j++){
		xcol = x.col(j);
		XX_XW2.col(j) = xcol % wpi_sqrt; 	
	}
	arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
	bool isinv = arma::inv_sympd (XX_covs, XX_Fisher); 
	if(!isinv){
		break;
	}	
	//}
	arma::vec delta = XX_covs * U_star;
	delta.replace(arma::datum::nan, 0);	

	double mx = arma::max(arma::abs(delta))/maxstep;
	if(mx > 1){
		delta = delta/mx;
	}
	evals = evals + 1;
	iter = iter + 1;
	beta = beta + delta;
	pi_0 = -x * beta - offset;
  	pi_0 = arma::exp(pi_0) + 1;
  	pi = 1/pi_0;
	if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		break;
	}
  }
	arma::mat var;
	if(XX_covs.has_nan()){
		var = XX_covs;
		beta = arma::datum::nan;
	}
	return beta;
}


