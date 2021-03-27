// [[Rcpp::depends(BH)]]
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <numeric>
#include <utility>
#include <stdexcept>
#include <memory>
#include <time.h>
#include <stdint.h>
#include <zlib.h>

#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/normal.hpp>
#include "ScoreTest.hpp"
#include "SPA.hpp"

#include <Rcpp.h>

typedef uLong uLongf;

//A global variable for the dosage file to test
FILE *m_fin;
std::vector <unsigned char> m_buf;
std::vector <unsigned char> m_zBuf;
uint m_zBufLens;
uint m_bufLens;

int gmtest_samplesize;
Rcpp::IntegerVector gm_sample_idx;

//genoToTest_bgenDosage
bool m_isQuery;
int m_N, m_N0, m_M, m_M0;
uint Nbgen;

std::vector< int > m_markerIndicesToInclude;
bool isReadVariantBgen = true;
//bool isOutputHetHomCountsinCaseCtrl = false;
bool m_isDropMissingDosages = false;
//double bgenMinMAF = 0;
//double bgenMinINFO = 0;
uint numSamples_bgen;
bool m_isSparseDosages = false;

ScoreClass ScoreTestObj;


double get_wall_time2(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


double get_cpu_time2(){
    return (double)clock() / CLOCKS_PER_SEC;
}


// [[Rcpp::export]]
void setIsDropMissingDosages_bgen (bool isdropmissingdosages){
  m_isDropMissingDosages = isdropmissingdosages;
}

// [[Rcpp::export]]
void setIsSparseDosage_bgen (bool isSparseDosage){
   m_isSparseDosages = isSparseDosage;
}


// [[Rcpp::export]]
void setMarkerIndicesToInclude (std::vector< int > & markerIndicesToInclude){
    m_markerIndicesToInclude = markerIndicesToInclude;
}

// [[Rcpp::export]]
int setgenoTest_bgenDosage(const std::string t_bgenFileName,
                 const std::string t_bgenFileIndex)
  {
    int mtotest;
    m_isQuery = false;
    if(t_bgenFileIndex == ""){
     m_isQuery = false;
     std::cout << "no index file for bgen is provided" << std::endl;
    }
    
   //else if(
   //  m_markerIndicesToInclude.size() == 0){	    
   //  std::cout << "no query list is provided" << std::endl;
   //  m_isQuery = false;
   // }else{
   //  m_isQuery = true;
   //  mtotest = m_markerIndicesToInclude.size();
   // }

    /****code from BOLT-LMM v2.3.4***/

    /********** READ HEADER v1.2**********/
    m_fin = fopen(t_bgenFileName.c_str(), "rb");
    uint offset; fread(&offset, 4, 1, m_fin); //cout << "offset: " << offset << endl;
    uint L_H; fread(&L_H, 4, 1, m_fin); //cout << "L_H: " << L_H << endl;
    uint m_M0; fread(&m_M0, 4, 1, m_fin); std::cout << "snpBlocks (m_M0): " << m_M0 << std::endl;
    assert(m_M0 != 0);
    if(!m_isQuery){mtotest = m_M0;};
    //uint Nbgen; fread(&Nbgen, 4, 1, m_fin); std::cout << "samples (Nbgen): " << Nbgen << std::endl;
    fread(&Nbgen, 4, 1, m_fin); std::cout << "samples (Nbgen): " << Nbgen << std::endl;
    numSamples_bgen = Nbgen; 
    //uint m_Nsample = t_SampleInBgen.size();
    //if (Nbgen != m_Nsample) {
    //  std::cerr << "ERROR: Number of samples in BGEN header does not match sample file" << std::endl;
    //  exit(1);
    //}
    char magic[5]; fread(magic, 1, 4, m_fin); magic[4] = '\0'; //cout << "magic bytes: " << string(magic) << endl;
    fseek(m_fin, L_H-20, SEEK_CUR); //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
    uint flags; fread(&flags, 4, 1, m_fin); //cout << "flags: " << flags << endl;
    uint CompressedSNPBlocks = flags&3; std::cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << std::endl;
    assert(CompressedSNPBlocks==1); // REQUIRE CompressedSNPBlocks==1
    uint Layout = (flags>>2)&0xf; std::cout << "Layout: " << Layout << std::endl;
    assert(Layout==1 || Layout==2); // REQUIRE Layout==1 or Layout==2
    fseek(m_fin, offset+4, SEEK_SET);
    return(m_M0);
  }





void Parse2(unsigned char *buf, uint bufLen, const unsigned char *zBuf, uint zBufLen,std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<int> & indexforMissing, double & info, std::vector<unsigned int> & iIndex, bool isFlip) {

    uLongf destLen = bufLen;

 //double wall2in = get_wall_time2();
// double cpu2in  = get_cpu_time2();

    if (uncompress(buf, &destLen, zBuf, zBufLen) != Z_OK || destLen != bufLen) {
      std::cerr << "ERROR: uncompress() failed" << std::endl;
      exit(1);
    }
 //double wall3in = get_wall_time2();
 //double cpu3in  = get_cpu_time2();
 //std::cout << "Wall Time in Parse2_1 = " << wall3in - wall2in << std::endl;
 //std::cout << "CPU Time  in Parse2_1= " << cpu3in - cpu2in  << std::endl;
      //}




    unsigned char *bufAt = buf;
    uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
    //std::cout << "N= " << N << std::endl;
    //std::cout << "Nbgen= " << Nbgen << std::endl; 
    if (N != Nbgen) {
      std::cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << std::endl;
      exit(1);
    }
    uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
    if (K != 2U) {
      std::cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << std::endl;
      exit(1);
    }
    uint Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
      std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << std::endl;
      exit(1);
    }
    uint Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
      std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << std::endl;
      exit(1);
    }
    const unsigned char *ploidyMissBytes = bufAt;
    for (uint i = 0; i < N; i++) {
      uint ploidyMiss = *bufAt; bufAt++;
      if (ploidyMiss != 2U && ploidyMiss != 130U) {
        std::cerr << "ERROR: " << snpName << " has ploidy/missingness byte = " << ploidyMiss
             << " (not 2 or 130)" << std::endl;
        exit(1);
      }
    }
    uint Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
      std::cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << std::endl;
      exit(1);
    }
    uint B = *bufAt; bufAt++;
    if (B != 8U) {
      std::cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << std::endl;
      exit(1);
    }

    //std::cout << "OKKKK" << std::endl;
    // Parse
    double lut[256];
    for (int i = 0; i <= 255; i++)
      lut[i] = i/255.0;

    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0, sum_fij_minus_eij2_sub = 0; // for INFO
    double p11,p10,p00,dosage,eij,fij, eijsub, fijsub;
    dosages.clear();
    dosages.reserve(m_N);
    if(!m_isSparseDosages){
      dosages.resize(m_N);
    }
    std::size_t missing_cnt = 0;

    for (uint i = 0; i < N; i++) {
     //if(i == 1){std::cout << "ploidyMissBytes[i] " << ploidyMissBytes[i] << std::endl;}
     if (ploidyMissBytes[i] != 130U){
      //bufAt += 2;
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      p00 = 1 - p11 - p10; //can remove
      dosage = 2*p11 + p10;

        eij = dosage;
        fij = 4*p11 + p10;
        sum_eij += eij;
        sum_fij_minus_eij2 += fij - eij*eij;
        if(gm_sample_idx[i] >= 0){
          if(!m_isSparseDosages){
              dosages[gm_sample_idx[i]] = 2 - dosage;
	      if(2 - dosage > 0){
                iIndex.push_back(gm_sample_idx[i]);
	      }	      
          }else{
              if(2 - dosage > 0){
                dosages.push_back(2 - dosage);
                iIndex.push_back(gm_sample_idx[i]);
              }
          }
          sum_eij_sub += eij;
        }
     }else if(ploidyMissBytes[i] == 130U){
        bufAt += 2;
        if(gm_sample_idx[i] >= 0){
          indexforMissing.push_back(gm_sample_idx[i]);
          ++missing_cnt;
          if(!m_isSparseDosages){
            dosages[gm_sample_idx[i]] = -1;
          }
        }
     }
    }
    //std::cout << "sum_eij_sub: " << sum_eij_sub << std::endl;
    AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;
    if(m_N == missing_cnt){
      AF = 0;
    }else{
      AF = AC/ 2/ ((double) (m_N - missing_cnt)) ;
    }

    double thetaHat = sum_eij / (2* (m_N - missing_cnt));
    //std::cout << "sum_eij " << sum_eij << std::endl;
    //std::cout << "missing_cnt " << sum_eij << std::endl;
    info = thetaHat==0 || thetaHat==1 ? 1 :
    1 - sum_fij_minus_eij2 / (2*(m_N - missing_cnt)*thetaHat*(1-thetaHat));

    isFlip = false;
    if(AF > 0.5){
      isFlip = true;
      iIndex.clear();
      AF = 1 - AF;
      uint dosagesSize = dosages.size();
      for (uint i = 0; i < dosagesSize; i++) {
	  if(dosages[i] != -1){
            dosages[i] = 2 - dosages[i];
	    if(dosages[i] > 0){
	      iIndex.push_back(gm_sample_idx[i]);
	    }
	  }
      }
      AC = 2*(double) (m_N - missing_cnt) - AC;      
    }

    //std::cout << "OKKKK1" << std::endl;
    if(missing_cnt > 0){
      //std::cout << "sample index with missing dosages for snpName " << snpName << " :";

      if(!m_isDropMissingDosages){
	//std::cout << "AC: " << AC << std::endl;
        double imputeDosage = 2*AF;
        for (unsigned int i = 0; i < indexforMissing.size(); i++)
        {
	  //std::cout << indexforMissing[i]+1 << ",";	
          if(!m_isSparseDosages){
            dosages[indexforMissing[i]] = imputeDosage;
	    iIndex.push_back(indexforMissing[i]);
          }else{
            dosages.push_back(imputeDosage);
            iIndex.push_back(indexforMissing[i]);
          }
            //std::cout << indexforMissing[i]+1 << ",";
          AC = AC + imputeDosage;
        }
	//std::cout << " " << std::endl;
	//std::cout << "AC new: " << AC << std::endl;
      }
    }

    //std::cout << "OKKKK2" << std::endl;
    //return(info);
  }






// [[Rcpp::export]]
Rcpp::List getOneMarker(int t_fileStartPos)
  {
    using namespace Rcpp; 
    //if(t_fileStartPos > 0){
   //	    m_isQuery = true;
    //}else{
//	    m_isQuery = false;
  //  }
    //std::cout << "m_isQuery: " << m_isQuery << std::endl;
    if(t_fileStartPos > 0){fseek(m_fin, t_fileStartPos, SEEK_SET);}
     std::string SNPID, RSID, chromosome, first_allele,second_allele ;
     uint position;
     std::vector< std::string > alleles ;
     std::vector< double > dosages;
     double AC, AF, info;
     std::vector< int > indexforMissing;
  std::vector< uint > iIndex;


    char snpID[65536], rsID[65536], chrStr[65536];
       uint maxLA = 65536, maxLB = 65536;
    char *allele1, *allele0;
   allele1 = (char *) malloc(maxLA+1);
    allele0 = (char *) malloc(maxLB+1);
    ushort LS; size_t numBoolRead = fread(&LS, 2, 1, m_fin); // cout << "LS: " << LS << " " << std::flush;
    bool isBoolRead;
    Rcpp::List result ;
    if ( numBoolRead > 0 ) {
      isBoolRead = true;
      fread(snpID, 1, LS, m_fin); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
      ushort LR; fread(&LR, 2, 1, m_fin); // cout << "LR: " << LR << " " << std::flush;
      fread(rsID, 1, LR, m_fin); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;
      RSID = std::string(rsID)=="." ? snpID : rsID;
      //std::string SNPID = string(snpID);

      ushort LC; fread(&LC, 2, 1, m_fin); // cout << "LC: " << LC << " " << std::flush;
      fread(chrStr, 1, LC, m_fin); chrStr[LC] = '\0';
      chromosome  = std::string(chrStr);

      uint physpos; fread(&physpos, 4, 1, m_fin); // cout << "physpos: " << physpos << " " << std::flush;
      position = physpos;
      ushort K; fread(&K, 2, 1, m_fin); //cout << "K: " << K << endl;
      if (K != 2) {
        std::cerr << "ERROR: Non-bi-allelic variant found: " << K << " alleles" << std::endl;
        exit(1);
      }

      uint LA; fread(&LA, 4, 1, m_fin); // cout << "LA: " << LA << " " << std::flush;
      if (LA > maxLA) {
        maxLA = 2*LA;
        free(allele1);
        allele1 = (char *) malloc(maxLA+1);
      }
      fread(allele1, 1, LA, m_fin); allele1[LA] = '\0';
      first_allele = std::string(allele1);
      uint LB; fread(&LB, 4, 1, m_fin); // cout << "LB: " << LB << " " << std::flush;
      if (LB > maxLB) {
        maxLB = 2*LB;
        free(allele0);
        allele0 = (char *) malloc(maxLB+1);
      }
      fread(allele0, 1, LB, m_fin); allele0[LB] = '\0';
      second_allele = std::string(allele0);

      uint C; fread(&C, 4, 1, m_fin); //cout << "C: " << C << endl;
      if (C > m_zBuf.size()) m_zBuf.resize(C-4);
      //std::cout << "m_zBuf.size() " << m_zBuf.size() << std::endl;
      uint D; fread(&D, 4, 1, m_fin); //cout << "D: " << D << endl;
      m_zBufLens = C-4; m_bufLens = D;
      fread(&m_zBuf[0], 1, C-4, m_fin);
      AC = 0;
      AF = 0;
      info = 0;
      //std::cout << m_bufLens << " m_bufLens" << std::endl;
      //std::cout << m_zBufLens << " m_zBufLens" << std::endl;
      //std::cout << C << " C" << std::endl;
      //std::cout << D << " D" << std::endl;
      //std::cout << first_allele << " first_allele" << std::endl;
      //std::cout << second_allele << " second_allele" << std::endl;
      //std::cout << chromosome << " chromosome" << std::endl;
      //std::cout << RSID << " RSID" << std::endl;
      if (m_bufLens > m_buf.size()) m_buf.resize(m_bufLens); //fix the length
         // double wall0in = get_wall_time2();
// double cpu0in  = get_cpu_time2();
      bool isFlip;
      Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens, RSID, Nbgen, dosages, AC, AF, indexforMissing, info, iIndex, isFlip);
      
      //double wall1in = get_wall_time2();
// double cpu1in  = get_cpu_time2();
// std::cout << "Wall Time in Parse2 = " << wall1in - wall0in << std::endl;
// std::cout << "CPU Time  in Parse2= " << cpu1in - cpu0in  << std::endl;
      //}
        Rcpp::DataFrame variants = Rcpp::DataFrame::create(
                Named("chromosome") = chromosome,
                Named("position") = position,
                Named("rsid") = RSID,
        //        Named("number_of_alleles") = number_of_allele,
                Named("allele0") = first_allele,
                Named("allele1") = second_allele,
                _["stringsAsFactors"] = false,
                Named("AC") = AC,
                Named("AF") = AF
               // Named("info") = info
                //Named("homN_cases") = homN_cases,
                //Named("hetN_cases") = hetN_cases,
                //Named("homN_ctrls") = homN_ctrls,
                //Named("hetN_ctrls") = hetN_ctrls
        );
    result[ "variants" ] = variants ;
    result["info"] = info;
    result[ "dosages" ] = dosages ;
    result["iIndex"] = iIndex;
    result["indexforMissing"] = indexforMissing;
    result["isFlip"] = isFlip;

    }else{
      isBoolRead = false;
      Rcpp::DataFrame variants = NULL;
      result["isBoolRead"] = isBoolRead;
    }

    dosages.clear();
    indexforMissing.clear();
    return(result);
  }





// [[Rcpp::export]]
bool getQueryStatus()
{
  return(m_isQuery);
}

// [[Rcpp::export]]
bool setQueryStatus(bool isQuery)
{
  m_isQuery = isQuery;
}

// [[Rcpp::export]]
bool getisReadVariantBgen()
{
  return(isReadVariantBgen);
}


/*
// [[Rcpp::export]]
double getMarkerInfo()
{
  return(markerInfo);
}
*/
/**************************************
        By SLEE 09/06/17
*************************************/

// [[Rcpp::export]]
void SetSampleIdx(Rcpp::IntegerVector sample_idx, int Ntest){
	gmtest_samplesize = Ntest;
        m_N = gmtest_samplesize;
	gm_sample_idx = sample_idx;
	//cc_idx = cc_index;
}


// [[Rcpp::export]]
void closetestGenoFile_bgenDosage() //needs further check
{
  
//  gm_sample_idx.erase();

  printf("closed the genofile!\n");

}


// [[Rcpp::export]]
int getSampleSizeinBgen(){
        return(numSamples_bgen);
}


// [[Rcpp::export]]
void assignforScoreTest_R(bool t_LOCO, std::vector<bool> & t_LOCOVec,  arma::mat & t_XVX, arma::mat & t_XXVX_inv,  arma::mat & t_XV, arma::mat & t_XVX_inv_XV, arma::mat & t_X, arma::vec & t_S_a,  arma::vec & t_res,  arma::vec & t_mu2, arma::vec & t_mu, double t_varRatio, arma::vec & t_tauvec, std::string t_traitType, bool t_isOutputAFinCaseCtrl, bool t_isOutputHetHomCountsinCaseCtrl, arma::vec & t_y){
	std::cout << "inside" << std::endl;
        ScoreTestObj.assignforScoreTest(t_LOCO, t_LOCOVec, t_XVX, t_XXVX_inv, t_XV, t_XVX_inv_XV, t_X, t_S_a, t_res, t_mu2, t_mu, t_varRatio, t_tauvec, t_traitType, t_isOutputAFinCaseCtrl, t_isOutputHetHomCountsinCaseCtrl, t_y); 
}

/*
// [[Rcpp::export]]
Rcpp::List getScoreTest(int t_fileStartPos) {
   Rcpp::List Glist;
   Glist = getOneMarker(t_fileStartPos);
   double t_Beta, t_seBeta, t_altFreq;
   std::string t_pval_str;
   arma::uvec iIndexvec = Rcpp::as<arma::uvec>(Glist["iIndex"]);
   std::cout << iIndexvec.n_elem << std::endl;
   //t_altFreq = Glist["variants"]["AF"];

   double wall2in = get_wall_time2();
   double cpu2in  = get_cpu_time2();

 //ScoreTestObj.scoreTest(Glist["dosages"], t_Beta, t_seBeta, t_pval_str, t_altFreq);
 double wall3in = get_wall_time2();
 double cpu3in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTest = " << wall3in - wall2in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTest= " << cpu3in - cpu2in  << std::endl;
 ScoreTestObj.scoreTestFast(Glist["dosages"], iIndexvec, t_Beta, t_seBeta, t_pval_str, t_altFreq);
 double wall4in = get_wall_time2();
 double cpu4in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTestFast = " << wall4in - wall3in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTestFast= " << cpu4in - cpu3in  << std::endl;
   Glist["Beta"] = t_Beta;
   Glist["se"] = t_seBeta;
   Glist["pval"] = t_pval_str;
   Glist["dosages"] = NULL;
   Glist["indexforMissing"] = NULL;
   Glist["iIndex"] = NULL;
   return(Glist);
}
*/
// [[Rcpp::export]]
Rcpp::List getScoreTest_SPA(int t_fileStartPos, std::string traitType) {
   Rcpp::List Glist;
   Glist = getOneMarker(t_fileStartPos);
   double t_Beta, t_seBeta, t_altFreq, t_Tstat, t_var1, t_var2;
   std::string t_pval_str;
   arma::uvec iIndexvec = Rcpp::as<arma::uvec>(Glist["iIndex"]);
   std::cout << iIndexvec.n_elem << std::endl;
   Rcpp::DataFrame variantsDF;
   variantsDF = Glist["variants"];
   t_altFreq = variantsDF["AF"];

   double wall2in = get_wall_time2();
   double cpu2in  = get_cpu_time2();
   bool isScoreFast;
   arma::vec dosages = Glist["dosages"];
   arma::vec t_gtilde;
   double wall3in,cpu3in, wall4in, cpu4in, wall5in, cpu5in, wall6in, cpu6in, wall7in, cpu7in, wall5inb, cpu5inb, wall5inc, cpu5inc;
 wall3in = get_wall_time2();
 cpu3in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTest = 3-2" << wall3in - wall2in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTest= " << cpu3in - cpu2in  << std::endl;
   if(t_altFreq > 0.05){
     isScoreFast = false;	   
     ScoreTestObj.scoreTest(dosages, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde);
    }else{
      isScoreFast = true;
      ScoreTestObj.scoreTestFast(dosages, iIndexvec, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2);
    }
 wall4in = get_wall_time2();
 cpu4in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTestFast = 4-3" << wall4in - wall3in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTestFast= " << cpu4in - cpu3in  << std::endl;
  double pval;
  bool isSPAConverge;
  double pval_noadj=std::stod(t_pval_str);
 if(pval_noadj < 0.05 && traitType != "quantitative"){
 //  SPA_binary_fast(m_mu)
 //  std::cout << "HERE SPA" << std::endl;
   arma::vec t_mu; 
   ScoreTestObj.get_mu(t_mu); 
   if(isScoreFast){
     t_gtilde = ScoreTestObj.getadjG(dosages);
   }
   double q, qinv, m1;
 wall5in = get_wall_time2();
 cpu5in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTestFast = 5-4" << wall5in - wall4in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTestFast= " << cpu5in - cpu4in  << std::endl;
   m1 = dot(t_mu, t_gtilde);
   arma::vec gNB = t_gtilde.elem(iIndexvec);
   arma::vec gNA = t_gtilde;
   gNA.shed_rows(iIndexvec);
   arma::vec muNB = t_mu.elem(iIndexvec);
   arma::vec muNA = t_mu;
   muNA.shed_rows(iIndexvec);
   double NAmu= m1-dot(gNB,muNB);
   double NAsigma;
   
   if(traitType == "binary"){
	   std::cout << "OKKKKKKK SPA" << std::endl;
		q = t_Tstat/sqrt(t_var1/t_var2) + m1;
			
		if((q-m1) > 0){
			qinv = -1 * std::abs(q-m1) + m1;
		}else if ((q-m1) == 0){
			qinv =  m1;
		}else{
			qinv = std::abs(q-m1) + m1;
		}
           NAsigma = t_var2- arma::sum(muNB % (1-muNB) % arma::pow(gNB,2));
	}else if(traitType == "survival"){
		q = t_Tstat/sqrt(t_var1/t_var2);
		qinv = -q;
           NAsigma = t_var2- arma::sum(muNB % arma::pow(gNB,2)); 
	}	

   //SPA(t_mu, t_gtilde, q, qinv, pval_noadj, 1e-5, false, traitType, pval, isSPAConverge);
    bool logp=false;
    SPA_fast(t_mu, t_gtilde, q, qinv, pval_noadj, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, 1e-5, traitType, pval, isSPAConverge);
    boost::math::normal ns;
    double t_qval;
    //if(!logp){
      t_qval = fabs(boost::math::quantile(ns, pval/2));
      t_seBeta = t_Beta/t_qval;
    //}else{
    //  t_qval = boost::math::quantile(ns, )
    //}	    

    wall6in = get_wall_time2();
 cpu6in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTestFast = 6-5" << wall6in - wall5in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTestFast= " << cpu6in - cpu5in  << std::endl;
//	  Rcpp::List SPA_binary_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol){
 }
 uint dosageSize=dosages.n_elem;
  double AF_case, AF_ctrl;
  int N_case_het, N_case_hom, N_ctrl_het, N_ctrl_homi, N_case, N_ctrl;   
  arma::vec dosage_case;
  arma::vec dosage_ctrl;
  arma::ivec N_case_ctrl_het_hom(4);
  arma::vec AF_case_ctrl(2);
  if(ScoreTestObj.m_isOutputAFinCaseCtrl || ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
    dosage_case = dosages.elem(ScoreTestObj.m_case_indices);
    dosage_ctrl = dosages.elem(ScoreTestObj.m_ctrl_indices);
    AF_case = arma::mean(dosage_case) /2 ;
    AF_ctrl = arma::mean(dosage_ctrl) /2 ;
    AF_case_ctrl[0] = AF_case;
    AF_case_ctrl[1] = AF_ctrl;
    Glist["AF_case_ctrl"] = AF_case_ctrl;
    N_case = dosage_case.n_elem;
    N_ctrl = dosage_ctrl.n_elem;

  }
  arma::uvec N_case_ctrl_het_hom0; 
  if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
    N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
    N_case_ctrl_het_hom[0] = N_case_ctrl_het_hom0.n_elem;
//N_case_ctrl_het_hom[0] = (arma::find(dosage_case <= 2 && dosage_case >=1.5)).n_elem;
    N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case > 0.5);
    N_case_ctrl_het_hom[1] = N_case_ctrl_het_hom0.n_elem;
//    N_case_ctrl_het_hom[1] = arma::find(dosage_case < 1.5 && dosage_case > 0.5);
    N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
   N_case_ctrl_het_hom[2] = N_case_ctrl_het_hom0.n_elem; 
    //N_case_ctrl_het_hom[2] = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
   N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl > 0.5);
   N_case_ctrl_het_hom[3] = N_case_ctrl_het_hom0.n_elem;
//   N_case_ctrl_het_hom[3] = arma::find(dosage_ctrl < 1.5 && dosage_ctrl > 0.5);	 
    Glist["N_case_ctrl_het_hom"] = N_case_ctrl_het_hom;
  }

  if(Glist["isFlip"]){
   t_Beta = -t_Beta;
   t_Tstat = -t_Tstat;
   double AF, AC;
   AF = variantsDF["AF"];
   AC = variantsDF["AC"];
   variantsDF["AF"] = 1-AF;
   variantsDF["AC"] = 2*double(dosageSize) - AC;
   Glist["variants"] = variantsDF;
   if(ScoreTestObj.m_isOutputAFinCaseCtrl || ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
     AF_case_ctrl = 1-AF_case_ctrl;
     Glist["AF_case_ctrl"] = AF_case_ctrl;
   }
   if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
     N_case_ctrl_het_hom[0] = N_case - N_case_ctrl_het_hom[0] - N_case_ctrl_het_hom[1];
     N_case_ctrl_het_hom[2] = N_ctrl - N_case_ctrl_het_hom[2] - N_case_ctrl_het_hom[3];     
     Glist["N_case_ctrl_het_hom"] = N_case_ctrl_het_hom;
   }	   

 }	 
   Glist["Beta"] = t_Beta;
   Glist["se"] = t_seBeta;
   Glist["pval"] = t_pval_str;
   Glist["Tstat"] = t_Tstat;
   Glist["var1"] = t_var1;
   Glist["var2"] = t_var2;
   Glist["dosages"] = NULL;
   Glist["indexforMissing"] = NULL;
   Glist["iIndex"] = NULL;
   if(traitType!="quantitative"){
   	if(pval_noadj < 0.05){
     		Glist["pval_SPA"] = pval;
   	}else{
		Glist["pval_SPA"] = t_pval_str;
	}	
     	Glist["isSPAConverge"] = isSPAConverge;
   }

   return(Glist);
}
