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
#include "getMem.hpp"

#include <Rcpp.h>

using namespace Rcpp;

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

double m_minInfo, m_minMAF;

ScoreClass ScoreTestObj;

// [[Rcpp::export]]
void set_minInfo_minMAF(double t_minInfo, double t_minMAF){
	m_minInfo = t_minInfo;
	m_minMAF = t_minMAF;
}	


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





void Parse2(unsigned char *buf, uint bufLen, const unsigned char *zBuf, uint zBufLen,std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<int> & indexforMissing, double & info, std::vector<unsigned int> & iIndex, bool & isFlip) {

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
	      //iIndex.push_back(gm_sample_idx[i]);
	      iIndex.push_back(i);
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
Rcpp::List getOneMarker_old(int t_fileStartPos)
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
      bool isFlip = false;
      Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens, RSID, Nbgen, dosages, AC, AF, indexforMissing, info, iIndex, isFlip);
      
      //double wall1in = get_wall_time2();
// double cpu1in  = get_cpu_time2();
// std::cout << "Wall Time in Parse2 = " << wall1in - wall0in << std::endl;
// std::cout << "CPU Time  in Parse2= " << cpu1in - cpu0in  << std::endl;
      //}
      /*
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
	*/
     result["chromosome"]= chromosome;
     result["position"] = position;
     result["rsid"] =RSID;
     result["allele0"] =first_allele;
     result["allele1"] =second_allele;
     result["AC"] = AC;
     result["AF"] = AF;
    //result[ "variants" ] = variants ;
    result["info"] = info;
    result[ "dosages" ] = dosages ;
    result["iIndex"] = iIndex;
    result["indexforMissing"] = indexforMissing;
    result["isFlip"] = isFlip;
    result["isBoolRead"] = isBoolRead;

    }else{
      isBoolRead = false;
      Rcpp::DataFrame variants = NULL;
      result["isBoolRead"] = isBoolRead;
      result["info"] = 1;
    }

    return(result);
  }


// [[Rcpp::export]]
void getOneMarker(int t_fileStartPos,
		     std::string & SNPID,
		     std::string & RSID, 
		     std::string & chromosome, 
		     std::string & first_allele,
		     std::string & second_allele,
     			uint & position,
     		    std::vector< double > & dosages,
     			double & AC,
		       	double & AF, 
			double & info,
     		std::vector< int > & indexforMissing,
  		std::vector< uint > & iIndex,
		bool & isFlip,
		bool & isBoolRead)
  {
    using namespace Rcpp;
    //if(t_fileStartPos > 0){
   //       m_isQuery = true;
    //}else{
//          m_isQuery = false;
  //  }
    //std::cout << "m_isQuery: " << m_isQuery << std::endl;
    if(t_fileStartPos > 0){fseek(m_fin, t_fileStartPos, SEEK_SET);}

    char snpID[65536], rsID[65536], chrStr[65536];
       uint maxLA = 65536, maxLB = 65536;
    char *allele1, *allele0;
   allele1 = (char *) malloc(maxLA+1);
    allele0 = (char *) malloc(maxLB+1);
    ushort LS; size_t numBoolRead = fread(&LS, 2, 1, m_fin); // cout << "LS: " << LS << " " << std::flush;
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
      free(allele1);	
      uint LB; fread(&LB, 4, 1, m_fin); // cout << "LB: " << LB << " " << std::flush;
      if (LB > maxLB) {
        maxLB = 2*LB;
        free(allele0);
        allele0 = (char *) malloc(maxLB+1);
      }
      fread(allele0, 1, LB, m_fin); allele0[LB] = '\0';
      second_allele = std::string(allele0);
      free(allele0);
      uint C; fread(&C, 4, 1, m_fin); //cout << "C: " << C << endl;
      if (C > m_zBuf.size()) m_zBuf.resize(C-4);
      //std::cout << "m_zBuf.size() " << m_zBuf.size() << std::endl;
      uint D; fread(&D, 4, 1, m_fin); //cout << "D: " << D << endl;
      m_zBufLens = C-4; m_bufLens = D;
      fread(&m_zBuf[0], 1, C-4, m_fin);
      AC = 0;
      AF = 0;
      info = 0;
      //std::cout << RSID << " RSID" << std::endl;
      if (m_bufLens > m_buf.size()) m_buf.resize(m_bufLens); //fix the length
         // double wall0in = get_wall_time2();
// double cpu0in  = get_cpu_time2();
      isFlip = false;
      Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens, RSID, Nbgen, dosages, AC, AF, indexforMissing, info, iIndex, isFlip);


    }else{
      isBoolRead = false;
      info = 1;
    }

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
void assignforScoreTest_R(bool t_LOCO, std::vector<bool> & t_LOCOVec,  arma::mat & t_XVX, arma::mat & t_XXVX_inv,  arma::mat & t_XV, arma::mat & t_XVX_inv_XV, arma::mat & t_X, arma::vec & t_S_a,  arma::vec & t_res,  arma::vec & t_mu2, arma::vec & t_mu, double t_varRatio, arma::vec & t_tauvec, std::string t_traitType, bool t_isOutputAFinCaseCtrl, bool t_isOutputHetHomCountsinCaseCtrl, bool t_isOutputNinCaseCtrl, arma::vec & t_y){
	std::cout << "inside" << std::endl;
        ScoreTestObj.assignforScoreTest(t_LOCO, t_LOCOVec, t_XVX, t_XXVX_inv, t_XV, t_XVX_inv_XV, t_X, t_S_a, t_res, t_mu2, t_mu, t_varRatio, t_tauvec, t_traitType, t_isOutputAFinCaseCtrl, t_isOutputHetHomCountsinCaseCtrl, t_isOutputNinCaseCtrl, t_y); 
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
Rcpp::List getScoreTest_SPA_old(int t_fileStartPos, std::string traitType) {
   Rcpp::List Glist;
   Glist = getOneMarker_old(t_fileStartPos);

   if(Glist["isBoolRead"]){


   double t_Beta, t_seBeta, t_altFreq, t_Tstat, t_var1, t_var2;
   std::string t_pval_str;
   //std::cout << iIndexvec.n_elem << std::endl;
   Rcpp::DataFrame variantsDF;
   variantsDF = Glist["variants"];
   t_altFreq = variantsDF["AF"];
    double m_info = Glist["info"];

   if(m_info >= m_minInfo && t_altFreq >= m_minMAF && t_altFreq <= 1-m_minMAF){


   bool isScoreFast;
/*
   double wall2in = get_wall_time2();
double cpu2in  = get_cpu_time2();
*/
   arma::vec dosages = Glist["dosages"];
   arma::vec t_gtilde;
   arma::uvec iIndexvec;
/*
   double wall3in,cpu3in, wall4in, cpu4in, wall5in, cpu5in, wall6in, cpu6in, wall7in, cpu7in, wall5inb, cpu5inb, wall5inc, cpu5inc;
wall3in = get_wall_time2();
cpu3in  = get_cpu_time2();
std::cout << "Wall Time in ScoreTestObj.scoreTest = 3-2 " << wall3in - wall2in << std::endl;
std::cout << "CPU Time  in ScoreTestObj.scoreTest= " << cpu3in - cpu2in  << std::endl;
*/
if(t_altFreq > 0.05){
     isScoreFast = false;	   
     ScoreTestObj.scoreTest(dosages, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde);
    }else{
      isScoreFast = true;
      iIndexvec = Rcpp::as<arma::uvec>(Glist["iIndex"]);
      ScoreTestObj.scoreTestFast(dosages, iIndexvec, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2);
    }
/*
wall4in = get_wall_time2();
 cpu4in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTestFast = 4-3 " << wall4in - wall3in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTestFast= " << cpu4in - cpu3in  << std::endl;
*/
 double pval;
  bool isSPAConverge;
  double pval_noadj=std::stod(t_pval_str);
 if(pval_noadj < 0.05 && traitType != "quantitative"){
 //  SPA_binary_fast(m_mu)
 //  std::cout << "HERE SPA" << std::endl;
   arma::vec t_mu; 
   ScoreTestObj.get_mu(t_mu); 
   if(isScoreFast){
     ScoreTestObj.getadjG(dosages, t_gtilde);
   }
   if(t_altFreq > 0.05){
	   iIndexvec = Rcpp::as<arma::uvec>(Glist["iIndex"]);
   }
   double q, qinv, m1;
/*
   wall5in = get_wall_time2();
 cpu5in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTestFast = 5-4 " << wall5in - wall4in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTestFast= " << cpu5in - cpu4in  << std::endl;
*/
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
//	   std::cout << "OKKKKKKK SPA" << std::endl;
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
/*
    wall6in = get_wall_time2();
 cpu6in  = get_cpu_time2();
 std::cout << "Wall Time in ScoreTestObj.scoreTestFast = 6-5  " << wall6in - wall5in << std::endl;
 std::cout << "CPU Time  in ScoreTestObj.scoreTestFast= " << cpu6in - cpu5in  << std::endl;
*/
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
   Glist["isTest"] = true;
  }else{
    Glist["isTest"] = false;
  }
   }else{
	Glist["isTest"] = false;
   }	   

   return(Glist);
}



// [[Rcpp::export]]
bool getScoreTest_SPA(int t_fileStartPos, std::string traitType,
		std::string & t_chromosome,
		uint & t_position,
		std::string & t_rsid,
		std::string & t_allele0,
		std::string & t_allele1,
		double & t_AF,
		double & t_AC,
		double & t_info,
		int & t_Ntest,
		double & t_Beta,
		double & t_se, 
		double & t_Tstat,
		double & t_var1, 
		double & t_var2, 
		std::string & t_noSPApval,
		double & t_SPApval,
		bool & t_isSPAConverge,
		double & t_AF_case,
		double & t_AF_ctrl,
		int & t_N_case_hom,
		int & t_N_case_het,
		int & t_N_ctrl_hom,
		int & t_N_ctrl_het,
		int & t_N_case, 
		int & t_N_ctrl){
          double vm2, rss2;
  
   std::string SNPID;
   std::vector< double > dosages;
   std::vector< int > indexforMissing;
   std::vector< uint > iIndex, iIndexComp;
   bool isFlip, isBoolRead;
   getOneMarker(t_fileStartPos, 
		SNPID,
	        t_rsid,
	        t_chromosome,
		t_allele0,
		t_allele1,
		t_position,
		dosages,
		t_AC,
		t_AF,
		t_info,
		indexforMissing,
		iIndex,
		isFlip, 
		isBoolRead);

  //std::cout << "isFlip " << isFlip << std::endl; 
//	  std::vector< uint > iIndex;
//	  std::vector< double > dosages;
  bool isTest = false;
//isBoolRead=false;
   if(isBoolRead){
   std::string t_pval_str;
   if(t_info >= m_minInfo && t_AF >= m_minMAF && t_AF <= 1-m_minMAF){

   bool isScoreFast;
   arma::vec dosagesVec(dosages);
   arma::vec t_gtilde;
   arma::uvec iIndexvec;
   iIndexvec = arma::conv_to<arma::uvec>::from(iIndex);
  
   t_Ntest = m_N; 
   if(t_AF > 0.05){
     	isScoreFast = false;	   
	//std::cout << "here0c1" << std::endl; 
     	ScoreTestObj.scoreTest(dosagesVec, t_Beta, t_se, t_pval_str, t_AF, t_Tstat, t_var1, t_var2, t_gtilde);
    }else{
      	isScoreFast = true;
	//std::cout << "here0c2" << std::endl; 
        ScoreTestObj.scoreTestFast(dosagesVec, iIndexvec, t_Beta, t_se, t_pval_str, t_AF, t_Tstat, t_var1, t_var2);
    }

  t_isSPAConverge = false;
  //std::cout << "t_pval_str: " << t_pval_str << std::endl;

    double pval_noadj;
    try {
        pval_noadj = std::stod(t_pval_str);
    } catch (const std::invalid_argument&) {
	pval_noadj = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
    } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
	pval_noadj = 0;
    }


  //double pval_noadj=std::stod(t_pval_str);
  t_noSPApval = t_pval_str;
  
 if(pval_noadj < 0.05 && traitType != "quantitative"){
   arma::vec t_mu; 
   ScoreTestObj.get_mu(t_mu); 
   if(isScoreFast){
     ScoreTestObj.getadjG(dosagesVec, t_gtilde);
     //std::cout << "isScoreFast" << std::endl;
   }
   double q, qinv, m1;
   arma::uvec iIndexComVec = arma::find(dosagesVec == 0);
   m1 = dot(t_mu, t_gtilde);
   arma::vec gNB = t_gtilde.elem(iIndexvec);
   arma::vec gNA = t_gtilde.elem(iIndexComVec);
   arma::vec muNB = t_mu.elem(iIndexvec);
   arma::vec muNA = t_mu.elem(iIndexComVec);
   //iIndexComVec.clear();
   double NAmu= m1-dot(gNB,muNB);
   double NAsigma;
   if(traitType == "binary"){
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

    bool logp=false;
    SPA_fast(t_mu, t_gtilde, q, qinv, pval_noadj, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, 1e-5, traitType, t_SPApval, t_isSPAConverge);
    boost::math::normal ns;
    double t_qval;
    //t_qval = fabs(boost::math::quantile(ns, t_SPApval/2));

        try { 
           t_qval = boost::math::quantile(ns, pval_noadj/2);
	   t_qval = fabs(t_qval);
	   t_se = fabs(t_Beta)/t_qval;
        }catch (const std::overflow_error&) {
          t_qval = std::numeric_limits<double>::infinity();
	  t_se = 0;
        }

 }
 uint dosageSize=dosagesVec.n_elem;
  int N_case, N_ctrl;   
  arma::vec dosage_case;
  arma::vec dosage_ctrl;
  if(ScoreTestObj.m_isOutputAFinCaseCtrl || ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
    dosage_case = dosagesVec.elem(ScoreTestObj.m_case_indices);
    dosage_ctrl = dosagesVec.elem(ScoreTestObj.m_ctrl_indices);
    t_AF_case = arma::mean(dosage_case) /2 ;
    t_AF_ctrl = arma::mean(dosage_ctrl) /2 ;
    t_N_case = dosage_case.n_elem;
    t_N_ctrl = dosage_ctrl.n_elem;

  }

  arma::uvec N_case_ctrl_het_hom0; 
  if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
    N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
    t_N_case_hom  = N_case_ctrl_het_hom0.n_elem;
    N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case > 0.5);
    t_N_case_het = N_case_ctrl_het_hom0.n_elem;
    N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
    t_N_ctrl_hom = N_case_ctrl_het_hom0.n_elem; 
    N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl > 0.5);
    t_N_ctrl_het= N_case_ctrl_het_hom0.n_elem;
    //N_case_ctrl_het_hom0.clear();
  }

  //if(ScoreTestObj.m_isOutputAFinCaseCtrl || ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
  // dosage_case.clear();
  // dosage_ctrl.clear(); 
  //}	 
  if(isFlip){
   t_Beta = -t_Beta;
   t_Tstat = -t_Tstat;
   t_AF = 1-t_AF;
   t_AC = 2*dosageSize - t_AC;
   if(ScoreTestObj.m_isOutputAFinCaseCtrl || ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
     t_AF_case = 1-t_AF_case;
     t_AF_ctrl = 1-t_AF_ctrl;
   }
   if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
     t_N_case_hom = N_case - t_N_case_het - t_N_case_hom;
     t_N_ctrl_hom = N_ctrl - t_N_ctrl_het - t_N_ctrl_hom;
   }	   

  }	 
   if(traitType!="quantitative"){
   	if(pval_noadj < 0.05){
     		t_SPApval = t_SPApval;
   	}else{
		t_SPApval = pval_noadj;
	}	
   }
   isTest= true;
  }else{
    isTest = false;
  }
}else{
  isTest = false;
}


   return(isTest);
}


// [[Rcpp::export]]
Rcpp::DataFrame getScoreTest_SPA_multi(int mth_start, int m_to_test, std::string traitType) {

    int markerIndex = 0;
  // set up output
  std::vector<std::string> chromVec(m_to_test);    // chromosome
  std::vector<uint> posVec(m_to_test);    // chromosome
  std::vector<std::string> rsIDVec(m_to_test);    // marker IDs
  std::vector<std::string> A1Vec(m_to_test);    // marker IDs
  std::vector<std::string> A2Vec(m_to_test);    // marker IDs
  std::vector<double> infoVec(m_to_test);      // marker information: CHR:POS:REF:ALT
  std::vector<double> ACVec(m_to_test);        // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altFreqVec(m_to_test);        // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> BetaVec(m_to_test);           // beta value for ALT allele
  std::vector<double> seBetaVec(m_to_test);         // seBeta value
  std::vector<std::string> pvalVec(m_to_test);           // p values
  std::vector<double> TstatVec(m_to_test);           // p values
  std::vector<double> var1Vec(m_to_test);           // p values
  std::vector<double> var2Vec(m_to_test);           // p values

  std::vector<int> NVec(m_to_test);           // p values
  std::vector<int> N_case_Vec(m_to_test);
  std::vector<int> N_ctrl_Vec(m_to_test);

  // p values
  //if(traitType == "binary"){
  std::vector<double> SPApvalVec(m_to_test);
  std::vector<bool> SPAConverge(m_to_test);
 //   if(ScoreTestObj.m_isOutputAFinCaseCtrl){ 
  std::vector<double> AFinCaseVec(m_to_test);
  std::vector<double> AFinCtrlVec(m_to_test);
 //   }

 //   if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
  std::vector<int> N_case_homVec(m_to_test);
  std::vector<int> N_ctrl_hetVec(m_to_test);
  std::vector<int> N_case_hetVec(m_to_test);
  std::vector<int> N_ctrl_homVec(m_to_test);
 //  }
  //}
		std::string t_chromosome;
                uint  t_position;
                std::string  t_rsid;
                std::string  t_allele0;
                std::string  t_allele1;
                double  t_AF;
                double  t_AC;
                double  t_info;
                int t_Ntest;
                double  t_Beta;
                double  t_se;
                double  t_Tstat;
                double  t_var1;
                double t_var2;
		std::string  t_noSPApval;
                double t_SPApval;
                bool t_isSPAConverge;
                double t_AF_case;
                double t_AF_ctrl;
                int t_N_case_hom;
                int t_N_case_het;
                int t_N_ctrl_hom;
                int t_N_ctrl_het;
                int t_N_case;
                int t_N_ctrl;

  
bool isTest;
    for (int i = 0; i < m_to_test; i++) {
      if(m_isQuery){
	      //std::cout << "mth_start " << mth_start << std::endl;
	      //std::cout << "i " << i << std::endl;
	      //std::cout << "mth_start+i " << mth_start+i << std::endl; 
        markerIndex = m_markerIndicesToInclude[mth_start+i];
      }else{
	      //std::cout << "m_isQuery " << m_isQuery << std::endl;
	markerIndex = 0;
      }	      
    isTest = getScoreTest_SPA(markerIndex, traitType, t_chromosome, t_position, t_rsid, t_allele0, t_allele1, t_AF, t_AC, t_info, t_Ntest, t_Beta, t_se, t_Tstat, t_var1, t_var2, t_noSPApval, t_SPApval, t_isSPAConverge, t_AF_case, t_AF_ctrl, t_N_case_hom, t_N_case_het, t_N_ctrl_hom, t_N_ctrl_het, t_N_case, t_N_ctrl);
    //std::cout << "isTest " << isTest << std::endl;


    if(isTest){
        chromVec.at(i) = t_chromosome;
	posVec.at(i) = t_position;
	rsIDVec.at(i) = t_rsid;
	A1Vec.at(i) = t_allele0;
	A2Vec.at(i) = t_allele1;
	altFreqVec.at(i) = t_AF;
	ACVec.at(i) = t_AC;
	infoVec.at(i) = t_info;
	NVec.at(i) = t_Ntest;
	BetaVec.at(i) = t_Beta;
        seBetaVec.at(i) = t_se;
        TstatVec.at(i) = t_Tstat;
	var1Vec.at(i) = t_var1;
	var2Vec.at(i) = t_var2;
	pvalVec.at(i) = t_noSPApval;
	if(traitType == "binary"){
		SPApvalVec.at(i) = t_SPApval;
		SPAConverge.at(i) = t_isSPAConverge;	
		if(ScoreTestObj.m_isOutputAFinCaseCtrl){
			AFinCaseVec.at(i) = t_AF_case;
			AFinCtrlVec.at(i) = t_AF_ctrl;	
		}
		if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
			N_case_homVec.at(i) = t_N_case_hom;
			N_case_hetVec.at(i) = t_N_case_het;
			N_ctrl_homVec.at(i) = t_N_ctrl_hom;
			N_ctrl_hetVec.at(i) = t_N_ctrl_het;
   		}
		if(ScoreTestObj.m_isOutputNinCaseCtrl){
			N_case_Vec.at(i) = t_N_case;
			N_ctrl_Vec.at(i) = t_N_ctrl;	
		}	
	}	

      }else{
        chromVec.at(i) = "-1";
      }	      
      
    }

    Rcpp::DataFrame OUT_DF = Rcpp::DataFrame::create(
                Named("CHR") = chromVec,
                Named("POS") = posVec,
                Named("rsid") =  rsIDVec,
                Named("Allele1") = A1Vec,
                Named("Allele2") = A2Vec,
                _["stringsAsFactors"] = false,
                Named("AC_Allele2") = ACVec,
                Named("AF_Allele2") = altFreqVec,
		Named("imputationInfo") = infoVec,
		Named("N") = NVec,
		Named("BETA") = BetaVec,
		Named("SE") = seBetaVec,
		Named("Tstat") = TstatVec

        );
	if(traitType=="binary"){
		OUT_DF["p.value"] = SPApvalVec;
		OUT_DF["p.value.NA"] = pvalVec;
		OUT_DF["Is.SPA.converge"] = SPAConverge;
		OUT_DF["varT"] = var1Vec;
		OUT_DF["varTstar"] = var2Vec;
		if(ScoreTestObj.m_isOutputNinCaseCtrl){
                        OUT_DF["N.Cases"] = AFinCaseVec;
                        OUT_DF["N.Controls"] = AFinCtrlVec;
                }
		if(ScoreTestObj.m_isOutputAFinCaseCtrl){
			OUT_DF["AF.Cases"] = N_case_Vec;
			OUT_DF["AF.Controls"] = N_ctrl_Vec;
                }



                if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
			OUT_DF["homN_Allele2_cases"] = N_case_homVec;
			OUT_DF["hetN_Allele2_cases"] = N_case_hetVec;
			OUT_DF["homN_Allele2_ctrls"] = N_ctrl_homVec;
			OUT_DF["hetN_Allele2_ctrls"] = N_ctrl_hetVec;
                }

	}else{
		OUT_DF["p.value"] = pvalVec;
		OUT_DF["varT"] = var1Vec;
                OUT_DF["varTstar"] = var2Vec;
	
	}

	return(OUT_DF);
}




// [[Rcpp::export]]
void getScoreTest_SPA_multi_new(int mth_start, int m_to_test, std::string traitType,
	       std::vector<std::string> & chromVec,
std::vector<uint> & posVec,
std::vector<std::string> & rsIDVec,
std::vector<std::string> & A1Vec,
std::vector<std::string> & A2Vec,
std::vector<double> & infoVec,
std::vector<double> & ACVec,
std::vector<double> & altFreqVec,
std::vector<double> & BetaVec,
std::vector<double> & seBetaVec,
std::vector<std::string> & pvalVec,
std::vector<double> & TstatVec,
std::vector<double> & var1Vec,
std::vector<double> & var2Vec,
std::vector<int> & NVec,
std::vector<double> & SPApvalVec,
std::vector<bool> & SPAConverge,
std::vector<double> & AFinCaseVec,
std::vector<double> & AFinCtrlVec,
std::vector<double> & N_case_homVec,
std::vector<double> & N_ctrl_hetVec,
std::vector<double> & N_case_hetVec,
std::vector<double> & N_ctrl_homVec	
		) {

    int markerIndex = 0;
      // set up output
               std::string t_chromosome;
                uint  t_position;
                std::string  t_rsid;
                std::string  t_allele0;
                std::string  t_allele1;
                double  t_AF;
                double  t_AC;
                double  t_info;
                int t_Ntest;
                double  t_Beta;
                double  t_se;
                double  t_Tstat;
                double  t_var1;
                double t_var2;
                std::string  t_noSPApval;
                double t_SPApval;
                bool t_isSPAConverge;
                double t_AF_case;
                double t_AF_ctrl;
                int t_N_case_hom;
                int t_N_case_het;
                int t_N_ctrl_hom;
                int t_N_ctrl_het;
                int t_N_case;
                int t_N_ctrl;


bool isTest;
      double vm, rss;
   //process_mem_usage(vm, rss);
   //std::cout << "VM: " << vm << "; RSS: " << rss << std::endl;
    for (int i = 0; i < m_to_test; i++) {
      if(m_isQuery){
              //std::cout << "mth_start " << mth_start << std::endl;
              //std::cout << "i " << i << std::endl;
              //std::cout << "mth_start+i " << mth_start+i << std::endl;
        markerIndex = m_markerIndicesToInclude[mth_start+i];
      }else{
              //std::cout << "m_isQuery " << m_isQuery << std::endl;
        markerIndex = 0;
      }
         //process_mem_usage(vm, rss);
	 //std::cout << "i " << i << std::endl; 
  // std::cout << "VM0: " << vm << "; RSS0: " << rss << std::endl;
    isTest = getScoreTest_SPA(markerIndex, traitType, t_chromosome, t_position, t_rsid, t_allele0, t_allele1, t_AF, t_AC, t_info, t_Ntest, t_Beta, t_se, t_Tstat, t_var1, t_var2, t_noSPApval, t_SPApval, t_isSPAConverge, t_AF_case, t_AF_ctrl, t_N_case_hom, t_N_case_het, t_N_ctrl_hom, t_N_ctrl_het, t_N_case, t_N_ctrl);
    //std::cout << "isTest " << isTest << std::endl;
   // std::cout << "t_noSPApval " << t_noSPApval << std::endl;
   //process_mem_usage(vm, rss);
   //std::cout << "VM1: " << vm << "; RSS1: " << rss << std::endl;

    if(isTest){
	    std::cout << "t_chromosome: " << t_chromosome << std::endl; 

	std::cout << "i " << i << std::endl;


        chromVec.at(i) = t_chromosome;
	std::cout << "chromVec.at(i) " << chromVec.at(i) << std::endl;
        posVec.at(i) = t_position;
        rsIDVec.at(i) = t_rsid;
        A1Vec.at(i) = t_allele0;
        A2Vec.at(i) = t_allele1;
        altFreqVec.at(i) = t_AF;
        ACVec.at(i) = t_AC;
        infoVec.at(i) = t_info;
        NVec.at(i) = t_Ntest;
        BetaVec.at(i) = t_Beta;
        seBetaVec.at(i) = t_se;
        TstatVec.at(i) = t_Tstat;
        var1Vec.at(i) = t_var1;
        var2Vec.at(i) = t_var2;
        pvalVec.at(i) = t_noSPApval;
        if(traitType == "binary"){
                SPApvalVec.at(i) = t_SPApval;
                SPAConverge.at(i) = t_isSPAConverge;
                if(ScoreTestObj.m_isOutputAFinCaseCtrl){
                        AFinCaseVec.at(i) = t_AF_case;
                        AFinCtrlVec.at(i) = t_AF_ctrl;
                }
                if(ScoreTestObj.m_isOutputHetHomCountsinCaseCtrl){
                        N_case_homVec.at(i) = t_N_case_hom;
                        N_case_hetVec.at(i) = t_N_case_het;
                        N_ctrl_homVec.at(i) = t_N_ctrl_hom;
                        N_ctrl_hetVec.at(i) = t_N_ctrl_het;
			                }
        }

      }else{
        chromVec.at(i) = "-1";
      }

    }

}
