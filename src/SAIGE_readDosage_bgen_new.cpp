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

#include <Rcpp.h>


//A global variable for the dosage file to test
FILE *m_fin;
unsigned char * m_buf;
unsigned char * m_zBuf;
uint m_zBufLens;
uint m_bufLens;



int gmtest_samplesize;
Rcpp::IntegerVector gm_sample_idx;

//genoToTest_bgenDosage
bool m_isQuery;
std::vector< int > m_markerIndicesToInclude;
bool isReadVariantBgen = true;
//bool isOutputHetHomCountsinCaseCtrl = false;
bool isDropMissingDosages_bgen = false;
//double bgenMinMAF = 0;
//double bgenMinINFO = 0;
uint numSamples_bgen;
bool isSparseDosage_bgen = false;

// [[Rcpp::export]]
void setIsDropMissingDosages_bgen (bool isdropmissingdosages){
  isDropMissingDosages_bgen = isdropmissingdosages;
}

// [[Rcpp::export]]
void setIsSparseDosage_bgen (bool isSparseDosage){
   isSparseDosage_bgen = isSparseDosage;
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
    }else if(
     m_markerIndicesToInclude.size() == 0){	    
     std::cout << "no query list is provided" << std::endl;
     m_isQuery = false;
    }else{
     m_isQuery = true;
     mtotest = m_markerIndicesToInclude.size();
    }

    /****code from BOLT-LMM v2.3.4***/
    mkl_set_num_threads(1);

    /********** READ HEADER v1.2**********/
    m_fin = fopen(t_bgenFileName.c_str(), "rb");
    uint offset; fread(&offset, 4, 1, m_fin); //cout << "offset: " << offset << endl;
    uint L_H; fread(&L_H, 4, 1, m_fin); //cout << "L_H: " << L_H << endl;
    uint m_M0; fread(&m_M0, 4, 1, m_fin); cout << "snpBlocks (Mbgen): " << m_M0 << endl;
    assert(Mbgen != 0);
    if(!m_isQuery){mtotest = m_M0};
    uint Nbgen; fread(&m_N0, 4, 1, m_fin); cout << "samples (Nbgen): " << m_N0 << endl;
    uint m_Nsample = t_SampleInBgen.size();
    if (Nbgen != m_Nsample) {
      cerr << "ERROR: Number of samples in BGEN header does not match sample file" << endl;
      exit(1);
    }
    char magic[5]; fread(magic, 1, 4, m_fin); magic[4] = '\0'; //cout << "magic bytes: " << string(magic) << endl;
    fseek(m_fin, L_H-20, SEEK_CUR); //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
    uint flags; fread(&flags, 4, 1, m_fin); //cout << "flags: " << flags << endl;
    uint CompressedSNPBlocks = flags&3; cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << endl;
    assert(CompressedSNPBlocks==1); // REQUIRE CompressedSNPBlocks==1
    uint Layout = (flags>>2)&0xf; cout << "Layout: " << Layout << endl;
    assert(Layout==1 || Layout==2); // REQUIRE Layout==1 or Layout==2
    fseek(m_fin, offset+4, SEEK_SET);
    return(m_M0);
  }





void Parse2(uchar *buf, uint bufLen, const uchar *zBuf, uint zBufLen,std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<unsigned int> & indexforMissing, double & info, std::vector< uint > & iIndex) {

    uLongf destLen = bufLen;
    if (uncompress(buf, &destLen, zBuf, zBufLen) != Z_OK || destLen != bufLen) {
      cerr << "ERROR: uncompress() failed" << endl;
      exit(1);
    }
    uchar *bufAt = buf;
    uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
    if (N != Nbgen) {
      cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << endl;
      exit(1);
    }
    uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
    if (K != 2U) {
      cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << endl;
      exit(1);
    }
    uint Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
      cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << endl;
      exit(1);
    }
    uint Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
      cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << endl;
      exit(1);
    }
    const uchar *ploidyMissBytes = bufAt;
    for (uint i = 0; i < N; i++) {
      uint ploidyMiss = *bufAt; bufAt++;
      if (ploidyMiss != 2U && ploidyMiss != 130U) {
        cerr << "ERROR: " << snpName << " has ploidy/missingness byte = " << ploidyMiss
             << " (not 2 or 130)" << endl;
        exit(1);
      }
    }
    uint Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
      cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << endl;
      exit(1);
    }
    uint B = *bufAt; bufAt++;
    if (B != 8U) {
      cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << endl;
      exit(1);
    }


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
     if (ploidyMissBytes[i] == 130U){
      bufAt += 2;
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      p00 = 1 - p11 - p10;
      dosage = 2*p11 + p10;

        eij = dosage;
        fij = 4*p11 + p10;
        sum_eij += eij;
        sum_fij_minus_eij2 += fij - eij*eij;
        if(gm_sample_idx[i] >= 0){
          if(!m_isSparseDosages){
              dosages[gm_sample_idx[i]] = 2 - dosage;
          }else{
              if(2 - dosage > 0){
                dosages.push_back(2 - dosage);
                iIndex.push_back(gm_sample_idx[i]+1);
              }
          }
          sum_eij_sub += eij;
        }
     }else if(ploidyMissBytes[i] == 2U){
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

    AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;

    double thetaHat = sum_eij / (2* (m_N - missing_cnt));
    info = thetaHat==0 || thetaHat==1 ? 1 :
    1 - sum_fij_minus_eij2 / (2*(m_N - missing_cnt)*thetaHat*(1-thetaHat));

    if(missing_cnt > 0){
      std::cout << "sample index with missing dosages for snpName " << snpName << " :";
      if(m_N == missing_cnt){
        AF = 0;
      }else{
        AF = AC/ 2/ ((double) (m_N - missing_cnt)) ;
      }

      if(!m_isDropMissingDosages){
        double imputeDosage = 2*AF;
        for (unsigned int i = 0; i < indexforMissing.size(); i++)
        {
          if(!m_isSparseDosages){
            dosages[indexforMissing[i]] = imputeDosage;
          }else{
            dosages.push_back(imputeDosage);
            iIndex.push_back(indexforMissing[i]+1);
          }
            //std::cout << indexforMissing[i]+1 << ",";
          AC = AC + imputeDosage;
        }
      }
    }

    //return(info);
  }


// [[Rcpp::export]]
Rcpp::List getOneMarker(int t_fileStartPos)
  {


    if(m_isQuery){fseek(m_fin, t_fileStartPos, SEEK_SET);}
     std::string SNPID, RSID, chromosome, first_allele,second_allele ;
     uint position;
     std::vector< std::string > alleles ;
     std::vector< double > dosages;
     double AC, AF,
     std::vector< int > indexforMissing;
  std::vector< uint > iIndex;


    char snpID[65536], rsID[65536], chrStr[65536];
    char *allele1, *allele0;
    ushort LS; size_t numBoolRead = fread(&LS, 2, 1, m_fin); // cout << "LS: " << LS << " " << std::flush;
    bool isBoolRead;
    if ( numBoolRead > 0 ) {
      isBoolRead = true;
      fread(snpID, 1, LS, m_fin); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
      ushort LR; fread(&LR, 2, 1, m_fin); // cout << "LR: " << LR << " " << std::flush;
      fread(rsID, 1, LR, fin); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;
      RSID = string(rsID)=="." ? snpID : rsID;
      //std::string SNPID = string(snpID);

      ushort LC; fread(&LC, 2, 1, m_fin); // cout << "LC: " << LC << " " << std::flush;
      fread(chrStr, 1, LC, m_fin); chrStr[LC] = '\0';
      chromosome  = string(chrStr);

      uint physpos; fread(&physpos, 4, 1, m_fin); // cout << "physpos: " << physpos << " " << std::flush;
      position = physpos;
      ushort K; fread(&K, 2, 1, m_fin); //cout << "K: " << K << endl;
      if (K != 2) {
        cerr << "ERROR: Non-bi-allelic variant found: " << K << " alleles" << endl;
        exit(1);
      }

      uint LA; fread(&LA, 4, 1, m_fin); // cout << "LA: " << LA << " " << std::flush;
      if (LA > maxLA) {
        maxLA = 2*LA;
        free(allele1);
        allele1 = (char *) malloc(maxLA+1);
      }
      fread(allele1, 1, LA, m_fin); allele1[LA] = '\0';
      second_allele = string(allele1);
      uint LB; fread(&LB, 4, 1, fin); // cout << "LB: " << LB << " " << std::flush;
      if (LB > maxLB) {
        maxLB = 2*LB;
        free(allele0);
        allele0 = (char *) malloc(maxLB+1);
      }
      fread(allele0, 1, LB, m_fin); allele0[LB] = '\0';
      first_allele = string(allele0);

      uint C; fread(&C, 4, 1, m_fin); //cout << "C: " << C << endl;
      if (C > m_zBuf.size()) m_zBuf.resize(C-4);
      uint D; fread(&D, 4, 1, m_fin); //cout << "D: " << D << endl;
      m_zBufLens = C-4; m_bufLens = D;
      fread(&m_zBuf[0], 1, C-4, m_fin);
      AC = 0;
      AF = 0;
      info = 0;
      if (m_bufLens > m_buf.size()) m_buf.resize(m_bufLens);
           Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens, RSID, m_N0, dosages, AC, AF, indexforMissing, info, iIndex);
      }
        DataFrame variants = DataFrame::create(
                Named("chromosome") = chromosome,
                Named("position") = position,
                Named("rsid") = RSID,
        //        Named("number_of_alleles") = number_of_allele,
                Named("allele0") = first_allele,
                Named("allele1") = second_allele,
                _["stringsAsFactors"] = false,
                Named("AC") = AC,
                Named("AF") = AF,
                Named("info") = info,
                //Named("homN_cases") = homN_cases,
                //Named("hetN_cases") = hetN_cases,
                //Named("homN_ctrls") = homN_ctrls,
                //Named("hetN_ctrls") = hetN_ctrls
        );

    }else{
      isBoolRead = false;
      DataFrame variants = NULL;
    }
    List result ;
    result[ "variants" ] = variants ;
    result[ "dosages" ] = dosages ;
    result["iIndex"] = iIndex;
    result["indexforMissing"] = indexforMissing;
    result["isBoolRead"] = isBoolRead;

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



// [[Rcpp::export]]
double getMarkerInfo()
{
  return(markerInfo);
}

/**************************************
        By SLEE 09/06/17
*************************************/

// [[Rcpp::export]]
void SetSampleIdx(Rcpp::IntegerVector sample_idx, int Ntest){
	gmtest_samplesize = Ntest;
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
