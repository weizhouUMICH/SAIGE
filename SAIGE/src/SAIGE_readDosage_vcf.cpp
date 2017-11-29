#include "VcfRecord.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h> 
#include <stdlib.h>
#include <cstring>


// Open the vcf file for reading.
VcfFileReader reader;
VcfHeader header;
VcfRecord record;
std::string testField;

Rcpp::IntegerVector gm_sample_idx_vcfDosage;
int gmtest_samplesize_vcfDosage;

using namespace std;

// [[Rcpp::export]]
void SetSampleIdx_vcfDosage(Rcpp::IntegerVector sample_idx, int Ntest){
  gmtest_samplesize_vcfDosage = Ntest;
  gm_sample_idx_vcfDosage = sample_idx;

}

//[[Rcpp::export]]
void setTestField(std::string testFieldInput){
  testField = testFieldInput;
}



//const std::string& s1 = std::string();

// [[Rcpp::export]]
bool setgenoTest_vcfDosage(const std::string& vcfFileName,  const std::string& vcfFileIndex, const std::string& ids_to_exclude_vcf, const std::string& ids_to_include_vcf, const std::string& chromNam, int32_t start = 0, int32_t end = 0){

  bool isVcfOpen = reader.open(vcfFileName.c_str(), header);
  if(isVcfOpen){
    std::cout << "Open VCF done" << std::endl;
    std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << header.getNumMetaLines() << endl;
    std::cout << "Number of samples in in the vcf file: " << header.getNumSamples() << endl;
  }else{
    std::cout << "WARNING: Open VCF failed" << std::endl;
  }

  bool isreadVCfIndex = reader.readVcfIndex(vcfFileIndex.c_str());
  if(isreadVCfIndex){
    std::cout << "Read VCF index done" << std::endl;
  }else{
    isVcfOpen = false;
    std::cout << "WARNING: Read VCF index failed" << std::endl;
  }

  if(ids_to_exclude_vcf.length() > 0){
    bool isexclude = reader.setExcludeIDs(ids_to_exclude_vcf.c_str());
    if(isexclude){
      std::cout << "Excluding variants done" << std::endl;
    }else{
      isVcfOpen = false;
      std::cout << "WARNING: Excluding variants failed" << std::endl;
    }
  }

  if(ids_to_include_vcf.length() > 0){
    bool isinclude = reader.setIncludeIDs(ids_to_include_vcf.c_str());
    if(isinclude){
      std::cout << "Including variants done" << std::endl;
    }else{
      isVcfOpen = false;
      std::cout << "WARNING: Including variants failed" << std::endl;
    }
  }
  
//  std::cout << "chromNam.length() " << chromNam.length() << std::endl; 
//  std::cout << "chromNam " << chromNam << std::endl; 
  //if(chromNam.length() > 0){
  std::string myZeroString ("0");

  if(chromNam  != myZeroString){
    if(end == 0){
      bool isReadSecion = reader.setReadSection(chromNam.c_str());
      if(isReadSecion){
        std::cout << "set to read chromosome " << chromNam << std::endl;
      }else{
	isVcfOpen = false;	
        std::cout << "WARNING: set to read chromosome " << chromNam << " failed" << std::endl;
      }
    }else{
      bool isReadBasedSecion = reader.set1BasedReadSection(chromNam.c_str(), start, end);
      if(isReadBasedSecion){
        std::cout << "set to read chromosome:start-end (nooverlap): " << chromNam << ":" << start << "-" << end << std::endl;
      }else{
        isVcfOpen = false;
        std::cout << "WARNING: set to read chromosome:start-end (nooverlap): " << chromNam << ":" << start << "-" << end  << " failed" << std::endl;
      }
    }
  }else{
    std::cout << "Set to read the entire vcf " << std::endl;
  }
  return(isVcfOpen);
}


// [[Rcpp::export]]
int getNumofSamples(){
  return(header.getNumSamples());
}


// [[Rcpp::export]]
std::vector< std::string > getSampleIDlist(){
  unsigned int numOfSamples = header.getNumSamples();
  std::vector< std::string > sampleIDList (numOfSamples);
  const char* sampleIDpre;
  for( std::size_t i = 0; i < numOfSamples; ++i ) {
    sampleIDpre =  header.getSampleName(i);

    std::string sampleID(sampleIDpre);
    sampleIDList.push_back(sampleIDpre);
  }
  return(sampleIDList);
}



// [[Rcpp::export]]
bool getGenoOfnthVar_vcfDosage_pre(){
  bool isReadVariant = reader.readRecord(record);
  return(isReadVariant);
}

//bool isGetDosage = FALSE;

/*
// [[Rcpp::export]]
bool getIsGetDosage(){
  return(isGetDosage);
}
*/

// [[Rcpp::export]]
Rcpp::List getGenoOfnthVar_vcfDosage(int mth) {
  //bool isGetDosage = TRUE;  
  using namespace Rcpp;
  std::string rsid(record.getIDStr());
  std::string chromosome(record.getChromStr());
  //int positionnum = record.get1BasedPosition();
  //std::string position = std::to_string(positionnum);
  int position = record.get1BasedPosition();
  int numAlt = record.getNumAlts();
  std::vector< std::string > alleles(numAlt);
  std::string alleles0(record.getRefStr());
  std::string alleles1(record.getAltStr()); 
  List result ;


//  if(alleles0 == "." || alleles1 == "."){
//    isGetDosage = FALSE;
//  }else if(numAlt != 1){
//    isGetDosage = FALSE;  
//  }else{
    VcfRecordGenotype & genoInfo = record.getGenotypeInfo();
    int numSamples = record.getNumSamples();
    std::vector< float > dosages;
    float dosage;
    const std::string * geno;
    float AC = 0;
   int genoi;
   int numGTs;
   //////////
   dosages.reserve(gmtest_samplesize_vcfDosage);
   dosages.resize(gmtest_samplesize_vcfDosage);
   for (uint i = 0; i < numSamples; i++) {
     if(testField == "DS"){  
       geno = genoInfo.getString(testField.c_str(),i);
       if(strcmp((*geno).c_str(),".")>0){
         dosage = atof((*geno).c_str());
       }else{
         dosage = -1; //missing dosages
       }
     }else if(testField == "GT"){
       numGTs = record.getNumGTs(i);
       genoi = 0;
       for(int j = 0; j < numGTs; j++)
       {
	 genoi = genoi + genoInfo.getGT(i,j);
       }
       dosage = genoi;
     }
 
     if(gm_sample_idx_vcfDosage[i] >= 0){
       dosages[gm_sample_idx_vcfDosage[i]] = dosage;
       AC = AC + dosage;
     }
  }
    float AF = AC/ 2/ ((float)gmtest_samplesize_vcfDosage) ;

    result[ "dosages" ] = dosages ;  
    dosages.clear();
  DataFrame variants = DataFrame::create(
                Named("chromosome") = chromosome,
                Named("position") = position,
                Named("rsid") = rsid,
        //        Named("number_of_alleles") = number_of_allele,
                Named("allele0") = alleles0,
                Named("allele1") = alleles1,
                _["stringsAsFactors"] = false,
		Named("AC") = AC,
                Named("AF") = AF
        ) ;


  //List result ;
  result[ "variants" ] = variants ;
  //result["isGetDosage"] = isGetDosage;
  //result[ "dosages" ] = dosages ;

//}
  
  return( result ) ;
}

