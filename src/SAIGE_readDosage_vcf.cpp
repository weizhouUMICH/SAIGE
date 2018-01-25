#include "../thirdParty/cget/include/savvy/reader.hpp"
#include "../thirdParty/cget/include/savvy/varint.hpp"
#include "../thirdParty/cget/include/savvy/sav_reader.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <stdlib.h>
#include <cstring>


// Open the vcf file for reading.
savvy::indexed_reader reader;
//VcfHeader header;
savvy::variant<std::vector<float>> record;
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



// [[Rcpp::export]]
bool setgenoTest_vcfDosage(const std::string& vcfFileName,  const std::string& vcfFileIndex, const std::string& vcfField, const std::string& ids_to_exclude_vcf, const std::string& ids_to_include_vcf, const std::string& chromNam, int32_t start = 0, int32_t end = 0){

  if(vcfField == "DS"){
    reader = savvy::indexed_reader(vcfFileName, {chromNam, std::uint32_t(start), std::uint32_t(end)}, savvy::fmt::dosage);
  }else if(vcfField == "GT"){
    reader = savvy::indexed_reader(vcfFileName, {chromNam, std::uint32_t(start), std::uint32_t(end)}, savvy::fmt::genotype);
  }

  bool isVcfOpen = reader.good();
  if(isVcfOpen){
    std::cout << "Open VCF done" << std::endl;
    cout << "To read the field " << vcfField << endl;
    std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << reader.headers().size() << endl;
    std::cout << "Number of samples in in the vcf file: " << reader.samples().size() << endl;
  }else{
    std::cout << "WARNING: Open VCF failed" << std::endl;
  }

  return(isVcfOpen);
}


// [[Rcpp::export]]
int getNumofSamples(){
  return(reader.samples().size());
}


// [[Rcpp::export]]
std::vector< std::string > getSampleIDlist(){
  std::vector< std::string > sampleIDList (reader.samples().begin(), reader.samples().end());
  return(sampleIDList);
}



// [[Rcpp::export]]
bool getGenoOfnthVar_vcfDosage_pre(){
  bool isReadVariant = reader >> record;
  return(isReadVariant);
}


// [[Rcpp::export]]
Rcpp::List getGenoOfnthVar_vcfDosage(int mth) {
  //bool isGetDosage = TRUE;
  using namespace Rcpp;
  std::string rsid(record.prop("ID"));
  std::string chromosome(record.chromosome());
  //int positionnum = record.get1BasedPosition();
  //std::string position = std::to_string(positionnum);
  int position = record.position();
  int numAlt = 1;
  std::vector< std::string > alleles(numAlt);
  std::string alleles0(record.ref());
  std::string alleles1(record.alt());
  List result ;


  int numSamples = reader.samples().size();
  std::vector< float > dosages;
  float dosage;
  const std::string * geno;
  float AC = 0;
  int genoi;
  int numGTs;
  //////////
  //dosages.resize(record.data().size());
  dosages.resize(gmtest_samplesize_vcfDosage);

  std::size_t missing_cnt = 0;
  std::vector< int > indexforMissing;
  for (std::vector<float>::iterator it = record.data().begin(); it != record.data().end(); ++it) {
    int i = std::distance(record.data().begin(), it);
    if(gm_sample_idx_vcfDosage[i] >= 0) {
      if (std::isnan(*it)) {
        dosages[gm_sample_idx_vcfDosage[i]] = -1;
        ++missing_cnt;
        indexforMissing.push_back(gm_sample_idx_vcfDosage[i]);
      } else {
        dosages[gm_sample_idx_vcfDosage[i]] = *it;
        AC = AC + *it;
      }
    }
  }

  float AF = AC / 2 / (float)(gmtest_samplesize_vcfDosage - missing_cnt) ;
  //set missing dosages to be the 2*AF
  if(missing_cnt > 0){
    float imputeDosage = 2*AF;
    for (unsigned int i = 0; i < indexforMissing.size(); i++)
    {
      dosages[indexforMissing[i]] = imputeDosage;
    }  
  }

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
  );


  //List result ;
  result[ "variants" ] = variants ;
  //result["isGetDosage"] = isGetDosage;
  //result[ "dosages" ] = dosages ;
  indexforMissing.clear();
  
//}
  

  return( result ) ;
}
