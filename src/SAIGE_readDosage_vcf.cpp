#include "savvy/reader.hpp"
#include "savvy/varint.hpp"
#include "savvy/sav_reader.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <stdlib.h>
#include <cstring>
#include <limits>

// Open the vcf file for reading.
savvy::reader reader{""};

//VcfHeader header;
savvy::variant record;
std::string testField;

Rcpp::IntegerVector gm_sample_idx_vcfDosage;
int gmtest_samplesize_vcfDosage;

int numSamples_vcf;
bool isDropMissingDosages_vcf = false;

using namespace std;




// [[Rcpp::export]]
void SetSampleIdx_vcfDosage(Rcpp::IntegerVector sample_idx, int Ntest){
  gmtest_samplesize_vcfDosage = Ntest;
  gm_sample_idx_vcfDosage = sample_idx;
}

//TODO: Why does this function exist? testField is set in setgenoTest_vcfDosage()
//[[Rcpp::export]]
void setTestField(std::string testFieldInput){
  testField = testFieldInput;
  if (std::find_if(reader.format_headers().begin(), reader.format_headers().end(),
    [](const savvy::header_value_details& h) { return h.id == testField; }) == reader.format_headers().end()) {
    if ((testField == "DS" || testField == "GT") && std::find_if(reader.format_headers().begin(), reader.format_headers().end(),
      [](const savvy::header_value_details& h) { return h.id == "HDS"; }) != reader.format_headers().end()) {
      testField = "HDS";
    }
  }
}

// [[Rcpp::export]]
void setIsDropMissingDosages_vcf (bool isdropmissingdosages){
  isDropMissingDosages_vcf = isdropmissingdosages;
}

// [[Rcpp::export]]
bool setgenoTest_vcfDosage(const std::string& vcfFileName,  const std::string& vcfFileIndex, const std::string& vcfField, const std::string& ids_to_exclude_vcf, const std::string& ids_to_include_vcf, const std::string& chromNam, int32_t start = 0, int32_t end = 0){
  testField = vcfField;  
  reader = savvy::reader(vcfFileName);
  
  if(chromNam != "" || start > 1 || std::uint32_t(end) < std::numeric_limits<std::int32_t>::max()){
    std::cout << "Setting genomic region " << chromNam << ":" << std::uint32_t(start) << ":" << std::uint32_t(end) << std::endl;
    reader.reset_bounds({chromNam, std::uint32_t(start), std::uint32_t(end)});
  }

  bool isVcfOpen = reader.good();
  if(isVcfOpen){
    std::cout << "Open VCF done" << std::endl;
    cout << "To read the field " << vcfField << endl;
    std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << reader.headers().size() << endl;
    std::cout << "Number of samples in in the vcf file: " << reader.samples().size() << endl;
    numSamples_vcf = reader.samples().size();

    if (std::find_if(reader.format_headers().begin(), reader.format_headers().end(),
      [](const savvy::header_value_details& h) { return h.id == testField; }) == reader.format_headers().end()) {
      if ((testField == "DS" || testField == "GT") && std::find_if(reader.format_headers().begin(), reader.format_headers().end(),
        [](const savvy::header_value_details& h) { return h.id == "HDS"; }) != reader.format_headers().end()) {
        testField = "HDS";
      } else {
        std::cerr << "ERROR: vcfField (" << testField << ") not present in genotype file." << std::endl;
        return false;
      }
    }
  } else {
    numSamples_vcf = -10;
    std::cerr << "WARNING: Open VCF failed" << std::endl;
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
  if (isReadVariant) {
    if (record.alts().size() != 1) { // TODO: support multiallelics
      std::cerr << "Error: variants must be biallelic" << std::endl;
      isReadVariant = false;
    }
  } else if (reader.bad()) {
    std::cerr << "Error: READ FAILED" << std::endl;
  }
  return(isReadVariant);
}


// [[Rcpp::export]]
Rcpp::List getGenoOfnthVar_vcfDosage(int mth) {
  //bool isGetDosage = TRUE;
  using namespace Rcpp;
  std::string rsid(record.id());
  std::string chromosome(record.chromosome());
  //int positionnum = record.get1BasedPosition();
  //std::string position = std::to_string(positionnum);
  int position = record.position();
  int numAlt = 1;
  std::vector< std::string > alleles(numAlt);
  std::string alleles0(record.ref());
  std::string alleles1(record.alts().empty() ? "" : record.alts()[0]);
  float variant_r2 = 1.f;
  record.get_info("R2", variant_r2);
  List result ;


  int numSamples = reader.samples().size();
  std::vector< float > dosages;
  savvy::compressed_vector<float> sparse_dosages;
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

  if (!record.get_format(testField, sparse_dosages)) {
    std::cerr << "Warning: missing FMT field (" << testField << ")" << std::endl;
    return result;
  }
  if (!sparse_dosages.size() || !reader.samples().size()) {
    std::cerr << "Warning: no sample level data available" << std::endl;
    return result;
  }

  std::size_t ploidy = sparse_dosages.size() / reader.samples().size();
  savvy::stride_reduce(sparse_dosages, ploidy);

  for (auto it = sparse_dosages.begin(); it != sparse_dosages.end(); ++it) {
    int i = it.offset();
    if(gm_sample_idx_vcfDosage[i] >= 0) {
      if (std::isnan(*it)) {
        dosages[gm_sample_idx_vcfDosage[i]] = -1.f;
        ++missing_cnt;
        indexforMissing.push_back(gm_sample_idx_vcfDosage[i]);
      } else {
        dosages[gm_sample_idx_vcfDosage[i]] = *it;
        AC = AC + *it;
      }
    }
  }
  float AF = 0.f;
  if(gmtest_samplesize_vcfDosage > missing_cnt){
    AF = AC / 2.f / (float)(gmtest_samplesize_vcfDosage - missing_cnt) ;
  }
  //set missing dosages to be the 2*AF
  if(missing_cnt > 0){
    float imputeDosage = 2*AF;
    for (unsigned int i = 0; i < indexforMissing.size(); i++)
    {
      dosages[indexforMissing[i]] = imputeDosage;
    } 
    if(!isDropMissingDosages_vcf){
      AC = AC + imputeDosage * missing_cnt;
    } 
  }

    result[ "dosages" ] = dosages;
    result["indexforMissing"] = indexforMissing;
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
    Named("AF") = AF,
    Named("markerInfo") = variant_r2
  );


  //List result ;
  result[ "variants" ] = variants ;
  //result["isGetDosage"] = isGetDosage;
  //result[ "dosages" ] = dosages ;
  indexforMissing.clear();
  
//}
  

  return( result ) ;
}


// [[Rcpp::export]]
void closetestGenoFile_vcfDosage() //needs further check
{
//  gm_sample_idx_vcfDosage.clear();
  printf("closed the genofile!\n");

}

// [[Rcpp::export]]
int getSampleSizeinVCF(){
  return(numSamples_vcf);
}
