#include "../thirdParty/cget/include/savvy/reader.hpp"
#include "../thirdParty/cget/include/savvy/region.hpp"
#include "../thirdParty/cget/include/savvy/variant_group_iterator.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <stdlib.h>
#include <cstring>

using namespace std;
//This file is revised from the Variant-Sample Matrix file from the savvy library by Jonathon LeFaive
//https://github.com/statgen/savvy/blob/develop/doc/variant_group_iterator.md

//std::ifstream marker_group_file("file_groups.txt", std::ios::binary);
//savvy::indexed_reader marker_file("marker_file.sav", {""}, savvy::fmt::dosage);
savvy::indexed_reader marker_file;
//std::string marker_group_line;
//std::size_t sample_size = marker_file.samples().size();
std::size_t sample_size;
std::vector<float> group_matrix;
std::vector< float > dosagesforOneMarker;

std::vector<int> genetest_sample_idx_vcfDosage;
int genetest_samplesize_vcfDosage;

float minMAF;
float maxMAF;


// [[Rcpp::export]]
void SetSampleIdx_forGenetest_vcfDosage(std::vector<int> sample_idx, int Ntest){
  genetest_samplesize_vcfDosage = Ntest;
  genetest_sample_idx_vcfDosage = sample_idx;
}

// [[Rcpp::export]]
void setMAFcutoffs(float minVal, float maxVal){
  minMAF = minVal;
  maxMAF = maxVal;
}

// [[Rcpp::export]]
bool setvcfDosageMatrix(const std::string& vcfFileName,  const std::string& vcfFileIndex, const std::string& vcfField){

  if(vcfField == "DS"){
    marker_file = savvy::indexed_reader(vcfFileName, {""}, savvy::fmt::ds);
  }else if(vcfField == "GT"){
    marker_file = savvy::indexed_reader(vcfFileName, {""}, savvy::fmt::ac);
  }
  bool isVcfOpen = marker_file.good();

  if(isVcfOpen){
    std::cout << "Open VCF done" << std::endl;
    cout << "To read the field " << vcfField << endl;
    std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << marker_file.headers().size() << endl;
    sample_size = marker_file.samples().size(); 
    std::cout << "Number of samples in in the vcf file: " << sample_size << endl;
  }else{
    std::cout << "WARNING: Open VCF failed" << std::endl;
  }
  return(isVcfOpen);
}

// [[Rcpp::export]]
Rcpp::List getGenoOfGene_vcf(std::string marker_group_line, float minInfo) {
  //bool isGetDosage = TRUE;
  using namespace Rcpp;
  List result ;
//  std::cout << "here0!" << std::endl; 
  //savvy::variant_group_iterator<savvy::compressed_vector<float>> it(marker_file, marker_group_line);
  savvy::variant_group_iterator<savvy::compressed_vector<float>> it(marker_file, marker_group_line);
//  std::cout << "here1!" << std::endl;
  savvy::variant_group_iterator<savvy::compressed_vector<float>> end{};
//  std::cout << "it.sites().size(): " << it.sites().size() << std::endl; 
//  std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alt();
//  std::cout << marker_id << std::endl;
  group_matrix.resize(0);
  std::cout << "std::size_t sample_size = marker_file.samples().size();" << marker_file.samples().size() << std::endl;

  if (it != end){
    //group_matrix.resize(it.sites().size() * sample_size);
    //group_matrix.reserve(it.sites().size() * genetest_samplesize_vcfDosage);
    std::size_t missing_cnt = 0;
    std::vector< int > indexforMissing;
    std::size_t cnt = 0;
    float AC = 0;
    std::vector<std::string> markerIDs;
    std::vector<float> markerAFs;
    markerIDs.clear();
    markerAFs.clear();

//  std::cout << "here2!" << std::endl; 
    for ( ; it != end; ++it)
    {
      AC = 0;
      missing_cnt = 0;
      indexforMissing.clear();
      dosagesforOneMarker.clear();
      dosagesforOneMarker.resize(genetest_samplesize_vcfDosage);


//  std::cout << "here3!" << std::endl; 
      it.group_id();
      it.sites();
 //     std::cout << "it.group_id(): " << it.group_id() << std::endl;	
 //     std::cout << "it.sites(): " << it.sites() << std::endl;	

      std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alt();
      //std::cout << "here4!" << std::endl; 
      //std::cout << it->prop("R2") << std::endl; 
      std::string markerInfo_str = it->prop("R2");
      float markerInfo = strtof((markerInfo_str).c_str(),0);
      for (auto dose_it = it->data().begin(); dose_it != it->data().end(); ++dose_it){
	int lengthi = std::distance(it->data().begin(), it->data().end());
	//std::cout << "lengthi: " << lengthi << std::endl;
        //std::cout << "dose_it.offset() " << dose_it.offset() << std::endl;
	//int i = std::distance(it->data().begin(), dose_it);
	//std::cout << "i " << i << std::endl;
	int i = dose_it.offset();	
        if(genetest_sample_idx_vcfDosage[i] >= 0) {
          if (std::isnan(*dose_it)) {
            dosagesforOneMarker[genetest_sample_idx_vcfDosage[i]] = -1;
            ++missing_cnt;
            indexforMissing.push_back(genetest_sample_idx_vcfDosage[i]);
          }else {
	      //if(i < 1000){
	      //std::cout << "*dose_it " << *dose_it << std::endl;	
	      //std::cout << "i " << i << std::endl;	
	    //}	
            dosagesforOneMarker[genetest_sample_idx_vcfDosage[i]] = *dose_it;
            AC = AC + *dose_it;
          }
        } 
      //group_matrix[cnt * sample_size + dose_it.offset()] = *dose_it;
      }

      float AF = (float)(AC) / 2 / (float)(genetest_samplesize_vcfDosage - missing_cnt) ;
      //check if the AF of the marker is within the required range
      float MAF;
      if(AF >= 0.5){
        MAF = 1 - AF;
      }else{
        MAF = AF;
      }
      if(MAF >= minMAF && MAF <= maxMAF && markerInfo >= minInfo){
        if(missing_cnt > 0){
	  std::cout << "missing_cnt > 0!" << std::endl;
          float imputeDosage = 2*AF;
          for (unsigned int i = 0; i < indexforMissing.size(); i++){
            dosagesforOneMarker[indexforMissing[i]] = imputeDosage;
          }
       }
        group_matrix.insert(std::end(group_matrix), std::begin(dosagesforOneMarker), std::end(dosagesforOneMarker));
        cnt = cnt + 1;
        markerIDs.push_back(marker_id);
        markerAFs.push_back(AF);
        //std::cout << "MAF: " << MAF << " minMAF: " << minMAF << " maxMAF: " << maxMAF << std::endl;
	//std::cout << "marker_id: " << marker_id << std::endl;
      }
    }  
    //group_matrix.resize(sample_size * cnt);
    result["dosages"] = group_matrix;
    result["markerIDs"] = markerIDs;
    result["markerAFs"] = markerAFs;
    result["cnt"] = cnt;
  }else{
    //if there is it === end
    result["cnt"] = 0;
  }

  return(result);
}


// [[Rcpp::export]]
void closevcfDosageFile() //needs further check
{
  group_matrix.clear();
  dosagesforOneMarker.clear();

  genetest_sample_idx_vcfDosage.clear();

//  gm_sample_idx_vcfDosage.clear();
//  printf("closed the genofile!\n");

}
