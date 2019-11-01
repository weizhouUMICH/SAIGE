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
std::vector<float > dosagesforOneMarker;
std::vector<int > jIndexforOneMarker;
std::vector<int > iIndexforOneMarker;
std::vector<int > iIndexVec;
std::vector<int > jIndexVec;

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
    std::cout << "Number of samples in the vcf file: " << sample_size << endl;
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
  savvy::variant_group_iterator<savvy::compressed_vector<float>> it(marker_file, marker_group_line);
  savvy::variant_group_iterator<savvy::compressed_vector<float>> end{};
  group_matrix.resize(0);
  std::cout << "std::size_t sample_size = marker_file.samples().size();" << marker_file.samples().size() << std::endl;
  std::vector< int > indexforMissingAll;

  if (it != end){
    int missing_cnt = 0;
    std::vector< int > indexforMissing;
    int cnt = 0;
    float AC = 0;
    std::vector<std::string> markerIDs;
    std::vector<float> markerAFs;
    std::vector<float> MACs;
    markerIDs.clear();
    markerAFs.clear();
    MACs.clear();
    for ( ; it != end; ++it)
    {

      AC = 0;
      missing_cnt = 0;
      indexforMissing.clear();
      dosagesforOneMarker.clear();
      iIndexforOneMarker.clear();
      jIndexforOneMarker.clear();


      it.group_id();
      it.sites();
      //std::cout << "cnt: " << cnt << std::endl;	

      std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alt();
      //std::cout << it->prop("R2") << std::endl; 
      std::string markerInfo_str = it->prop("R2");
      float markerInfo = strtof((markerInfo_str).c_str(),0);
      for (auto dose_it = it->data().begin(); dose_it != it->data().end(); ++dose_it){
	int lengthi = std::distance(it->data().begin(), it->data().end());
	int i = dose_it.offset();	
	//std::cout << "i " << i << std::endl;
        if(genetest_sample_idx_vcfDosage[i] >= 0) {
	  //std::cout << "genetest_sample_idx_vcfDosage[i] " << genetest_sample_idx_vcfDosage[i] << std::endl;		
          if (std::isnan(*dose_it)) {

      //      dosagesforOneMarker[genetest_sample_idx_vcfDosage[i]] = float(1);
            ++missing_cnt;
            indexforMissing.push_back(genetest_sample_idx_vcfDosage[i]);
	    indexforMissingAll.push_back(genetest_sample_idx_vcfDosage[i]);
          }else {
	    if(*dose_it > 0){	
		dosagesforOneMarker.push_back(*dose_it);
		jIndexforOneMarker.push_back(cnt+1);
		iIndexforOneMarker.push_back(genetest_sample_idx_vcfDosage[i]+1);
            		//dosagesforOneMarker[genetest_sample_idx_vcfDosage[i]] = *dose_it;
            	AC = AC + *dose_it;
	    }	

          }
        } 
      //group_matrix[cnt * sample_size + dose_it.offset()] = *dose_it;
      }

      float AF = (float)(AC) / 2 / (float)(genetest_samplesize_vcfDosage - missing_cnt) ;
      //check if the AF of the marker is within the required range
      float MAF;
      float MAC;
      if(AF >= 0.5){
        MAF = 1 - AF;
        MAC = (float)(genetest_samplesize_vcfDosage - missing_cnt) *2 - AC;
      }else{
        MAF = AF;
	MAC = AC;
      }
      if(MAF >= minMAF && MAF <= maxMAF && markerInfo >= minInfo){
        if(missing_cnt > 0){
	  std::cout << "missing_cnt > 0!" << std::endl;
          float imputeDosage = 2*AF;
          for (unsigned int i = 0; i < indexforMissing.size(); i++){			
		dosagesforOneMarker.push_back(imputeDosage);
		jIndexforOneMarker.push_back(cnt+1);
		iIndexforOneMarker.push_back(indexforMissing[i]+1);
            //dosagesforOneMarker[indexforMissing[i]] = imputeDosage;
          }
       }

//	if(AF > 0.5){
//		std::cout << marker_id << " has AF > 0.5, so the alleles are flipped to use the dosages for minor allele";
//		dosagesforOneMarker = 2 - dosagesforOneMarker;
		//AF = AF - 1;
//	}
        group_matrix.insert(std::end(group_matrix), std::begin(dosagesforOneMarker), std::end(dosagesforOneMarker));
        iIndexVec.insert(std::end(iIndexVec), std::begin(iIndexforOneMarker), std::end(iIndexforOneMarker));
        jIndexVec.insert(std::end(jIndexVec), std::begin(jIndexforOneMarker), std::end(jIndexforOneMarker));
        cnt = cnt + 1;
        markerIDs.push_back(marker_id);
        markerAFs.push_back(AF);
        MACs.push_back(MAC);
        //std::cout << "MAF: " << MAF << " minMAF: " << minMAF << " maxMAF: " << maxMAF << std::endl;
	//std::cout << "marker_id: " << marker_id << std::endl;
      }
    } //for ( ; it != end; ++it)  
    //group_matrix.resize(sample_size * cnt);
    result["dosages"] = group_matrix;
    result["indexforMissing"] = indexforMissingAll;
    result["markerIDs"] = markerIDs;
    result["markerAFs"] = markerAFs;
    result["cnt"] = cnt;
    result["iIndex"] = iIndexVec;
    result["jIndex"] = jIndexVec;
  }else{
    //if there is it === end
    result["cnt"] = 0;
  }

  iIndexVec.clear();
  jIndexVec.clear();
  group_matrix.clear();

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
