#include "savvy/reader.hpp"
#include "savvy/region.hpp"
#include "variant_group_iterator.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <stdlib.h>
#include <cstring>
#include <limits>

using namespace std;
//This file is revised from the Variant-Sample Matrix file from the savvy library by Jonathon LeFaive
//https://github.com/statgen/savvy/blob/develop/doc/variant_group_iterator.md

//std::ifstream marker_group_file("file_groups.txt", std::ios::binary);
//savvy::indexed_reader marker_file("marker_file.sav", {""}, savvy::fmt::dosage);
savvy::reader marker_file{""};
std::string fmtField;
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

  marker_file = savvy::reader(vcfFileName);
  bool isVcfOpen = marker_file.good();

  if(isVcfOpen){
    fmtField = vcfField;
    std::cout << "Open VCF done" << std::endl;
    cout << "To read the field " << vcfField << endl;
    std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << marker_file.headers().size() << endl;
    sample_size = marker_file.samples().size(); 
    std::cout << "Number of samples in the vcf file: " << sample_size << endl;
    
    if (std::find_if(marker_file.format_headers().begin(), marker_file.format_headers().end(),
      [](const savvy::header_value_details& h) { return h.id == fmtField; }) == marker_file.format_headers().end()) {
      if ((fmtField == "DS" || fmtField == "GT") && std::find_if(marker_file.format_headers().begin(), marker_file.format_headers().end(),
        [](const savvy::header_value_details& h) { return h.id == "HDS"; }) != marker_file.format_headers().end()) {
        fmtField = "HDS";
      } else {
        std::cerr << "ERROR: vcfField (" << fmtField << ") not present in genotype file." << std::endl;
        return false;
      }
    }    
  } else {
    std::cerr << "WARNING: Open VCF failed" << std::endl;
  }
  return(isVcfOpen);
}

// [[Rcpp::export]]
int getNumofSamples_Matrix(){
  return(marker_file.samples().size());
}

// [[Rcpp::export]]
std::vector< std::string > getSampleIDlist_vcfMatrix(){
  std::vector< std::string > sampleIDList (marker_file.samples().begin(), marker_file.samples().end());
  return(sampleIDList);
}



// [[Rcpp::export]]
Rcpp::List getGenoOfGene_vcf(std::string marker_group_line, float minInfo) {
  //bool isGetDosage = TRUE;
  using namespace Rcpp;
  List result ;
  variant_group_iterator it(marker_file, marker_group_line);
  variant_group_iterator end{};
  savvy::compressed_vector<float> variant_dosages;
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
    std::vector<int> positions;
    std::vector<float> MACs;
    markerIDs.clear();
    markerAFs.clear();
    positions.clear();
    MACs.clear();
    for ( ; it != end; ++it)
    {
      if (it->alts().size() != 1)
      {
        std::cerr << "Warning: skipping multiallelic variant" << std::endl;
        continue;
      }

      AC = 0;
      missing_cnt = 0;
      indexforMissing.clear();
      dosagesforOneMarker.clear();
      iIndexforOneMarker.clear();
      jIndexforOneMarker.clear();


      it.group_id();
      it.sites();
      //std::cout << "cnt: " << cnt << std::endl;	

      std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alts()[0];
      //std::cout << it->prop("R2") << std::endl; 
      //std::cout << "marker_id: " << marker_id  << std::endl; 
      float markerInfo = 1.f;
      it->get_info("R2", markerInfo);

      it->get_format(fmtField, variant_dosages);
      std::size_t ploidy = sample_size ? variant_dosages.size() / sample_size : 1;
      savvy::stride_reduce(variant_dosages, ploidy);

      for (auto dose_it = variant_dosages.begin(); dose_it != variant_dosages.end(); ++dose_it) {
	int i = dose_it.offset();	
	//std::cout << "i " << i << std::endl;
        if(genetest_sample_idx_vcfDosage[i] >= 0) {
	  //std::cout << "genetest_sample_idx_vcfDosage[i] " << genetest_sample_idx_vcfDosage[i] << std::endl;		
          if (std::isnan(*dose_it)) {
		//std::cout << "here" << std::endl;
      //      dosagesforOneMarker[genetest_sample_idx_vcfDosage[i]] = float(1);
            ++missing_cnt;
            indexforMissing.push_back(genetest_sample_idx_vcfDosage[i]);
	    indexforMissingAll.push_back(genetest_sample_idx_vcfDosage[i]);
          }else {   	
	    dosagesforOneMarker.push_back(*dose_it);
	    jIndexforOneMarker.push_back(cnt+1);
	    iIndexforOneMarker.push_back(genetest_sample_idx_vcfDosage[i]+1);
            //dosagesforOneMarker[genetest_sample_idx_vcfDosage[i]] = *dose_it;
            AC = AC + *dose_it;
          }
        }
	//std::cout << "missing_cnt: " << missing_cnt << std::endl; 
      //group_matrix[cnt * sample_size + dose_it.offset()] = *dose_it;
      }
      float AF;
      if(genetest_samplesize_vcfDosage == missing_cnt){
        AF = 0;
      }else{	
        AF = (float)(AC) / 2 / (float)(genetest_samplesize_vcfDosage - missing_cnt) ;
      }

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
	//std::cout << "MAF: " << MAF <<  std::endl;
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
	positions.push_back(it->position());
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
    result["MACs"] = MACs; 	
    result["positions"] = positions;
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
}
