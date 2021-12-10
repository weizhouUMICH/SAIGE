// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
 
#include "../thirdParty/cget/include/savvy/reader.hpp"
#include "../thirdParty/cget/include/savvy/varint.hpp"
#include "../thirdParty/cget/include/savvy/sav_reader.hpp"
#include "../thirdParty/cget/include/savvy/variant_group_iterator.hpp"
 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <stdlib.h>
#include <cstring>
 
#include "VCF.hpp"


 
namespace VCF {
 
  VcfClass::VcfClass(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::string t_chr,
            int32_t t_start,
            int32_t t_end,
            bool t_isDropMissingDosageInVcf,
	    bool t_isSparseDosageInVcf,
            std::vector<std::string> t_SampleInModel)
   {
     setVcfObj(t_vcfFileName,
              t_vcfFileIndex,
              t_vcfField,
              t_chr,
              t_start,
              t_end);
 
     setPosSampleInVcf(t_SampleInModel);
     setIsDropMissingDosagesInVcf(t_isDropmissingdosagesInVcf);
     setIsSparseDosageInVcf(t_isSparseDosageInVcf); 
   }
 
 
 
 
 
   void VcfClass::setVcfObj(const std::string t_vcfFileName,
                  const std::string t_vcfFileIndex,
                  const std::string t_vcfField,
                  const std::string t_chr,
                  int32_t t_start,
                  int32_t t_end)
   {
     if(t_vcfField == "DS"){
       m_reader = savvy::indexed_reader(t_vcfFileName, {t_chr, std::uint32_t(t_start), std::uint32_t(t_end)}, savvy::fmt::ds);
 
     }else if(t_vcfField == "GT"){
       m_reader = savvy::indexed_reader(t_vcfFileName, {t_chr, std::uint32_t(t_start), std::uint32_t(t_end)}, savvy::fmt::ac);	    
     }
 
     bool isVcfOpen = m_reader.good();
 
     if(isVcfOpen){
       std::cout << "Open VCF done" << std::endl;
       m_vcfField = t_vcfField;
       std::cout << "To read the field " << m_vcfField << std::endl;
       std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << m_reader.headers().size() << std::endl;
       std::cout << "Number of samples in in the vcf file: " << m_reader.samples().size() << std::endl;
       m_SampleInVcf = m_reader.samples();
       m_N0 = m_SampleInVcf.size();
     }else{
       std::cout << "WARNING: Open VCF failed" << std::endl;
     }
 
   }
 
   void VcfClass::setPosSampleInVcf(std::vector<std::string> & t_SampleInModel)
   {
     std::cout << "Setting position of samples in VCF files...." << std::endl;
     m_N = t_SampleInModel.size();

  // updated by BWJ on 03/14/2021

     Rcpp::CharacterVector SampleInVcf(m_N0);
     for(uint32_t i = 0; i < m_N0; i++)
       SampleInVcf(i) = m_SampleInVcf.at(i);

     Rcpp::CharacterVector SampleInModel(m_N);
     for(uint32_t i = 0; i < m_N; i++)
       SampleInModel(i) = t_SampleInModel.at(i);

     Rcpp::IntegerVector posSampleInVcf = Rcpp::match(SampleInModel, SampleInVcf);
     for(uint32_t i = 0; i < m_N; i++){
       if(Rcpp::IntegerVector::is_na(posSampleInVcf.at(i)))
          Rcpp::stop("At least one subject requested is not in VCF file.");
     }

     Rcpp::IntegerVector posSampleInModel = Rcpp::match(SampleInVcf, SampleInModel);
     m_posSampleInModel.resize(m_N0);
     for(uint32_t i = 0; i < m_N0; i++){
       if(Rcpp::IntegerVector::is_na(posSampleInModel.at(i))){
         m_posSampleInModel.at(i) = -1;
       }else{
         m_posSampleInModel.at(i) = posSampleInModel.at(i) - 1;   // convert "starting from 1" to "starting from 0"
       }
     }
  }


 
 
   // get dosages/genotypes of one marker   
   bool VcfClass::getOneMarker(uint64_t t_gIndex,        // different meanings for different genoType
                                  std::string& t_ref,       // REF allele
                                  std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                  std::string& t_marker,    // marker ID extracted from genotype file
                                  uint32_t& t_pd,           // base position
                                  std::string& t_chr,       // chromosome
                                  double& t_altFreq,        // frequency of ALT allele
                                  double& t_altCounts,      // counts of ALT allele
                                  double& t_missingRate,    // missing rate
                                  double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                                  bool t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                                  std::vector<uint32_t>& t_indexForMissingforOneMarker,     // index of missing genotype data
                                  bool t_isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                  std::vector<uint32_t>& t_indexForNonZero,
				  //if true, the marker has been read successfully
				  bool t_isReadVariant,
				  std::vector<double> & t_dosage) 
   {	  
     m_isReadVariant = m_reader >> m_record; 
     if(m_isReadVariant){
       t_marker = m_record.prop("ID");
       t_chr = m_record.chromosome();
       t_pd = m_record.position();
       int numAlt = 1;
       t_ref = m_record.ref();
       t_alt = m_record.alt();
       double dosage;
       double AC = 0;
       t_dosage.clear();
       t_dosage.reserve(m_N);

       if(!t_isOnlyOutputNonZero){
         t_dosage.resize(m_N);
       }
       //arma::vec & t_dosage;
       //t_dosage.set_size(m_N);
       std::size_t missing_cnt = 0;
       int i;
       for (std::vector<float>::iterator it = m_record.data().begin(); it != m_record.data().end(); ++it) {
         i = std::distance(m_record.data().begin(), it);
         if(m_posSampleInModel.at(i) != -1) {
           if (std::isnan(*it)) {
             t_dosage(m_posSampleInModel.at(i)) = -1;
             ++missing_cnt;
             t_indexForMissingforOneMarker.push_back(m_posSampleInModel.at(i));
           }else {
	     dosage = *it;	   
             if(!t_isOnlyOutputNonZero){		     
	       t_dosage(m_posSampleInModel.at(i)) = dosage;
	     }else{
               if(dosage > 0){
                 t_dosage.push_back(dosage);
		 t_indexForNonZero.push_back(m_posSampleInModel.at(i) + 1);
	       }	       
             }		     
             AC = AC + dosage;
           }
         }
       }
 
       if(missing_cnt > 0){
         if(missing_cnt == m_N){
           t_altFreq = 0;
         }else{
           t_altFreq = AC / 2 / (double)(m_N - missing_cnt);
         }
         if(!m_isDropMissingDosages){	
  	   double imputeDosage = 2*t_altFreq;
           for (unsigned int i = 0; i < t_indexForMissingforOneMarker.size(); i++){
             t_dosage(t_indexForMissingforOneMarker.at(i)) = imputeDosage;
 	     AC = AC + imputeDosage;
           }         	
 	 }  
       }
 
     }else{
       std::cout << "Reach the end of the vcf file" << std::endl;
     }
     t_isReadVariant = m_isReadVariant;
     return(t_isReadVariant);    	  
   }
 
   bool VcfClass::getMultiMarker(std::string& t_markerline,		
                   std::vector<std::string> & t_refMulti,
                   std::vector<std::string> & t_altMulti,
                   std::vector<std::string> & t_markerMulti,
                   std::vector<uint32_t> & t_pdMulti,
                   std::vector<uint8_t> & t_chrMulti,
                   std::vector<double> & t_altFreqMulti,
                   std::vector<double> & t_altCountsMulti,
                   std::vector<double> & t_macMulti,
		   std::vector<double> & t_missingRateMulti,
		   std::vector<double> & t_imputeInfoMulti,
                   bool t_isOutputIndexForMissing,
                   std::vector<uint32_t>& t_indexForMissingforMulti,
                   std::vector<uint32_t> & t_iIndexMulti,
                   std::vector<uint32_t> & t_jIndexMulti,
		   std::vector<double> & t_dosageMulti){
	   
     savvy::variant_group_iterator<savvy::compressed_vector<double>> it(m_reader, t_markerline);
     savvy::variant_group_iterator<savvy::compressed_vector<double>> end{};
     std::vector<double> t_dosageMulti;
     t_dosageMulti.clear();
     std::vector<double> t_dosageforOneMarker;
     std::vector<uint32_t> iIndexforOneMarker;
     std::vector<uint32_t> jIndexforOneMarker;
     std::vector<uint32_t> t_indexForMissingforOneMarker;
     double freq, maf, mac;

     if (it != end){
       m_isReadVariant = true;	    
       int missing_cnt = 0;
       int cnt = 0;
       double AC = 0;

       for ( ; it != end; ++it)
       {
         AC = 0;
         missing_cnt = 0;
         t_indexForMissingforOneMarker.clear();
         t_dosageforOneMarker.clear();
         iIndexforOneMarker.clear();
         jIndexforOneMarker.clear();
         it.group_id();
         it.sites();
         std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alt();
         std::string markerInfo_str = it->prop("R2");
         double markerInfo = strtod((markerInfo_str).c_str(),0);
         for (auto dose_it = it->data().begin(); dose_it != it->data().end(); ++dose_it){
           int lengthi = std::distance(it->data().begin(), it->data().end());
           int i = dose_it.offset();
  
           if(m_posSampleInModel.at(i) != -1) {
             if (std::isnan(*dose_it)) {
               ++missing_cnt;
 	       t_indexForMissingforOneMarker.push_back(m_posSampleInModel.at(i));
               t_posMissingGeno.push_back(m_posSampleInModel.at(i));
             }else {
               if(*dose_it > 0){   		    
                 t_dosageforOneMarker.push_back(*dose_it);
 		 jIndexforOneMarker.push_back(cnt+1);
 		 iIndexforOneMarker.push_back(m_posSampleInModel.at(i)+1); //1-based for R
                 AC = AC + *dose_it;
               }
 	     }
 	   }    
         }
 
         if(missing_cnt > 0){
           if(missing_cnt == m_N){
             freq = 0;
           }else{
             freq = AC / 2 / (double)(m_N - missing_cnt);
           }
         }
         maf = freq;
	 mac = AC;
         if(freq > 0.5){
	   maf = 1-maf;
	   mac = 2*(double)(m_N - missing_cnt) - AC;
	 }	
 
  	 if(maf >= g_marker_minMAF_cutoff && maf <= g_marker_maxMAF_cutoff && markerInfo >= g_marker_minINFO_cutoff){
           if(missing_cnt > 0 && !m_isDropMissingDosages){
             //std::cout << "missing_cnt > 0!" << std::endl;
             double imputeDosage = 2*freq;
             for (unsigned int i = 0; i < missing_cnt; i++){
               t_dosageforOneMarker.push_back(imputeDosage);
               jIndexforOneMarker.push_back(cnt+1);
               iIndexforOneMarker.push_back(t_indexForMissingforOneMarker.at(i)+1); //1-based for R
 	      AC = AC + imputeDosage;
             }
	     mac = AC;
	     if(freq > 0.5){
               mac = 2*(double)(m_N) - AC;
             }
 	   }

           t_dosageMulti.insert(std::end(t_dosageMulti), std::begin(t_dosageforOneMarker), std::end(t_dosageforOneMarker));
           t_iIndexMulti.insert(std::end(t_iIndexMulti), std::begin(iIndexforOneMarker), std::end(iIndexforOneMarker));
 	   t_jIndexMulti.insert(std::end(t_jIndexMulti), std::begin(jIndexforOneMarker), std::end(jIndexforOneMarker));
           cnt = cnt + 1;
           t_markerMulti.push_back(marker_id);
           t_altFreqMulti.push_back(freq);
 	   t_altCountsMulti.push_back(AC);
	   t_macMulti.push_back(mac);
           t_pdMulti.push_back(it->position());
        }      
      }
 
     }else{
       m_isReadVariant = false;     
     }
     return(m_dosageMulti);
   }

void VcfClass::setIsSparseDosageInVcf (bool t_isSparseDosageInVcf){
  m_isSparseDosagesInVcf = t_isSparseDosageInVcf;
}

void VcfClass::setIsDropMissingDosagesInVcf (bool t_isDropmissingdosagesInVcf){
  m_isDropMissingDosagesInVcf = t_isDropmissingdosagesInVcf;
}

}
