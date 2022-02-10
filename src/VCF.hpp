#ifndef VCF_HPP
#define VCF_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "savvy/reader.hpp"
#include "savvy/region.hpp"
#include "variant_group_iterator.hpp" 



namespace VCF {
 
  class VcfClass{
  private:

 // vcf file and the index file (.tbi)
   std::string m_vcfFile, m_vcfIndexFile;
   //field in vcf to test. GT or DS 
   std::string m_chr;
   int32_t m_start = 0;
   int32_t m_end = 0;
   std::string m_markerIDsToExclude;
   std::string m_markerIDsToInclude; 
   std::vector<int32_t> m_posSampleInVcf;
   std::vector<int32_t> m_SampleInModel, m_posSampleInModel;
   bool m_isDropMissingDosages;
   bool m_isSparseDosagesInVcf; 

   uint32_t m_M0, m_M;
   uint32_t m_N0, m_N;
   std::vector<std::string> m_MarkerInVcf;     // Variant identifier  
   std::vector<std::string> m_SampleInVcf;
 
   bool m_isVcfOpen; 
    
  public:
 
   //savvy objects for streaming vcf files
   //savvy::reader m_marker_file{""};
   //VcfHeader header;
   //savvy::variant<std::vector<float>> m_record;

   //savvy::variant_group_iterator m_it_;
   std::string m_fmtField;
   savvy::reader m_marker_file{""};
   variant_group_iterator m_it_;

   VcfClass(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            bool t_isSparseDosageInVcf,
            std::vector<std::string> t_SampleInModel);
 
 
   // setup VcfClass
   bool setVcfObj(const std::string & t_vcfFileName,
                  const std::string & t_vcfFileIndex,
                  const std::string & t_vcfField);
   //set up the iterator 
   void set_iterator(std::string & variantList);
   void set_iterator(std::string & chrom, int & beg_pd, int & end_pd);

   void move_forward_iterator(int i);

   bool check_iterator_end();
   void setPosSampleInVcf(std::vector<std::string> & t_SampleInModel);
   void setIsSparseDosageInVcf (bool t_isSparseDosageInVcf);
 
   // get t_dosage/genotypes of one marker

   void getOneMarker(
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
                                  bool & t_isBoolRead,
                                  arma::vec & dosages,
				  bool t_isImputation);


    void getSampleIDlist_vcfMatrix(); 
    uint32_t getN0(){return m_N0;}
    uint32_t getN(){return m_N;}
    uint32_t getM0(){return m_M0;}
    uint32_t getM(){return m_M;}
 };
 
}

#endif
