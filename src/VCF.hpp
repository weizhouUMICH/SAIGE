#ifndef VCF_HPP
#define VCF_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "../thirdParty/cget/include/savvy/reader.hpp"
#include "../thirdParty/cget/include/savvy/varint.hpp"
#include "../thirdParty/cget/include/savvy/sav_reader.hpp"
 
 
namespace VCF {
 
  class VcfClass{
  private:

 // vcf file and the index file (.tbi)
   std::string m_vcfFile, m_vcfIndexFile;
   std::string m_vcfField;
   std::string m_chr;
   int32_t m_start = 0;
   int32_t m_end = 0;
   std::string m_markerIDsToExclude;
   std::string m_markerIDsToInclude; 
   std::vector<uint32_t> m_posSampleInVcf;
   std::vector<uint32_t> m_SampleInModel, m_posSampleInModel;
   bool m_isDropMissingDosages;
 
   //savvy objects for streaming vcf files
   savvy::indexed_reader m_reader;
   //VcfHeader header;
   savvy::variant<std::vector<float>> m_record;
   //field in vcf to test. GT or DS 
   uint32_t m_M0, m_M;
   uint32_t m_N0, m_N;
   std::vector<std::string> m_MarkerInVcf;     // Variant identifier  
   std::vector<std::string> m_SampleInVcf;
 
   bool m_isVcfOpen; 
   bool m_isReadVariant;
 
   double m_minMAF,m_maxMAF;
   double m_minInfo;
   
  public:
 
 
   VcfClass(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::string t_chr,
            int32_t t_start,
            int32_t t_end,
            bool t_isDropMissingDosages,
            double t_minMAF,
            double t_minInfo,
            double t_maxMAF,
 	    std::vector<std::string> t_SampleInModel);
 
 
   // setup VcfClass
   void setVcfObj(const std::string t_vcfFileName,
                  const std::string t_vcfFileIndex,
                  const std::string t_vcfField,
                  const std::string t_chr,
                  int32_t t_start,
                  int32_t t_end,
                  bool t_isDropMissingDosages,
                  double t_minMAF,
                  double t_minInfo,
                  double t_maxMAF);
   
   
   void setPosSampleInVcf(std::vector<std::string> & t_SampleInModel);
 
   // get t_dosage/genotypes of one marker
   bool getOneMarker(uint32_t t_gIndex,
 		     double& t_altFreq, 
                     double& t_missingRate,
                     std::vector<uint32_t>& t_posMissingGeno,
                     std::string& t_ref,
                     std::string& t_alt,
                     std::string& t_marker,
                     uint32_t& t_pd,
                     uint8_t& t_chr,
 	             std::string& t_chrMarker,
                     arma::vec & t_dosage);
 
 
 
    bool getMultiMarker(std::string& t_markerline,
                   std::vector<uint32_t>& t_posMissingGeno,
                   std::vector<double> & t_altFreqMulti,
                   std::vector<double> & t_altCountsMulti,
                   std::vector<std::string> & t_refMulti,
                   std::vector<std::string> & t_altMulti,
                   std::vector<std::string> & t_markerMulti,
                   std::vector<uint32_t> & t_pdMulti,
                   std::vector<uint8_t> & t_chrMarkerMulti,
                   std::vector<double> & t_dosageMulti,
                   std::vector<uint32_t> & t_iIndexMulti,
                   std::vector<uint32_t> & t_jIndexMulti
                   );
 
 
 
 
    uint32_t getN0(){return m_N0;}
    uint32_t getN(){return m_N;}
    uint32_t getM0(){return m_M0;}
    uint32_t getM(){return m_M;}
 };
 
}
 
#endif
