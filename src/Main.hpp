#ifndef MAIN_HPP
#define MAIN_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

void setMarker_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
                               unsigned int t_omp_num_threads,
                               bool t_isOutputMoreDetails,
			       int t_marker_chunksize,
			       double g_dosage_zerod_cutoff,
                               double g_dosage_zerod_MAC_cutoff);


void setRegion_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               arma::vec t_max_maf_region,
                               unsigned int t_max_markers_region,
                               unsigned int t_omp_num_threads,
                               std::string t_method_to_CollapseUltraRare,
                               double t_MACCutoff_to_CollapseUltraRare,
                               double t_DosageCutoff_for_UltraRarePresence,
			       double g_dosage_zerod_cutoff,
                               double g_dosage_zerod_MAC_cutoff);

Rcpp::DataFrame mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::string & t_traitType,
                           arma::ivec & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation);

void Unified_getOneMarker(std::string & t_genoType,   // "PLINK", "BGEN"
                               uint32_t & t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool & t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint>& t_indexForMissing,     // index of missing genotype data
                               bool & t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint>& t_indexForNonZero,
			       arma::vec & t_GVec);

void Unified_getMarkerPval(
                           arma::vec & t_GVec,
                           bool t_isOnlyOutputNonZero,
                           arma::uvec & t_indexForNonZero_vec,
                           arma::uvec & t_indexForZero_vec,
                           double& t_Beta,
                           double& t_seBeta,
                           double& t_pval,
                           double& t_pval_noSPA,
                           double& t_Tstat,
			   double& t_gy,
                           double& t_varT,
                           double t_altFreq,
                           bool & t_isSPAConverge,
                           arma::vec & t_gtilde,
                           bool & is_gtilde,
			   bool  is_region,
			   arma::vec & t_P2Vec,
			   bool t_isCondition, 
			   double& t_Beta_c,
                           double& t_seBeta_c,
                           double& t_pval_c,
                           double& t_pval_noSPA_c,
                           double& t_Tstat_c,
                           double& t_varT_c,
			   arma::rowvec & t_G1tilde_P_G2tilde_Vec);


void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> & t_SampleInModel,
                      std::string t_AlleleOrder);



void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> & t_SampleInBgen,
                     std::vector<std::string> & t_SampleInModel,
                     std::string t_AlleleOrder);

void setSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
        arma::mat & t_X,
        arma::vec &  t_S_a,
        arma::vec & t_res,
        arma::vec & t_mu2,
        arma::vec & t_mu,
        arma::vec & t_varRatio,
        arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::vec & t_tauvec,
        std::string t_traitType,
        arma::vec & t_y,
        std::string t_impute_method,
        bool t_flagSparseGRM,
        arma::umat & t_locationMat,
        arma::vec & t_valueVec,
        int t_dimNum,
        bool t_isCondition,
        std::vector<uint32_t> & t_condition_genoIndex);


void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint32_t> & t_genoIndex,
                           unsigned int t_n);

void assign_conditionMarkers_factors_binary_region(
                           arma::vec & scalefactor_G2_cond);

#endif
