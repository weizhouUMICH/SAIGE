#ifndef SCORE_HPP
#define SCORE_HPP

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <sstream>
#include <time.h>
#include <Rcpp.h>
#include <stdint.h>

//namespace SCORE {
	
class ScoreClass {
    private:
      std::vector<bool> m_LOCOvec;
      bool m_LOCO;
      arma::mat m_XXVX_inv;
      arma::mat m_XV;
      arma::mat m_XVX;
      arma::mat m_XVX_inv_XV;
      arma::mat m_X;
      arma::vec m_res;
      arma::vec m_mu;
      //arma::vec m_y;
      arma::vec m_mu2;
      arma::vec m_tauvec;
      double m_varRatio;
      arma::vec  m_S_a;
      std::string m_traitType;

    public:
      arma::vec m_y;
      bool m_isOutputAFinCaseCtrl;
      bool m_isOutputHetHomCountsinCaseCtrl;
      arma::uvec m_case_indices;
      arma::uvec m_ctrl_indices;
      arma::uvec m_case_hom_indices;
      arma::uvec m_case_het_indices;
      arma::uvec m_ctrl_hom_indices;
      arma::uvec m_ctrl_het_indices;

      
      void assignforScoreTest(bool t_LOCO, std::vector<bool> & t_LOCOVec,  arma::mat & t_XVX, arma::mat & t_XXVX_inv,  arma::mat & t_XV, arma::mat & t_XVX_inv_XV, arma::mat & t_X, arma::vec & t_S_a,  arma::vec & t_res,  arma::vec & t_mu2, arma::vec & t_mu, double t_varRatio, arma::vec & t_tauvec, std::string t_traitType, bool t_isOutputAFinCaseCtrl, bool t_isOutputHetHomCountsinCaseCtrl, arma::vec & t_y);
      void scoreTest(arma::vec t_GVec,
                     double& t_Beta, 
                     double& t_seBeta, 
                     std::string& t_pval_str, 
                     double t_altFreq,
		     double &t_Tstat, 
		     double &t_var1, 
		     double &t_var2,
		     arma::vec & t_gtilde);
      void scoreTestFast(arma::vec t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
		     double &t_Tstat,
                     double &t_var1,
                     double &t_var2);

      void get_mu(arma::vec & t_mu);
      arma::vec getadjG(arma::vec t_GVec);

};

//}
#endif
