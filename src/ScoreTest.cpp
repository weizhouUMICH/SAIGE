// [[Rcpp::depends(BH)]]

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
#include <boost/math/distributions/chi_squared.hpp>

#include "ScoreTest.hpp"

//namespace SCORE {

void ScoreClass::assignforScoreTest(bool t_LOCO, std::vector<bool> & t_LOCOVec, arma::mat & t_XVX, arma::mat & t_XXVX_inv,  arma::mat & t_XV, arma::mat & t_XVX_inv_XV, arma::mat & t_X, arma::vec &  t_S_a,  arma::vec & t_res,  arma::vec & t_mu2, arma::vec & t_mu, double t_varRatio, arma::vec & t_tauvec, std::string t_traitType, bool t_isOutputAFinCaseCtrl, bool t_isOutputHetHomCountsinCaseCtrl, arma::vec & t_y){
    m_LOCO = t_LOCO;
    m_LOCOvec = t_LOCOVec;
    m_XVX = t_XVX;
    m_XV = t_XV;
    m_XXVX_inv = t_XXVX_inv;
    m_XVX_inv_XV = t_XVX_inv_XV;
    m_X = t_X;
    m_S_a = t_S_a;
    m_res = t_res;
    m_mu2 = t_mu2;
    m_mu = t_mu;
    m_varRatio = t_varRatio;
    m_tauvec = t_tauvec;  
    m_traitType = t_traitType;
    m_isOutputAFinCaseCtrl = t_isOutputAFinCaseCtrl;
    m_isOutputHetHomCountsinCaseCtrl = t_isOutputHetHomCountsinCaseCtrl;
    m_y = t_y;

    //arma::uvec m_case_indices;
    //arma::uvec m_ctrl_indices;

    if(m_traitType == "binary"){
      if(m_isOutputAFinCaseCtrl){
        m_case_indices = arma::find(m_y == 1);
        m_ctrl_indices = arma::find(m_y == 0);
      }
    }	    
}	


void ScoreClass::scoreTest(arma::vec t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde){

    //arma::mat g, g1, g2;
    //arma::mat Sm, var2m;
    arma::vec g, g1, g2;
    arma::vec Sm, var2m;
    double S, var2;
    g = m_XV * t_GVec;
    g1 = t_GVec - m_XXVX_inv * g;
    t_gtilde = g1;
    S = dot(g1, m_res);
    S = S/m_tauvec[0];
    g2 = arma::square(g1);
    var2m = dot(g2 , m_mu2)*m_tauvec[0];
    var2 = var2m(0,0);
    double var1 = var2 * m_varRatio;
    double stat = S*S/var1;
    double t_pval;
    //std::cout << "min double " <<  std::numeric_limits<double>::min() << std::endl;
    if (var1 < std::numeric_limits<double>::min()){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }
    char pValueBuf[100];
    if (t_pval != 0)	    
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr; 

    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(stat);
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "S: " << S << std::endl;
    //std::cout << "var1: " << var1 << std::endl;
    //std::cout << "var2: " << var2 << std::endl;
}


void ScoreClass::scoreTestFast(arma::vec t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){
    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    arma::vec Z = A1.t() * g1;
    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    if(m_traitType == "binary"){
      mu21  = m_mu2.elem(t_indexForNonZero);
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);	
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2; 
    }	    
    var1 = var2 * m_varRatio;
    S1 = dot(res1, g1_tilde);
    arma::mat res1X1_temp = (res1.t()) * X1;
    arma::vec res1X1 = res1X1_temp.t();
    S_a2 = m_S_a - res1X1;
    S2 = - arma::dot(S_a2,  Z);
    S = S1 + S2;
    S = S/m_tauvec[0];	    
    double stat = S*S/var1;
    double t_pval;
    if (var1 < std::numeric_limits<double>::min()){
	t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1); 
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }
    char pValueBuf[100];
    if (t_pval != 0)	    
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr; 
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(stat);
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "S: " << S << std::endl;
    //std::cout << "var1: " << var1 << std::endl;
    //std::cout << "var2: " << var2 << std::endl;
    //std::cout << "t_Tstat: " << t_Tstat << std::endl;
}

void ScoreClass::get_mu(arma::vec & t_mu){
    t_mu = m_mu;
}

arma::vec ScoreClass::getadjG(arma::vec t_GVec){
    arma::vec g;
    g = m_XV * t_GVec;
    g = t_GVec - m_XXVX_inv * g;
    return g;
}
//}

//};
//};
