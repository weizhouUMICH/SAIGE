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
#include "../thirdParty/bgen/3rd_party/boost_1_55_0/boost/math/distributions/chi_squared.hpp"


class ScoreTest {
    private:
      std::vector<bool> LOCOvec;
      bool LOCO;
      arma::fmat XXVX_inv_noLOCO;
      arma::fmat XV_inv_noLOCO;
      arma::fvec res_noLOCO;
      arma::fvec mu2_noLOCO;
      double varRatio;
    public:
      void assignforScoreTest(bool LOCO_ext, std::vector<bool> & LOCOVec_ext, arma::fmat & XXVX_inv_noLOCO_ext,  arma::fmat & XV_inv_noLOCO_ext,  arma::fvec & res_noLOCO_ext,  arma::fvec & mu2_noLOCO_ext, double varRatio_ext);
      //Rcpp::List Scoretest(const arma::mat& G, int chri);
	Rcpp::List Scoretest(const arma::fvec& G, int chri);
//      void Scoretest(const arma::fvec& G, int chri);
};
/*
    void ScoreTest::assignforScoreTest(bool LOCO_ext, std::vector<bool> & LOCOVec_ext, arma::fmat & XXVX_inv_noLOCO_ext,  arma::fmat & XV_inv_noLOCO_ext,  arma::fvec & res_noLOCO_ext,  arma::fvec & mu2_noLOCO_ext, double varRatio_ext) {
    LOCO = LOCO_ext;
    LOCOvec = LOCOVec_ext;
    XV_inv_noLOCO = XV_inv_noLOCO_ext;
    XXVX_inv_noLOCO = XXVX_inv_noLOCO_ext;
    res_noLOCO = res_noLOCO_ext;
    mu2_noLOCO = mu2_noLOCO_ext;
    varRatio = varRatio_ext;

}	
/*

  Rcpp::List Scoretest(const arma::sp_mat& G, int chri) {
    arma::colvec g, g1, g2;
    double S, var2;
    //if(!LOCO || !LOCOvec[chri-1]){	  
      g = XV_inv_noLOCO * G;
      g1 = g - XXVX_inv_noLOCO * g;
      S = arma::dot(g1, res_noLOCO);
      g2 = arma::square(g1);
      var2 = arma::dot(g2, mu2_noLOCO);
    //}else{
    //  g = XV_inv_list.slice(chri-1) * G;
    //  g1 = g - XXVX_inv_list.slice(chri-1) * g;
    //  S = arma::dot(g1, res_list.col(chri-1));
    //  g2 = arma::square(g1);
    //  var2 = arma::dot(g2, mu2_list.col(chri-1));
    //}

    double var1 = var2 * varRatio;
    double stat = S*S/var1;
    boost::math::chi_squared chisq_dist(1); 
    double pValue = boost::math::cdf(complement(chisq_dist, stat));

    char pValueBuf[100];
    //if (pValue != 0)
    //    sprintf(pValueBuf, "%.1E", pValue);
    //else {
    //    double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
    //    int exponent = floor(log10p);
    //    double fraction = pow(10.0, log10p - exponent);
    //    if (fraction >= 9.95) {
    //      fraction = 1;
   //       exponent++;
   //     }
   //     sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    //}

    double beta = S/var1;
    double se = fabs(beta) / sqrt(stat);

    return Rcpp::List::create(Rcpp::Named("BETA") = beta,
                              Rcpp::Named("SE") = se,
                              Rcpp::Named("Tstat")  = S,
                              Rcpp::Named("pval.noadj")  = pValue,
                              Rcpp::Named("is.converge")=TRUE,
                              Rcpp::Named("var1")=var1,
                              Rcpp::Named("var2")=var2);
    }
*/

//};
