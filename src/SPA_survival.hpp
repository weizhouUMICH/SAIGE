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

	
double Korg_Poi(double t1, arma::vec & mu, arma::vec & g);
double K1_adj_Poi(double t1, arma::vec & mu, arma::vec & g, double q);
double K2_Poi(double t1, arma::vec & mu, arma::vec & g);
Rcpp::List getroot_K1_Poi(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter = 1000);
Rcpp::List Get_Saddle_Prob_Poi(double zeta,  arma::vec & mu, arma::vec & g, double q, bool logp =false);
Rcpp::List SPA_survival(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp=false);
double Korg_fast_Poi(double t1, arma::vec & mu, arma::vec & g,  arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);
double K1_adj_fast_Poi(double t1, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma);
double K2_fast_Poi(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma);
Rcpp::List getroot_K1_fast_Poi(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, int maxiter = 1000);
Rcpp::List Get_Saddle_Prob_fast_Poi(double zeta,  arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, bool logp=false);
Rcpp::List SPA_survival_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol);
