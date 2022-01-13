//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


void SPA(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge);

void SPA_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge);
