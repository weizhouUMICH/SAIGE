// Important: this definition ensures Armadillo enables SuperLU
#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat mult_sp_sp_to_sp(const arma::sp_mat& a, const arma::sp_mat& b) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a * b);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    	
    return result;
}

// [[Rcpp::export]]
arma::sp_mat mult_sp_den_to_sp(const arma::sp_mat& a, const arma::mat& b) {
    // sparse x dense -> sparse
    arma::sp_mat result(a * b);
    return result;
}

// [[Rcpp::export]]
arma::sp_mat mult_den_sp_to_sp(const arma::mat& a, const arma::sp_mat& b) {
    // dense x sparse -> sparse
    arma::sp_mat result(a * b);
    return result;
}


// [[Rcpp::export]]
arma::sp_mat gen_sp(const arma::sp_mat& a) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;

    return result;
}
