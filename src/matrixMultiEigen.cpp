// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>


// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> & A, Eigen::Map<Eigen::MatrixXd> & B){
     Eigen::MatrixXd C = A * B;
     return Rcpp::wrap(C);
}

