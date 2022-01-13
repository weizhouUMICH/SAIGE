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
#include "utils.hpp"


// [[Rcpp::export]]
double sum_arma1(arma::vec& X) {
    double sum = 0;
    for (int i = 0; i < X.n_elem; ++i) {
        if (arma::is_finite(X(i)))
            sum += X(i);
    }
    return sum;
}


// [[Rcpp::export]]
double add_logp(double p1, double p2)
{
        using namespace Rcpp;
        p1 = std::abs(p1);
        p2 = std::abs(p2);
        double maxp = std::max(p1,p2);
        double  minp = std::min(p1,p2);
        double result = maxp+std::log(1+std::exp(minp-maxp));
        return(result);
}

// [[Rcpp::export]]
arma::vec arma_sub_cond(arma::vec x, arma::uvec iu) {

	arma::vec ids = x.elem(iu); // Find indices


    return ids;
}
