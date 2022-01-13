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

double sum_arma1(arma::vec& X);

double add_logp(double p1, double p2);

arma::vec arma_sub_cond(arma::vec x, arma::uvec iu);
