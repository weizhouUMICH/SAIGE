
#ifndef UTIL_HPP
#define UTIL_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <sys/time.h>

const static std::unordered_map<std::string,int> string_to_case{
   {"best_guess",1},
   {"mean",2},
   {"minor",3}
};

double getWeights(std::string t_kernel, 
                  double t_freq, 
                  arma::vec t_wBeta);

void imputeGeno(arma::vec& GVec, 
                double freq, 
                std::vector<uint32_t> posMissingGeno);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);


bool imputeGenoAndFlip(arma::vec& t_GVec,
                       double & t_altFreq,
		       double & t_altCount,
                       std::vector<uint32_t> &  t_indexForMissing,
                       std::string t_impute_method,
                       double t_dosage_zerod_cutoff,
                       double t_dosage_zerod_MAC_cutoff,
                       double & t_MAC,
		       std::vector<uint> & t_indexZero,
                       std::vector<uint> & t_indexNonZero);

arma::vec getTime();

void printTime(arma::vec t1, arma::vec t2, std::string message);

double getinvStd(double t_freq);

// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
// void set_seed(unsigned int seed) {
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(seed);  
// };

// 
// arma::vec nb(int n){
//   return(Rcpp::rbinom(n,1,0.5));
// }

arma::vec nb(unsigned int n);

double sum_arma1(arma::vec& X);

double add_logp(double p1, double p2);

#endif

