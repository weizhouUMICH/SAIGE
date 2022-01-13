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
#include <boost/math/distributions/normal.hpp>
#include "SPA_binary.hpp"
#include "SPA_survival.hpp"
#include "utils.hpp"
#include "getMem.hpp"

// [[Rcpp::export]]
void SPA(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge){
        using namespace Rcpp;
        List result;
        double p1, p2;
        bool Isconverge = true;
	Rcpp::List outuni1;
	Rcpp::List outuni2;
        if( traitType == "binary"){
          outuni1 = getroot_K1_Binom(0, mu, g, q, tol);
          outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol);
        }else if(traitType == "timeToEvent"){
          outuni1 = getroot_K1_Poi(0, mu, g, q, tol);
          outuni2 = getroot_K1_Poi(0, mu, g, qinv, tol);
        }

        //double outuni1root = outuni1["root"];
        //double outuni2root = outuni2["root"];
        //bool Isconverge1 = outuni1["Isconverge"];
        //bool Isconverge2 = outuni2["Isconverge"];

        //std::cout << "outuni1root" << outuni1root << std::endl;
        //std::cout << "outuni2root" << outuni2root << std::endl;
        //std::cout << "Isconverge1" << Isconverge1 << std::endl;
        //std::cout << "Isconverge2" << Isconverge2 << std::endl;


        Rcpp::List getSaddle;
        Rcpp::List getSaddle2;
        if(outuni1["Isconverge"]  && outuni2["Isconverge"])
        {
                //std::cout << "q is " << q << " 3 qinv is " << qinv << std::endl;
                if( traitType == "binary"){
                  getSaddle = Get_Saddle_Prob_Binom(outuni1["root"], mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Binom(outuni2["root"], mu, g, qinv, logp);
                }else if(traitType == "timeToEvent"){
                  getSaddle = Get_Saddle_Prob_Poi(outuni1["root"], mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Poi(outuni2["root"], mu, g, qinv, logp);
                }

                if(getSaddle["isSaddle"]){
                        p1 = getSaddle["pval"];
                }else{

                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
                if(getSaddle2["isSaddle"]){
                        p2 = getSaddle2["pval"];
                }else{
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }
                //std::cout << "p1_nofast " << p1 << "p2 " << p2 << std::endl;

                if(logp)
                {
                        pval = add_logp(p1,p2);
                } else {
                        pval = std::abs(p1)+std::abs(p2);
                }
                Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        Isconverge=false;
                }
        isSPAConverge = Isconverge;
	//result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}


// [[Rcpp::export]]
void SPA_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge){

        using namespace Rcpp;
        List result;
        double p1, p2;
        bool Isconverge = true;
	Rcpp::List outuni1;
        Rcpp::List outuni2;
	double mem1, mem2;
	process_mem_usage(mem1, mem2);
   std::cout << "VM 5 a: " << mem1/1000000 << "; RSS 5 a: " << mem2/1000000 << std::endl;
        if( traitType == "binary"){
          outuni1 = getroot_K1_fast_Binom(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
          //double qinv = -1*q;
          outuni2 = getroot_K1_fast_Binom(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
          //std::cout << "outuni1root" << outuni1["root"] << std::endl;
          //std::cout << "outuni2root" << outuni2["root"] << std::endl;
        }else if(traitType == "timeToEvent"){
          outuni1 = getroot_K1_fast_Poi(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
          outuni2 = getroot_K1_fast_Poi(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
        }
	process_mem_usage(mem1, mem2);
   std::cout << "VM 5 b: " << mem1/1000000 << "; RSS 5 b: " << mem2/1000000 << std::endl;

        Rcpp::List getSaddle;
        Rcpp::List getSaddle2;
        if(outuni1["Isconverge"]  && outuni2["Isconverge"])
        {
          if( traitType == "binary"){
                getSaddle  = Get_Saddle_Prob_fast_Binom(outuni1["root"], mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
                getSaddle2 = Get_Saddle_Prob_fast_Binom(outuni2["root"], mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
	process_mem_usage(mem1, mem2);
   std::cout << "VM 5 c: " << mem1/1000000 << "; RSS 5 c: " << mem2/1000000 << std::endl;
          }else if(traitType == "timeToEvent"){
                getSaddle  = Get_Saddle_Prob_fast_Poi(outuni1["root"], mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
                getSaddle2 = Get_Saddle_Prob_fast_Poi(outuni2["root"], mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
          }
                if(getSaddle["isSaddle"]){
                        p1 = getSaddle["pval"];
                }else{

                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
	process_mem_usage(mem1, mem2);
   std::cout << "VM 5 d: " << mem1/1000000 << "; RSS 5 d: " << mem2/1000000 << std::endl;

                if(getSaddle2["isSaddle"]){
                        p2 = getSaddle2["pval"];

                }else{
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }
	process_mem_usage(mem1, mem2);
   std::cout << "VM 5 e: " << mem1/1000000 << "; RSS 5 e: " << mem2/1000000 << std::endl;

                if(logp){
                        pval = add_logp(p1,p2);
                } else {
                        pval = std::abs(p1)+std::abs(p2);
                        //std::cout << "p1 " << p1 << "p2 " << p2 << std::endl;
                }
                Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        Isconverge=false;
                }
	isSPAConverge = Isconverge;
        //result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}
