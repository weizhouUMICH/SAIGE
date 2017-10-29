# SAIGE

SAIGE is an R package that implements the Scalable and Accurate Implementation of Generalized mixed model that uses the 
saddlepoint approximation (SPA)(mhof, J. P. , 1961; Kuonen, D. 1999; Dey, R. et.al 2017) 
and large scale optimization techniques to calibrate case-control ratios in logistic mixed model score tests
(Chen, H. et al. 2016) in large PheWAS. SAIGE can take dosage files in either bgen format and plain txt format

*This R package is still under development

## Installation

Installation from the binary file in linux

    R CMD INSTALL SAIGE_0.13_R_x86_64-pc-linux-gnu.tar.gz

The following R pakages need to be installed for running SAIGE:

Rcpp, RcppArmadillo, RcppParallel, data.table, SPAtest, RcppEigen, Matrix, methods

## Running SAIGE for genetic association analysis




