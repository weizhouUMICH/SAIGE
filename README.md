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

*Rcpp, RcppArmadillo, RcppParallel, data.table, SPAtest, RcppEigen, Matrix, methods*

## Running SAIGE

SAIGE contains 2 main steps:

1. Fitting the null logistic mixed model to estiamte variance component and other model parameters

    Run the **fitNULLGLMM** function for step 1
    
2. Testing for association between each genetic variant and phenotypes by applyting SPA to the score test
    
    Run the **SPAGMMATtest** function for step 2
    
### Examples

Examplary data and script can be found in ./extdata. Run 

    bash cmd.sh

to run the 2 steps. 

The R package optparse is required to run this script




