# SAIGE

SAIGE is an R package that implements the Scalable and Accurate Implementation of Generalized mixed model that uses the 
saddlepoint approximation (SPA)(mhof, J. P. , 1961; Kuonen, D. 1999; Dey, R. et.al 2017) 
and large scale optimization techniques to calibrate case-control ratios in logistic mixed model score tests
(Chen, H. et al. 2016) in large PheWAS. 

SAIGE can take dosage files in bgen, plain text (gzipped file is supported), or VCF format.

*This R package is still under development

The SAIGE manuscript can be found in the bioRxiv https://www.biorxiv.org/content/early/2017/11/01/212357.article-metrics

## Installation

Installation from the binary file in linux

    R CMD INSTALL SAIGE_0.16_R_x86_64-pc-linux-gnu.tar.gz

The following R pakages need to be installed for running SAIGE:

*Rcpp, RcppArmadillo, RcppParallel, data.table, SPAtest, RcppEigen, Matrix, methods*

The source code is SAIGE_0.16.tar.gz

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

## UK Biobank GWAS Results

The GWAS results for binary phenotypes in UK Biobank using SAIGE are currently available for public download at
https://www.dropbox.com/sh/wuj4y8wsqjz78om/AAACfAJK54KtvnzSTAoaZTLma?dl=0

We will continue to populate the public download repository with results for all UK Biobank phenotypes (> 1,600) with the PheCodes3 constructed based on ICD codes. 
*This research has been conducted using the UK Biobank Resource under application number 24460. 





