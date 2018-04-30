Table of Contents
=================

   * [Introduction](#introduction)
   * [Citation](#citation)
   * [Installation](#installation)
   * [Log for fixing bugs](#log-for-fixing-bugs)
   * [Notes for users before running jobs](#notes-for-users-before-running-jobs)
   * [Running SAIGE](#running-saige)
      * [Examples](#examples)
      * [Input files](#input-files)
         * [Step 1: Genotype file (for contructing the genetic relationship matrix)](#step-1-genotype-file-for-contructing-the-genetic-relationship-matrix)
         * [Step 1: Phenotype file](#step-1-phenotype-file)
         * [Step 2: Dosage file containing dosages for genetic variates to be tested](#step-2-dosage-file-containing-dosages-for-genetic-variates-to-be-tested)
         * [Step 2: Two files output by step 1](#step-2-two-files-output-by-step-1)
   * [UK Biobank GWAS Results](#uk-biobank-gwas-results)

# Introduction

SAIGE is an R package that implements the Scalable and Accurate Implementation of Generalized mixed model that uses the saddlepoint approximation (SPA)(mhof, J. P. , 1961; Kuonen, D. 1999; Dey, R. et.al 2017) 
and large scale optimization techniques to calibrate case-control ratios in logistic mixed model score tests
(Chen, H. et al. 2016) in large-scale GWAS 

SAIGE can take dosage files in bgen, plain text (gzipped file is supported), or VCF format.

*This R package is still under development

# Citation
The SAIGE manuscript can be found in the bioRxiv https://www.biorxiv.org/content/early/2017/11/01/212357.article-metrics

# Installation

Installation from the binary file in linux

    R CMD INSTALL SAIGE_0.XX_R_x86_64-pc-linux-gnu.tar.gz

The following R pakages need to be installed for running SAIGE:

*Rcpp, RcppArmadillo, RcppParallel, data.table, SPAtest, RcppEigen, Matrix, methods*

Source code for all versions are here https://www.dropbox.com/sh/zmlu1llpxd66pjl/AADFqdssvOBjbWZch6Q9zYNaa?dl=0

# Log for fixing bugs
* 0.29:
```
1. The colSums() error when there is no covariate has been fixed. 
2. BETA and Tstat are now for the alt allele for both quantitative and binary traits. Note that in version <= 0.26, for binary traits, BETA is for alt allele and for quantitative traits, BETA is for minor allele
3. Options for leave-one-chromosome-out (LOCO), cutoffs for the coefficient of variation (CV) for trace estimates and variance ratio estimates have been added, but these three options have not been extensively tested. CV is mainly for automatically determining whether the number of random markers selected is sufficient or not. If not, the number will be increased until the CV is lower than the specified cutoff.  
```
* 0.26: fixed a bug for the Tstat in the output
* 0.25: allow models with no covariates and GRM contruction using a large number of genetic markers (> 600,000)
* 0.24: centerVariable is no longer needed. QR transformation of the covariate matrix is automatically performed. Supports the dosage files in the VCF,BCF and SAV formats using the SAVVY library 

# Notes for users before running jobs
1. Since the SPA test always provides close to 0 p-values for variants with MAC < 3, please use at least minMAC = 3 to filter out the results
2. When the query is used for bgen files, please make sure there is no duplicate SNP ids in the list
3. If the error message "Error in setgeno(genofile, subSampleInGeno, memoryChunk) :
  vector::_M_range_check", try use a smaller memeoryChunk, such as 2
4. In version <= 0.26, for binary traits, BETA is for alt allele and for quantitative traits, BETA is for minor allele 
5. Please note that LOCO only works for autosomal genetic variants. For non-autosomal genetic variants, please leave LOCO=FALSE in step 2.

# Running SAIGE

SAIGE contains 2 main steps:

1. Fitting the null logistic mixed model to estiamte variance component and other model parameters

    Run the **fitNULLGLMM** function for step 1
    
2. Testing for association between each genetic variant and phenotypes by applyting SPA to the score test
    
    Run the **SPAGMMATtest** function for step 2
    
## Examples

Examplary data and script can be found in ./extdata. Run 

    bash cmd.sh

to run the 2 steps. 

The R package optparse is required to run this script

## Input files

### Step 1: Genotype file (for contructing the genetic relationship matrix)

You can use *plinkFile* to specify the genotype file. 

SAIGE takes the PLINK binary file for the genotypes and assumes the file prefix is the same one for .bed, .bim. and .fam

### Step 1: Phenotype file
You can use *phenoFile* to sepcify the phenotype file, use *phenoCol* to specify the column name for the phentoype (e.g. y), use *sampleIDColinphenoFile* to specify the column name for the sample ids (e.g. IID), use *covarColList* to specify the column names for covariates (e.g. c("x1","x2")) 

*Note: Current version of SAIGE does not support categorical covariates that have more than two categories*

### Step 2: Dosage file containing dosages for genetic variates to be tested
SAIGE takes dosage files in plain text (gzipped file is supported), BGEN, VCF, BCF and [SAV](https://github.com/statgen/savvy) format.

### Step 2: Two files output by step 1
*GMMATmodelFile* and *varianceRatioFile*

# UK Biobank GWAS Results
The GWAS results for binary phenotypes in UK Biobank using SAIGE are currently available for public download at
https://www.dropbox.com/sh/wuj4y8wsqjz78om/AAACfAJK54KtvnzSTAoaZTLma?dl=0

We will continue to populate the public download repository with results for all UK Biobank phenotypes (> 1,600) with the PheCodes3 constructed based on ICD codes. 
*This research has been conducted using the UK Biobank Resource under application number 24460. 





