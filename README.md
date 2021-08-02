
Table of Contents
=================

   * [Introduction](#introduction)
   * [Citation](#citation)
   * [How to install SAIGE and SAIGE-GENE](#how-to-install-and-run-saige-and-saige-gene)
   * [Notes for users before running jobs](#notes-for-users-before-running-jobs)
   * [UK Biobank GWAS Results](#uk-biobank-gwas-results)
   * [Log for fixing bugs](#log-for-fixing-bugs)
   

# Introduction

##Current version is 0.44.6.2 (Updated on August 2, 2021) - 0.44.6.2 add extdata/extractNglmm.R to extract the effective sample size without running Step 1. extdata/cmd_extractNeff.sh has the pipeline. The effective sample size (Nglmm) is differently calculated than the previous versions.

## Previous version is 0.44.6.1 (Updated on July 16, 2021) - SAIGE-GENE+: for group tests, collpasing ultra-rare variants with MAC <= 10. Set --method_to_CollapseUltraRare="absence_or_presence" as default to collpase ultra-rare varaints with MAC <= 10. SAIGE-GENE+ has well controlled type I error rates when the maximum MAF cutoff (maxMAFforGroupTest) is lower than 1%, e.g. 0.01% or 0.1%. Tests with multiple MAF cutoffs and variant annotations can be combined using the Cauchy combination (function CCT)

## Please re-install 0.44.2 if you installed this verion on March 31. 

## For BGEN input, 8 bits are required. 

## For BGEN input in step 2 with missing dosages, Please use version 0.38 or later.


SAIGE is an R package with Scalable and Accurate Implementation of Generalized mixed model (Chen, H. et al. 2016). It accounts for sample relatedness and is feasible for genetic association tests in large cohorts and biobanks (N > 400,000).

SAIGE performs single-variant association tests for binary traits and quantitative taits. For binary traits, SAIGE uses the saddlepoint approximation (SPA)(mhof, J. P. , 1961; Kuonen, D. 1999; Dey, R. et.al 2017) to account for case-control imbalance.

SAIGE-GENE (implemented in the SAIGE R package) performs gene- or region-based association tests (Burde, SKAT, SKAT-O) for binary traits and quantitative traits. Note: SAIGE-GENE accounts for case-control imbalance in gene-based tests (>= 0.35.8.5)


# Citation
The SAIGE manuscript:
Wei Zhou, Jonas B. Nielsen, Lars G. Fritsche, Maiken B. Elvestad, Brooke Wolford, Maoxuan Lin, Kristian Hveem, Hyun Min Kang, Goncalo R. Abecasis, Cristen J. Willer*, Seunggeun Lee* “Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies.” Nature Genetics 50, 1335–1341 (2018)

The SAIGE-GENE pre-print:
https://www.biorxiv.org/content/10.1101/583278v2


# How to install and run SAIGE and SAIGE-GENE


## Install SAIGE/SAIGE-GENE

### List of dependencies: 

* R-3.6.1, gcc >= 5.4.0, cmake 3.14.1, [cget](https://cget.readthedocs.io/en/latest/src/intro.html#installing-cget)
* R packages: "R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "SKAT","MetaSKAT"
* /extdata/install_packages.R can be used to install the R packages
* SAIGE v0.39.2 depends on the SPAtest v3.1.2
* MetaSKAT is currently not available on CRAN. Please install it from github using R
  ``` 
   devtools::install_github("leeshawn/MetaSKAT") 
  ```

###  Install SAIGE from conda

![r-saige](https://anaconda.org/bioconda/r-saige/badges/version.svg)
![latest_update](https://anaconda.org/bioconda/r-saige/badges/latest_release_date.svg)

To install saige from conda simply create environment with latest version of R and saige:
```
conda create -n saige -c conda-forge -c bioconda "r-base>=4.0" r-saige
conda activate saige
```

More info on [r-saige conda package](https://anaconda.org/bioconda/r-saige) and available versions can be found in the [issue #272](https://github.com/weizhouUMICH/SAIGE/issues/272).

###  Install SAIGE using the conda environment

1. Create a conda environment using 
     ([conda environment file](https://github.com/weizhouUMICH/SAIGE/blob/master/conda_env/environment-RSAIGE.yml)) 
     Here is a link to download the [conda environment file](https://raw.githubusercontent.com/weizhouUMICH/SAIGE/master/conda_env/environment-RSAIGE.yml)

     After downloading environment-RSAIGE.yml, run following command
     ```
       conda env create -f environment-RSAIGE.yml
   ```

2. Activate the conda environment RSAIGE

     ```
       conda activate RSAIGE
       FLAGPATH=`which python | sed 's|/bin/python$||'`
       export LDFLAGS="-L${FLAGPATH}/lib"
       export CPPFLAGS="-I${FLAGPATH}/include"
     ```
Please make sure to set up the LDFLAGS and CPPFLAGS using export (the last two command lines), so libraries can be linked correctly when the SAIGE source code is compiled. Note: [Here](https://github.com/weizhouUMICH/SAIGE/blob/master/conda_env/createCondaEnvSAIGE_steps.txt) are the steps to create the conda environment file 

3. Open R, run following script to install the MetaSKAT R library.
   
     ```
       devtools::install_github("leeshawn/MetaSKAT") 
     ```

4. Install SAIGE from the source code. 

     Method 1: 

     ```
       src_branch=master
       repo_src_url=https://github.com/weizhouUMICH/SAIGE
       git clone --depth 1 -b $src_branch $repo_src_url

       R CMD INSTALL --library=path_to_final_SAIGE_library SAIGE
     ```
     
     When call SAIGE in R, set lib.loc=path_to_final_SAIGE_library   

     ```
       library(SAIGE, lib.loc=path_to_final_SAIGE_library)
     ```

    Method 2: 

    Open R. Run

    ```
      devtools::install_github("weizhouUMICH/SAIGE")
    ```

### Run SAIGE using a docker image 

Thanks to Juha Karjalainen for sharing the Dockerfile. 
The docker image can be pulled

```
docker pull wzhou88/saige:0.44.6.2
```

Functions can be called
```
step1_fitNULLGLMM.R --help
step2_SPAtests.R --help
createSparseGRM.R --help
```


## Run SAIGE for single-variant association tests and SAIGE-GENE for gene- or region-based tests

Here is a wiki page containg tutorial to run SAIGE and SAIGE-GENE
  https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE
  
### Examples

Example data and script can be found in ./extdata. Run

    bash cmd.sh

to run single-variant and gene-based association tests


## extract effectize sample size v0.44.6.2)
```
  SAIGE_extractNeff.R --help
  bash cmd_extractNeff.sh
```    


# Notes before running jobs

### FAQ can be found  [here](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE#Frequently-asked-questions)

### More notes
1. Since the SPA test always provides close to 0 p-values for variants with MAC < 3, please use at least minMAC = 3 to filter out the results
2. When query is used for bgen files, please make sure there are no duplicate SNP ids in the list
3. If the error message "Error in setgeno(genofile, subSampleInGeno, memoryChunk) :
  vector::_M_range_check", try use a smaller memeoryChunk, such as 2
4. IMPORTANT:In version <= 0.26, for binary traits, BETA is for alt allele and for quantitative traits, BETA is for minor allele 
5. Please note that LOCO only works for autosomal genetic variants. For non-autosomal genetic variants, please leave LOCO=FALSE in step 2.
6. SAIGE-GENE 0.36.3 and 0.36.3.1 now output an effect size for burden tests with the option IsOutputBETASEinBurdenTest in step2. Please note that the magnitude of the effect size is difficult to interpret. 
7. We haven't throughly tested the program on a small sample size. All simulation studies were done using 10,000 samples. Similar to BOLT-LMM, SAIGE uses asymptotic approaches to for feasibility on large samples. Based on our previous real-data analysis, we saw the performance on 3,000 samples were fine. 

# UK Biobank GWAS Results
1. The GWAS results for binary phenotypes in UK Biobank (1,283 phenotypes) using SAIGE are currently available for public download at

https://www.leelabsg.org/resources

Pheweb browser for the UK Biobank results

http://pheweb.sph.umich.edu/SAIGE-UKB/


*This research has been conducted using the UK Biobank Resource under application number 24460.

2. The exome-wide gene-based association results for quantitative traits in UK Biobank (53 traits) using SAIGE-GENE are currently available for public download at

https://www.leelabsg.org/resources

*This research has been conducted using the UK Biobank Resource under application number 45227.


# Log for fixing bugs
* 0.44.6.2 (August-2-2021). add extdata/extractNglmm.R to extract the effective sample size without running Step 1. extdata/cmd_extractNeff.sh has the pipeline. The effective sample size (Nglmm) is differently calculated than the previous versions. 

* 0.44.6.1 (July-16-2021). add the function CCT to perform Cauchy combination to combine multipel tests

* 0.44.6 (July-13-2021). Set --method_to_CollapseUltraRare="absence_or_presence" as default to collpase ultra-rare varaints with MAC <= 10. We call this version SAIGE-GENE+. SAIGE-GENE+ has well controlled type I error rates when the maximum MAF cutoff (maxMAFforGroupTest) is lower than 1%, e.g. 0.01% or 0.1%.  

* 0.44.5 (April-21-2021). 1. re-write code for leave-one-chromosome-out in Step 1 to have more efficient parallel computation in Step 1. 2. Speed up the single-variant association tests when running gene-based tests

* 0.44.2 (March-31-2021) 1.add an option useSparseGRMtoFitNULL to allow for fitting the null model using the sparse GRM and 2. add options to collapse the ultra-rare variants in the set-based tests. --method_to_CollapseUltraRare, --MACCutoff_to_CollapseUltraRare, --DosageCutoff_for_UltraRarePresence

* 0.44.1 (Feb-16-2021) 1. Fixed the error " X %*% Z : non-conformable arguments" for monomorphic variants. 2. merged Jonathon's codes to update savvy to savvy 2.0. For markers in VCF or SAV files without imputation info R2 values, the imputationInfo column will be 1 in the output file, so the markers will not but removed by minInfo 

* 0.44 (January-11-2021) 1. Fixed the error "Phi_ccadj[-indexNeg, -indexNeg]"; 2.  inverse normalization is only performed for quantitative traits; 3. For step 2, bgen input requires the sample file. vcf input does not require a seperate sample file. If sample file is not provided, sample ids will be read from vcf file

* 0.43.3 (January-05-2021)  error "FALis_rewrite_XnonPAR_forMalesSE not found" has been fixed

* 0.43.2 (December-13-2020)  add scripts to calcuate the effectize sample size in Step 1 for binary traits

* 0.43.1. with LOCO=TRUE, remove model results for other chromosomes to save memory usage for Step 2. 

* 0.43 (November-21-2020) Further modify the sparse version of the score test for quantitative traits. This causes slight different assoc tests for variants with MAF < 0.05 for quantitative traits. Set LOCO = TRUE to the default values for step 1 and step 2. In step 2, --chrom needs to be specified for LOCO=TRUE.

* 0.42.1 (September-21-2020) uncomment isSparse=FALSE for quantitative traits. This was commented out for testing in 0.42

* 0.42 (September-16-2020) fix a bug for variance ratio adjustion when account for case-control imbalance for gene-based tests. minMAC is set to 1/(2*N) instead of 0 if is_rewrite_XnonPAR_forMales=TRUE

* 0.41 (August-30-2020) improve the LOCO feature, implement LOCO for gene- and region- based tests (require --chrom to be specified), and with minInfo cutoff, if the input VCF files do not contain info scores, info will be output as NA and markers won't be filtered out. fixed an issue when subsetting pre-calcuated terms (regress X out of G) to drop missing dosages. Use sparse matrices for genotypes/dosages in gene- and region- based tests, so memory usage is dramatically decreased

* 0.39.4 (August-11-2020) use sparse matrix to represent genotype matrix for gene-based tests to save memory

* 0.39.3 (August-6-2020)  add five options --sexCol, --FemaleCode, --FemaleOnly, --MaleCode, --MaleOnly to perform sex-specific Step 1.

* 0.39.2 (July-27-2020)
** add three options --sampleFile_male, --X_PARregion, --is_rewrite_XnonPAR_forMales for chromosome X association tests, in which genotypes/dosages of non-PAR region of males will be multiplied by 2 

* 0.39.1 (July-27-2020)
** add an option --IsOutputlogPforSingle to output log(P) for single-variant assoc tests. v0.39.1 requires SPAtest 3.1.2.  

* 0.39 (May-27-2020)
** fixed an error when conditional analysis is conducted based on vcf input (introduced in 0.38)

* 0.38 (May-4-2020)
** further fixed the bug for output the allele 2 when bgen input with missing dosages was used and missing dosages were dropped. 
** sampleFile is no longer needed if VCF file is used in Step 2
** add --IsOverwriteVarianceRatioFile in step 1 to overwrite the variance ratio file

* 0.37 (May-1-2020)
** fixed an issue with AC values when bgen input is used with missing dosages to be mean imputed (default setting).
 
* 0.36.6 (April-15-2020)
** add an option IsOutputHetHomCountsinCaseCtrl to output the heterozygous and homozygous counts in cases and controls

* 0.36.5.1 (March-29-2020)
** add the option SPAcutoff, If the test statistic lies within the standard deviation cutoff of the mean, p-value based on traditional score test is returned. Otherwise, SPA will be applied. Default value of SPAcutoff is 2 (corresponding p.value.NA 0.05

* 0.36.5 (March-29-2020)
** Fix a typo to extract p.value. 0.36.5: fix an issue for LOCO=TRUE. This issue was introduced when the option minMAFforGRM was introduced.

* 0.36.4.2 (March-20-2020)
** Fix a bug by unlist(p.value), which was introduced in 0.36.4

* 0.36.4.1 (March-18-2020)
** Trying to fix a bug when minMAFforGRM is set and LOCO=TRUE


* 0.36.4 (March-18-2020)
** add an option includeNonautoMarkersforVarRatio in step 1. If TRUE, non-autosomal markers are also used for variance ratio estimation, which will make the algorithm more appropriate for assoc tests for non-autosomal markers; use the new function with sparse sigma for p-values for single variants in gene-based tests; assign AF to be 0 if all samples have missing genotypes or dosages


* 0.36.3.2 (February-25-2020)
** Bug fixed: 1. fixed a bug for gene-based conditioning tests with multiple conditioning markers 2. add codes to re-check markers after dropping samples with missing dosages/genotypes in gene-based tests

* 0.36.3.1 (February-04-2020):
** Note: in v0.36.3.1, uses SPAtest 3.0.2

* 0.36.3 (January-05-2020):
** Note: in v0.36.3, an option IsOutputBETASEinBurdenTest in step 2 is added to output effect sizes for burden tests
Bugs fixed: the header in output files from conditional analysis in gene or reigon-based tests is corrected.  

* 0.36.2 (November-23-2019):
** Note: in v0.36.2, users can specify customized weights for markers in gene- or region-based tests by adding a weight for each marker in the group file

Bugs fixed: 1. The option weights.beta.common is not fully correctly developed, so we make weights.beta.common equal to weights.beta.rare for now. 2. Instead of output NA for SKAT-O p values when the function SKAT:::Met_SKAT_Get_Pvalue failed, output 2*min(SKAT p, Burden p, 0.5).


* 0.36.1 (November-12-2019): 

** Note: in v0.36.1, plain text dosage files are no longer allowed as input in step 2 to get rid of the dependence of the boost_iostream library

Bugs fixed: 1. fixed the freq calculation for mean impute for missing genotypes in  plinkFile 2. Diagonal elements of GRM are now estimated using markers in plinkFile with MAF >= minMAFforGRM 3. Conditional analysis for gene- or region-based test for binary traits is now accounting for case-control imbalance 4. plain dosage files are no longer supported for step 2 so no external boost_iostream library is needed

** minMAFforGRM is added as a parameter in step 0 and 1, so only markers in the plinkFile with MAF >= minMAFforGRM will be used for GRM
** weights.beta.rare, weights.beta.common, weightMAFcutoff, dosageZerodCutoff, IsOutputPvalueNAinGroupTestforBinary, IsAccountforCasecontrolImbalanceinGroupTest are added as new parameters in step 2

* 0.35.8.8: Fixes a matrix inversion issue in the null model and adds an optional argument for the null computation to remove binary covariates with low counts by juhis

* 0.35.8.8 (August-27-2019): Fixes a matrix inversion issue in the null model and adds an optional argument for the null computation to remove binary covariates with low counts by juhis

* 0.35.8.7 (August-15-2019): fixed the bug when there is no covariate specified, added an argument IsOutputNinCaseCtrl for step 2 to allow for output sample sizes in cases and controls for binary traits in the output file, fixed the out of boundary bug for LOCO

* 0.35.8.6 (August-13-2019): fixed the output bug when the genotype matrix has rank 1 for binary phenotypes and add an argument minMAFtoConstructGRM for step 0 and step 1 to allow users to specify the minumum MAF of markers used to construct GRM (default: 1%)

* 0.35.8.5 (June-29-2019): account for case control imbalance for binary traits in gene-based tests

* 0.35.8.3 (May-14-2019): fix a bug in the function getCovM_nopcg, which affected the conditional analysis for binary traits. Merge hyacz/master to use cget to manage superlu 

* 0.35.8.2 (April-16-2019): minor changes include fix error message, change MAC to MAF, add a line to check if the chomosome in plink file is numeric or not, add rsid to the header when input file is bgen

* 0.35.8.1: fix some errors in documentation and the warning message for case-control imbalance of binary traits when running SAIGE-GENE

* 0.35.8 merge changes in the master-gene branch to master

* 0.35.7 merge changes in 0.29.6 and 0.29.7 from master

* 0.35.6 merge in 0.29.5 from master

* 0.35.5 (fix a bug for updating predicted values in the model fit for binary traits. Added a function to create a sparse GRM only for a data set)

* 0.35.3 (this is a clean version for single-variance assoc tests, gene-based tests, and conditional analysis)

* 0.35.2.3 (this version works with the conditonal analysis and gene-based tests)

* 0.29.4.2 (this version works with R-3.5.1)

* 0.29.4 (this version works with R-3.4.4) update SAIGE as a bug for reading vcf and sav files was fixed in the savvy library

* 0.29.3.2 this version works with R-3.5.1

* 0.29.3: update SAIGE step 1 to use the updated R libary SPAtest 3.0.0

* 0.29.2: update SAIGE to use the updated R library SPAtest 3.0.0

* 0.29:
```
1. The colSums() error when there is no covariate has been fixed. 
2. BETA and Tstat are now for the alt allele for both quantitative and binary traits. Note that in version <= 0.26, for binary traits, BETA is for alt allele and for quantitative traits, BETA is for minor allele
3. Options for leave-one-chromosome-out (LOCO), cutoffs for the coefficient of variation (CV) for trace estimates and variance ratio estimates have been added, but these three options have not been extensively tested. CV is mainly for automatically determining whether the number of random markers selected is sufficient or not. If not, the number will be increased until the CV is lower than the specified cutoff.  
```
* 0.26: fixed a bug for the Tstat in the output
* 0.25: allow models with no covariates and GRM contruction using a large number of genetic markers (> 600,000)
* 0.24: centerVariable is no longer needed. QR transformation of the covariate matrix is automatically performed. Supports the dosage files in the VCF,BCF and SAV formats using the SAVVY library 
