Table of Contents
=================

   * [Introduction](#introduction)
   * [Citation](#citation)
   * [How to install SAIGE and SAIGE-GENE](#how-to-install-and-run-saige-and-saige-gene)
   * [UK Biobank GWAS Results](#uk-biobank-gwas-results)
   * [Log for fixing bugs](#log-for-fixing-bugs)
   * [Notes for users before running jobs](#notes-for-users-before-running-jobs)

# Introduction


## Current version is 0.36.3.2 (Updated on February 25, 2020)

SAIGE is an R package with Scalable and Accurate Implementation of Generalized mixed model (Chen, H. et al. 2016). It accounts for sample relatedness and is feasible for genetic association tests in large cohorts and biobanks (N > 400,000).

SAIGE performs single-variant association tests for binary traits and quantitative taits. For binary traits, SAIGE uses the saddlepoint approximation (SPA)(mhof, J. P. , 1961; Kuonen, D. 1999; Dey, R. et.al 2017) to account for case-control imbalance.

SAIGE-GENE (implemented in the SAIGE R package) performs gene- or region-based association tests (Burde, SKAT, SKAT-O) for binary traits and quantitative traits. Note: SAIGE-GENE accounts for case-control imbalance in gene-based tests (>= 0.35.8.5)

*This R package is still under development

# Citation
The SAIGE manuscript:
Wei Zhou, Jonas B. Nielsen, Lars G. Fritsche, Maiken B. Elvestad, Brooke Wolford, Maoxuan Lin, Kristian Hveem, Hyun Min Kang, Goncalo R. Abecasis, Cristen J. Willer*, Seunggeun Lee* “Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies.” Nature Genetics 50, 1335–1341 (2018)

The SAIGE-GENE pre-print:
https://www.biorxiv.org/content/10.1101/583278v2


# How to install and run SAIGE and SAIGE-GENE


## Install SAIGE/SAIGE-GENE

SAIGE can be installed in 4 ways. The first 3 ways require installing dependencies first. 

* R-3.6.1, gcc >= 5.4.0, cmake 3.14.1, [cget](https://cget.readthedocs.io/en/latest/src/intro.html#installing-cget)
* R packages: "R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "SKAT","MetaSKAT"
* /extdata/install_packages.R can be used to install the R packages


1. The binary install file can be downloaded from [SAIGE releases](https://github.com/weizhouUMICH/SAIGE/releases) or from the master branch. Then using the following command line to install SAIGE

```
R CMD INSTALL SAIGE_XX_R_x86_64-pc-linux-gnu.tar.gz
```
2. Installing using the R library devtools from github
```
devtools::install_github("weizhouUMICH/SAIGE") 

```
3. Installing from the source code. 
```
src_branch=master
repo_src_url=https://github.com/weizhouUMICH/SAIGE
git clone --depth 1 -b $src_branch $repo_src_url

R CMD INSTALL SAIGE
```

4. Using a docker image (Thanks to Juha Karjalainen for sharing the Dockerfile). The docker image can be pulled

```
docker pull wzhou88/saige:0.36.3.2
```
[Dockerfile for creating a docker image for SAIGE](https://github.com/weizhouUMICH/Docker/tree/master/SAIGE)

Functions can be called
```
step1_fitNULLGLMM.R --help
step2_SPAtests.R --help
createSparseGRM.R --help
```

5. Using a conda environment (Thanks to [Wallace(Minxian) Wang](https://github.com/weizhouUMICH/SAIGE/issues/118))

a) create a conda environment using 
 ([conda environment file](https://github.com/weizhouUMICH/SAIGE/blob/master/conda_env/environment-RSAIGE.yml)) 

```
conda env create -f environment-RSAIGE.yml
conda activate RSAIGE
FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"
```

Note: [Here](https://github.com/weizhouUMICH/SAIGE/blob/master/conda_env/createCondaEnvSAIGE_steps.txt) are the steps to create the conda environment file 


b) Using method 3
Open R and install package MetaSKAT

```
install.packages('MetaSKAT')
```

exit R and run command
```
src_branch=master
repo_src_url=https://github.com/weizhouUMICH/SAIGE
git clone --depth 1 -b $src_branch $repo_src_url
R CMD INSTALL SAIGE

```

c) Or using method 2
Open R and run (choose 3 no update any packages): 
```
devtools::install_github("weizhouUMICH/SAIGE")
```

## Run SAIGE for single-variant association tests and SAIGE-GENE for gene- or region-based tests

Here is a wiki page containg tutorial to run SAIGE and SAIGE-GENE
  https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE
  
### Examples

Examplary data and script can be found in ./extdata. Run

    bash cmd.sh

to run single-variant and gene-based association tests


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

# Notes for users before running jobs
1. Since the SPA test always provides close to 0 p-values for variants with MAC < 3, please use at least minMAC = 3 to filter out the results
2. When query is used for bgen files, please make sure there are no duplicate SNP ids in the list
3. If the error message "Error in setgeno(genofile, subSampleInGeno, memoryChunk) :
  vector::_M_range_check", try use a smaller memeoryChunk, such as 2
4. IMPORTANT:In version <= 0.26, for binary traits, BETA is for alt allele and for quantitative traits, BETA is for minor allele 
5. Please note that LOCO only works for autosomal genetic variants. For non-autosomal genetic variants, please leave LOCO=FALSE in step 2.
6. SAIGE-GENE 0.36.3 and 0.36.3.1 now output an effect size for burden tests with the option IsOutputBETASEinBurdenTest in step2. Please note that the magnitude of the effect size is difficult to interpret. 
7. We haven't throughly tested the program on a small sample size. All simulation studies were done using 10,000 samples. Similar to BOLT-LMM, SAIGE uses asymptotic approaches to for feasibility on large samples. Based on our previous real-data analysis, we saw the performance on 3,000 samples were fine. 

