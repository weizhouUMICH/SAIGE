Table of Contents
=================

   * [Introduction](#introduction)
   * [Citation](#citation)
   * [How to install and run GATE](#how-to-install-and-run-saige-and-saige-gene)
   * [UK Biobank GWAS Results](#uk-biobank-gwas-results)
   * [Log for fixing bugs](#log-for-fixing-bugs)
   * [Notes for users before running jobs](#notes-for-users-before-running-jobs)

# Introduction


## Current version is 0.40.2

GATE (Genetic Analysis of Time-to-Event phenotypes) is an R package with Scalable and accurate genome-wide association analysis of censored survival data in large scale biobanks using frailty models. 

GATE performs single-variant association tests for time-to-event endpoints. GATE uses uses the saddlepoint approximation (SPA)(mhof, J. P. , 1961; Kuonen, D. 1999; Dey, R. et.al 2017) to account for heavy censoring rates. 

GATE is based on joint work by Rouank Dey and Wei Zhou. 


*This R package is still under development using the SAIGE github repository and will be moved to a new github repository soon. 

# Citation


# How to install and run GATE (Similar to SAIGE and SAIGE-GENE) 

  https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE

The github branch name is **SAIGE_homN_hetN_surv**
The docker image can be found in the docker hub **wzhou88/saige.survival:0.40.2**

# Notes for users before running jobs
* After installation, the package needs to be called as SAIGE (will update)

* Currently, the only difference between this and the regular SAIGE job is that for step1, two additional arguments are used 
```
eventTimeCol = “”,
eventTimeBinSize = 1
```
and 

```
traitType = “survival”
```
eventTimeCol is the column name for the event time, e.g. age of diagnosis
eventTimeBinSize is used to set the bin size for evene times. eventTimeBinSize=1 means the bin size will be 1 and if  eventTimeBinSize is not specified, raw event time values will be used





