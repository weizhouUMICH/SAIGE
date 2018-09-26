options(stringsAsFactors=F)

## load R libraries
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.29.4.R.3.5.1-perfectSep")
library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.29.4.R.3.5.1")
require(optparse) #install.packages("optparse")

## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile", type="character",default="",
    help="path to plink file for creating the genetic relationship matrix (GRM)"),
  make_option("--phenoFile", type="character", default="",
    help="path to the phenotype file. The file can be either tab or space delimited"),
  make_option("--phenoCol", type="character", default="",
    help="coloumn name for phenotype in the phenotype file, e.g CAD"),
    make_option("--traitType", type="character", default="binary",
    help="binary/quantitative [default=binary]"),
  make_option("--invNormalize", type="logical",default=FALSE,
    help="if quantitative, whether asking SAIGE to perform the inverse normalization for the phenotype [default='FALSE']"),
  make_option("--covarColList", type="character", default="",
    help="list of covariates (comma separated)"),
  make_option("--sampleIDColinphenoFile", type="character", default="IID",
    help="Column name of the IDs in the phenotype file"),
  make_option("--numMarkers", type="integer", default=30,
    help="An integer greater than 0. Number of markers to be used for estimating the variance ratio [default=30]"),
  make_option("--nThreads", type="integer", default=16,
    help="Number of threads (CPUs) to use"),
  make_option("--skipModelFitting", type="logical", default=FALSE,
    help="whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"),
  make_option("--traceCVcutoff", type="numeric", default=1,
    help="The threshold for coefficient of variation (CV) for the trace estimator. Number of runs for trace estimation will be increased until the CV is below the threshold. By default 1. suggested: 0.0025. This option has not been extensively tested."),
  make_option("--ratioCVcutoff", type="numeric", default=1,
    help="The threshold for coefficient of variation (CV) for estimating the variance ratio. The number of randomly selected markers will be increased until the CV is below the threshold. By default 1. suggested 0.001. This option has not been extensively tested."),
  make_option("--LOCO", type="logical", default=FALSE,
    help="Whether to apply the leave-one-chromosome-out (LOCO) approach. By default, FALSE. This option has not been extensively tested."), 
  make_option("--outputPrefix", type="character", default="~/",
    help="path and prefix to the output files [default='~/']")
)

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

covars <- strsplit(opt$covarColList,",")[[1]]

#set seed
set.seed(1)


fitNULLGLMM(plinkFile=opt$plinkFile,
            phenoFile = opt$phenoFile,
            phenoCol = opt$phenoCol,
            traitType = opt$traitType,
            invNormalize = opt$invNormalize,
            covarColList = covars,
            qCovarCol = NULL,
            sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
            nThreads = opt$nThreads,
            numMarkers = opt$numMarkers,
            skipModelFitting = opt$skipModelFitting,
	    traceCVcutoff = opt$traceCVcutoff,
	    ratioCVcutoff = opt$ratioCVcutoff,
            LOCO = opt$LOCO,
            outputPrefix = opt$outputPrefix)


