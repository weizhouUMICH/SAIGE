options(stringsAsFactors=F)

## load R libraries
library(SAIGE)
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
    help="whether to skip model fitting and only estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"),
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
            outputPrefix = opt$outputPrefix)


