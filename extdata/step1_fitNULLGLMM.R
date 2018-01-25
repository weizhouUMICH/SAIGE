options(stringsAsFactors=F)

## load R libraries
library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/01232018/SAIGE_install")
require(optparse) #install.packages("optparse")

## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile", type="character",default="",
    help="path to plink file to be used for the kinship matrix"),
  make_option("--phenoFile", type="character", default="",
    help="path to the phenotype file, a column 'IID' is required"),
  make_option("--phenoCol", type="character", default="",
    help="coloumn name for phenotype in phenotype file, a column 'IID' is required"),
  make_option("--covarColList", type="character", default="",
    help="list of covariates (comma separated)"),
  make_option("--sampleIDColinphenoFile", type="character", default="IID",
    help="Column name of the IDs in the phenotype file"),
  make_option("--centerVariables", type="character", default="",
    help="Covariates that should be centered (comma separated)"),
  make_option("--skipModelFitting", type="logical", default=FALSE,
    help="skip model fitting, [default='FALSE']"),
  make_option("--traitType", type="character", default="binary",
    help="binary/quantitative [default=binary]"),
  make_option("--outputPrefix", type="character", default="~/",
    help="path to the output files [default='~/']"),
  make_option("--numMarkers", type="integer", default=30,
    help="An integer greater than 0 Number of markers to be used for estimating the variance ratio [default=30]"),
  make_option("--nThreads", type="integer", default=16,
    help="Number of threads"),
  make_option("--invNormalize", type="logical",default=FALSE,
    help="inverse normalize [default='FALSE']")
)
## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

covars <- strsplit(opt$covarColList,",")[[1]]
centervars <- strsplit(opt$centerVariables,",")[[1]]


fitNULLGLMM(plinkFile=opt$plinkFile,
            phenoFile = opt$phenoFile,
            phenoCol = opt$phenoCol,
            traitType = opt$traitType,
            invNormalize = opt$invNormalize,
            covarColList = covars,
            qCovarCol = NULL,
            sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
            centerVariables=centervars,
            tol=0.02,
            maxiter=20,
            tolPCG=1e-5,
            maxiterPCG=500,
            nThreads = opt$nThreads,
            Cutoff = 2,
            numMarkers = opt$numMarkers,
            skipModelFitting = opt$skipModelFitting,
            outputPrefix = opt$outputPrefix)


