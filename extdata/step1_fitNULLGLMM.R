options(stringsAsFactors=F)

## load R libraries
#library(SAIGE)
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.35.2.mmSKAT.debugged.R-3.5.1.test2_subsetSparseSigma_speedup_test2")
library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.35.3.2")

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
  make_option("--tol", type="numeric", default=0.02,
    help="The tolerance for fitting the null GLMMM to converge. By default, 0.02."),
  make_option("--tolPCG", type="numeric", default=1e-5,
    help="The tolerance for PCG to converge"),
  make_option("--maxiterPCG", type="integer", default=500,
    help="The maximum number of iterations for PCG. By default, 500"),		
    make_option("--nThreads", type="integer", default=16,
    help="Number of threads (CPUs) to use"),
  make_option("--SPAcutoff", type="numeric", default=2,
    help="The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2."),
  make_option("--numMarkers", type="integer", default=30,
    help="An integer greater than 0. Number of markers to be used for estimating the variance ratio [default=30]"),
  make_option("--skipModelFitting", type="logical", default=FALSE,
    help="whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"),
  make_option("--memoryChunk", type="numeric", default=2,
   help="The size (Gb) for each memory chunk. By default, 2"),
  make_option("--tauInit", type="character", default="0,0",
   help="Unitial values for tau. [default=0,0]"),
  make_option("--LOCO", type="logical", default=FALSE,
    help="Whether to apply the leave-one-chromosome-out (LOCO) approach. By default, FALSE. This option has not been extensively tested."),
  make_option("--traceCVcutoff", type="numeric", default=0.0025,
    help="The threshold for coefficient of variation (CV) for the trace estimator. Number of runs for trace estimation will be increased until the CV is below the threshold. By default 0.0025."),
  make_option("--ratioCVcutoff", type="numeric", default=0.001,
    help="The threshold for coefficient of variation (CV) for estimating the variance ratio. The number of randomly selected markers will be increased until the CV is below the threshold. By default 0.001. "),
  make_option("--outputPrefix", type="character", default="~/",
    help="path and prefix to the output files [default='~/']"),
  make_option("--IsSparseKin", type="logical", default=FALSE,
    help="Whether to use sparse kinship for association test [default='FALSE']"),
  make_option("--sparseSigmaFile", type="character", default="",
   help="Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed [default='']"),
  make_option("--sparseSigmaSampleIDFile", type="character", default="",
   help="Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM [default='']"),
  make_option("--numRandomMarkerforSparseKin", type="integer", default=500,
    help="number of randomly selected markers (MAF >= 1%) to be used to identify related samples for sparse GRM [default=500]"),
  make_option("--isCateVarianceRatio", type="logical", default=FALSE,
    help="Whether to estimate variance ratio based on different MAC categories. If yes, variance ratio will be estiamted for multiple MAC categories corresponding to cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude. Currently, if isCateVarianceRatio=TRUE, then LOCO=FALSE [default=FALSE]"),
  make_option("--relatednessCutoff", type="numeric", default=0.125,
    help="The threshold to treat two samples as unrelated if IsSparseKin is TRUE [default=0.125]"),
  make_option("--cateVarRatioMinMACVecExclude",type="character", default="0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5",
    help="vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5']"),
  make_option("--cateVarRatioMaxMACVecInclude",type="character", default="1.5,2.5,3.5,4.5,5.5,10.5,20.5",
    help="vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='1.5,2.5,3.5,4.5,5.5,10.5,20.5']"),    
  make_option("--isCovariateTransform", type="logical", default=TRUE,
    help="Whether use qr transformation on non-genetic covariates [default='TRUE']."),
  make_option("--isDiagofKinSetAsOne", type="logical", default=FALSE,
    help="Whether to set the diagnal elements in GRM to be 1 [default='FALSE'].")
)



## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

covars <- strsplit(opt$covarColList,",")[[1]]
tauInit <- as.numeric(strsplit(opt$tauInit, ",")[[1]])
cateVarRatioMinMACVecExclude <- as.numeric(strsplit(opt$cateVarRatioMinMACVecExclude,",")[[1]])
cateVarRatioMaxMACVecInclude <- as.numeric(strsplit(opt$cateVarRatioMaxMACVecInclude,",")[[1]])

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
	    tol=opt$tol,
	    maxiter=opt$maxiter,
            tolPCG=opt$tolPCG,
            maxiterPCG=opt$maxiterPCG,
            nThreads = opt$nThreads,
            SPAcutoff = opt$SPAcutoff,
            numMarkers = opt$numMarkers,
            skipModelFitting = opt$skipModelFitting,
            memoryChunk = opt$memoryChunk,
            tauInit = tauInit,
            LOCO = opt$LOCO,
            traceCVcutoff = opt$traceCVcutoff,
            ratioCVcutoff = opt$ratioCVcutoff,
            outputPrefix = opt$outputPrefix,
            IsSparseKin = opt$IsSparseKin,
            sparseSigmaFile=opt$sparseSigmaFile,
            sparseSigmaSampleIDFile=opt$sparseSigmaSampleIDFile,
            numRandomMarkerforSparseKin = opt$numRandomMarkerforSparseKin,
            relatednessCutoff = opt$relatednessCutoff,
            isCateVarianceRatio = opt$isCateVarianceRatio,
            cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
            cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
            isCovariateTransform = opt$isCovariateTransform,
            isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne)	
