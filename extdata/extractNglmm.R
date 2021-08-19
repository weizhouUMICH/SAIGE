#!/usr/bin/env Rscript

options(stringsAsFactors=F)

## load R libraries
#library(SAIGE, lib.loc="../../install_dir/0.38")
#library(SAIGE, lib.loc="../../install_dir/0.36.6")
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.44.6_Neff")
library(SAIGE)
require(optparse) #install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile", type="character",default="",
    help="Path to plink file for creating the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in he plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios"),
  make_option("--phenoFile", type="character", default="",
    help="Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null GLMM. Please note that covariates to be used in the NULL GLMM need to specified using the argument covarColList"),
  make_option("--phenoCol", type="character", default="",
    help="Coloumn name for phenotype to be tested in the phenotype file, e.g CAD"),
    make_option("--traitType", type="character", default="binary",
    help="binary/quantitative [default=binary]"),
  make_option("--invNormalize", type="logical",default=FALSE,
    help="if quantitative, whether asking SAIGE to perform the inverse normalization for the phenotype [default='FALSE']"),
  make_option("--covarColList", type="character", default="",
    help="List of covariates (comma separated)"),
  make_option("--sampleIDColinphenoFile", type="character", default="IID",
    help="Column name of sample IDs in the phenotype file, e.g. IID"),
  make_option("--tol", type="numeric", default=0.02,
    help="Tolerance for fitting the null GLMM to converge [default=0.02]."),
  make_option("--maxiter", type="integer", default=20,
    help="Maximum number of iterations used to fit the null GLMM [default=20]."),
  make_option("--tolPCG", type="numeric", default=1e-5,
    help="Tolerance for PCG to converge [default=1e-5]."),
  make_option("--maxiterPCG", type="integer", default=500,
    help="Maximum number of iterations for PCG [default=500]."),		
  make_option("--nThreads", type="integer", default=1,
    help="Number of threads (CPUs) to use [default=1]."),
  make_option("--memoryChunk", type="numeric", default=2,
   help="Size (Gb) for each memory chunk [default=2]"),
  make_option("--tauInit", type="character", default="0,0",
   help="Initial values for tau. [default=0,0]"),
  make_option("--outputPrefix", type="character", default="~/",
    help="Path and prefix of the output files [default='~/']"),
  make_option("--outputPrefix_varRatio", type="character", default=NULL,
    help="Path and prefix of the output the variance ratio file [default=NULL]. if NULL, it will be the same as the outputPrefix"),
  make_option("--IsSparseKin", type="logical", default=FALSE,
    help="Whether to use sparse kinship for association test [default='FALSE']"),
  make_option("--sparseGRMFile", type="character", default=NULL,
   help="Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]"),
  make_option("--sparseGRMSampleIDFile", type="character", default=NULL,
   help="Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM [default=NULL]"),
  make_option("--numRandomMarkerforSparseKin", type="integer", default=2000,
    help="Number of randomly selected markers to be used to identify related samples for sparse GRM [default=2000]"),
  make_option("--relatednessCutoff", type="numeric", default=0.125,
    help="Threshold to treat two samples as unrelated if IsSparseKin is TRUE [default=0.125]"),
  make_option("--isCovariateTransform", type="logical", default=TRUE,
    help="Whether use qr transformation on non-genetic covariates [default='TRUE']."),
  make_option("--useSparseSigmaConditionerforPCG", type="logical", default=FALSE,
    help="Whether to sparse GRM to speed up the PCG. Current this option is deactivated. [default='FALSE']."),
  make_option("--useSparseSigmaforInitTau", type="logical", default=FALSE,
    help="Whether to use sparse Sigma to estiamte initial tau [default='FALSE']."),
  make_option("--useSparseGRMtoFitNULL", type="logical", default=FALSE,
    help="Whether to use sparse GRM to fit the null model [default='FALSE']."),	
  make_option("--minMAFforGRM", type="numeric", default=0.01,
    help="Minimum MAF of markers used for GRM"),
  make_option("--minCovariateCount", type="numeric", default=-1,
    help="If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues [default=-1] (no covariates will be excluded)."),
  make_option("--FemaleOnly", type="logical", default=FALSE,
    help="Whether to run Step 1 for females only [default=FALSE]. if TRUE, --sexCol and --FemaleCode need to be specified"), 
  make_option("--MaleOnly", type="logical", default=FALSE,
    help="Whether to run Step 1 for males only [default=FALSE]. if TRUE, --sexCol and --MaleCode need to be specified"),   
  make_option("--sexCol", type="character", default="",
   help="Coloumn name for sex in the phenotype file, e.g Sex"),
  make_option("--FemaleCode", type="character", default="1",
   help="Values in the column for sex in the phenotype file are used for females [default, '1']"),
  make_option("--MaleCode", type="character", default="0",
   help="Values in the column for sex in the phenotype file are used for males [default, '0']"),
  make_option("--noEstFixedEff", type="logical", default=FALSE,
   help="Whether to estimate fixed effect coeffciets. [default, 'FALSE']")
)
#!/usr/bin/env Rscript


## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

covars <- strsplit(opt$covarColList,",")[[1]]

convertoNumeric = function(x,stringOutput){
        y= tryCatch(expr = as.numeric(x),warning = function(w) {return(NULL)})
        if(is.null(y)){
                stop(stringOutput, " is not numeric\n")
        }else{
                cat(stringOutput, " is ", y, "\n")
        }
        return(y)
}

tauInit <- convertoNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")


#set seed
set.seed(1)


getNeff(plinkFile=opt$plinkFile,
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
            memoryChunk = opt$memoryChunk,
            tauInit = tauInit,
            outputPrefix = opt$outputPrefix,
	    outputPrefix_varRatio = opt$outputPrefix_varRatio,
            IsSparseKin = opt$IsSparseKin,
            sparseGRMFile=opt$sparseGRMFile,
            sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
            numRandomMarkerforSparseKin = opt$numRandomMarkerforSparseKin,
            relatednessCutoff = opt$relatednessCutoff,
            isCovariateTransform = opt$isCovariateTransform,
	    useSparseSigmaConditionerforPCG = opt$useSparseSigmaConditionerforPCG,
	    useSparseSigmaforInitTau = opt$useSparseSigmaforInitTau,
	    minMAFforGRM = opt$minMAFforGRM,
	    minCovariateCount=opt$minCovariateCount,
	    sexCol=opt$sexCol,
    	    FemaleCode=opt$FemaleCode,
	    FemaleOnly=opt$FemaleOnly,
	    MaleCode=opt$MaleCode,
	    MaleOnly=opt$MaleOnly,
	    noEstFixedEff=FALSE,
	    useSparseGRMtoFitNULL=opt$useSparseGRMtoFitNULL 
	)	
