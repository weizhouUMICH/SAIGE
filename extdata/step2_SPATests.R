options(stringsAsFactors=F)

library(SAIGE)
library(optparse)

option_list <- list(
  make_option("--dosageFile", type="character",default="",
    help="path to dosage file. Each line contains dosages for a marker to be tested"),
  make_option("--dosageFileNrowSkip", type="numeric", default=0,
    help="Number of lines to be skiped in the dosage file"), 
  make_option("--dosageFileNcolSkip", type="numeric", default=5,
    help="Number of columns to be skiped in the dosage file"),    
  make_option("--dosageFilecolnamesSkip", type="character", default="",
    help="list of names for the columns that are skipped (comma separated). These names are not necessarily exactly the same as in the dosage file. They will be used in the header in output file"),
  make_option("--minMAF", type="numeric", default=0,
    help="minimum minor allele frequency for markers to be tested"),
  make_option("--sampleFile", type="character",default="",
    help="File contains one column for IDs of samples in the dosage file"),
  make_option("--phenoFile", type="character", default="",
    help="path to the phenotype file, a column 'IID' is required"),   
  make_option("--phenoCol", type="character", default="",
    help="path to the phenotype file, a column 'IID' is required"),   
  make_option("--traitType", type="character", default="binary",
    help="binary/quantitative [default=binary]"),
  make_option("--covarColList", type="character", default="",
    help="list of covariates (comma separated)"),   
  make_option("--sampleIDColinphenoFile", type="character", default="",
    help="Column name of the IDs in the phenotype file"),   
  make_option("--centerVariables", type="character", default="",
    help="Covariates that should be centered (comma separated)"),   	
  make_option("--GMMATmodelFile", type="character",default="",
    help="path to the input file containing the glmm model"),
  make_option("--varianceRatioFile", type="character",default="",
    help="path to the input file containing the variance ratio"),
  make_option("--SPAGMMAToutputFile", type="character", default="",
    help="path to the output file containing the SPAGMMAT test results"),
  make_option("--invNormalize", type="logical",default=FALSE,
    help="inverse normalize [default='FALSE']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

try(if(length(which(opt == "")) > 0) stop("Missing arguments"))

ids_to_include <- readLines(opt$sampleFile)
covars <- strsplit(opt$covarColList,",")[[1]]
centervars <- strsplit(opt$centerVariables,",")[[1]]
dosageFilecolnames <- strsplit(opt$dosageFilecolnamesSkip,",")[[1]]

SPAGMMATtest(dosageFile=opt$dosageFile,
             dosageFileNrowSkip=opt$dosageFileNrowSkip,
             dosageFileNcolSkip=opt$dosageFileNcolSkip,
	     dosageFilecolnamesSkip=dosageFilecolnames,	
             sampleFile=opt$sampleFile,
             phenoFile=opt$phenoFile,
	     phenoCol=opt$phenoCol,
             traitType=opt$traitType,
             covarColList=covars,
             sampleIDColinphenoFile=opt$sampleIDColinphenoFile,
             centerVariables=centervars,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SPAGMMAToutputFile=opt$SPAGMMAToutputFile,
             ids_to_include=ids_to_include,
             invNormalize = opt$invNormalize,
	     minMAF = 0
)

