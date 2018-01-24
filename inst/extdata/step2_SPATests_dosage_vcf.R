options(stringsAsFactors=F, digits=3)

library(SAIGE)
library(optparse)

option_list <- list(
  make_option("--vcfFile", type="character",default="",
    help="path to vcf file. Each line contains dosages for a marker to be tested"),
  make_option("--vcfFileIndex", type="character",default="",
    help="path to vcf index file. Indexed by tabix"),
  make_option("--vcfField", type="character",default="DS",
    help="DS or GT, [default=DS]"),
  make_option("--chrom", type="character",default="0",
    help="chromosome in vcf to be tested. If not specified, all markers in the vcf will be tested"),
  make_option("--start", type="numeric",default=0,
    help="start genome position in the vcf to be tested"),
  make_option("--end", type="numeric",default=0,
    help="end genome position in the vcf to be tested. If not specified, the whole genome will be tested"),
  make_option("--minMAF", type="numeric", default=0,
    help="minimum minor allele frequency for markers to be tested"),
  make_option("--minMAC", type="numeric", default=0,
    help="minimum minor allele count for markers to be tested. Note final threshold will be the greater one between minMAF and minMAC"),
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
  make_option("--SAIGEOutputFile", type="character", default="",
    help="path to the output file containing the SAIGE test results"),
  make_option("--invNormalize", type="logical",default=FALSE,
    help="inverse normalize [default='FALSE']"),
  make_option("--numLinesOutput", type="numeric",default=10000,
    help="output results for every n markers [default=10000]"),
  make_option("--IsOutputAFinCaseCtrl", type="logical",default=FALSE,
    help="whether to output allele frequency in cases and controls for dichotomous traits [default=FALSE]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

try(if(length(which(opt == "")) > 0) stop("Missing arguments"))

covars <- strsplit(opt$covarColList,",")[[1]]
centervars <- strsplit(opt$centerVariables,",")[[1]]


SPAGMMATtest(vcfFile=opt$vcfFile,
             vcfFileIndex=opt$vcfFileIndex,
	     vcfField=opt$vcfField,
             chrom=opt$chrom,
             start=opt$start,
             end=opt$end,
	     sampleFile=opt$sampleFile,
             phenoFile=opt$phenoFile,
             phenoCol=opt$phenoCol,
             traitType=opt$traitType,
             covarColList=covars,
             sampleIDColinphenoFile=opt$sampleIDColinphenoFile,
             centerVariables=centervars,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,
             invNormalize = opt$invNormalize,
             minMAF = opt$minMAF,
	     minMAC = opt$minMAC,
             numLinesOutput = opt$numLinesOutput,
	     IsOutputAFinCaseCtrl = opt$IsOutputAFinCaseCtrl
	
)



