options(stringsAsFactors=F)

#library(SAIGE, lib="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/tempInstall")
library(SAIGE)
library(optparse)

option_list <- list(
  make_option("--vcfFile", type="character",default="",
    help="path to vcfFile."),
  make_option("--vcfIndexFile", type="character",default="",
    help="path to vcfIndexFile."), 
  make_option("--UKBMarkerListFile", type="character", default=0,
    help="path to the file for UKB marker list for testing"), 
  make_option("--sampleFile", type="character",default="",
    help="File contains one column for IDs of samples in the bgen file, no header"),
  make_option("--phenoFile", type="character", default="",
    help="path to the phenotype file, a column 'IID' is required"),   
  make_option("--phenoCol", type="character", default="",
    help="path to the phenotype file, a column 'IID' is required"),   
  make_option("--traitType", type="character", default="binary",
    help="binary/quantitative [default=binary]"),
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


SPAGMMATtest(vcfFile=opt$vcfFile,
	     vcfFileIndex=opt$vcfIndexFile,	
             sampleFile=opt$sampleFile,
             phenoFile=opt$phenoFile,
	     phenoCol=opt$phenoCol,
             traitType=opt$traitType,
	     covarColList = c("Sex", "birthYear","PC1.wb","PC2.wb","PC3.wb","PC4.wb"),
             sampleIDColinphenoFile="IID",
             centerVariables="birthYear",
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SPAGMMAToutputFile=opt$SPAGMMAToutputFile,
)


library(SAIGE)
library(SAIGE, lib="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/tempInstall")
vcfFile="realdata.vcf.gz"
vcfFileIndex="realdata.vcf.gz.tbi"
#sampleFile="readdata.vcf.samplelist.txt"
sampleFile="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/gz_quantitaive_08022017/UKbgen/ukbgen.sample"
GMMATmodelFile="/net/dumbo/home/zhowei/projects/UKBIOBANK/SPAGMMAT/step1/output/CAD.rda"
varianceRatioFile="/net/dumbo/home/zhowei/projects/UKBIOBANK/SPAGMMAT/step1/output/CAD.varianceRatio.txt"
SPAGMMAToutputFile="./X041"
phenoCol="CAD"
phenoFile="/net/dumbo/home/zhowei/projects/UKBIOBANK/pheno/pheno.CAD"
covarColList = c("Sex", "birthYear","PC1.wb","PC2.wb","PC3.wb","PC4.wb")
sampleIDColinphenoFile="IID"
centerVariables="birthYear"


SPAGMMATtest(vcfFile=vcfFile,
             vcfFileIndex=vcfFileIndex,
             sampleFile=sampleFile,
             phenoFile=phenoFile,
             phenoCol=phenoCol,
             covarColList = c("Sex", "birthYear","PC1.wb","PC2.wb","PC3.wb","PC4.wb"),
             sampleIDColinphenoFile="IID",
             centerVariables="birthYear",
             GMMATmodelFile=GMMATmodelFile,
             varianceRatioFile=varianceRatioFile,
             SPAGMMAToutputFile=SPAGMMAToutputFile
)




