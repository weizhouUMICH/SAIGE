options(stringsAsFactors=F, digits=3)
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
  make_option("--vcfFile", type="character",default="",
    help="path to vcf file. Each line contains dosages for a marker to be tested"),
  make_option("--vcfFileIndex", type="character",default="",
    help="path to vcf index file. Indexed by tabix"),
  make_option("--vcfField", type="character",default="DS",
    help="DS or GT, [default=DS]"),
  make_option("--bgenFile", type="character",default="",
    help="path to bgen file. Each line contains dosages for a marker to be tested"),
  make_option("--bgenFileIndex", type="character",default="",
    help="Path to the .bgi file (index of the bgen file)"),
  make_option("--savFile", type="character",default="",
    help="Path to the sav file."),
  make_option("--savFileIndex", type="character",default="",
    help="Path to the .s1r file (index of the sav file)."),
  make_option("--chrom", type="character",default="",
    help="chromosome in vcf to be tested. If not specified, all markers in the vcf will be tested"),
  make_option("--start", type="numeric",default=1,
    help="start genome position in the vcf to be tested"),
  make_option("--end", type="numeric",default=250000000,
    help="end genome position in the vcf to be tested. If not specified, the whole genome will be tested"),
  make_option("--minMAF", type="numeric", default=0,
    help="minimum minor allele frequency for markers to be tested"),
  make_option("--minMAC", type="numeric", default=0,
    help="minimum minor allele count for markers to be tested. Note final threshold will be the greater one between minMAF and minMAC"),
  make_option("--sampleFile", type="character",default="",
    help="File contains one column for IDs of samples in the dosage file"),
  make_option("--GMMATmodelFile", type="character",default="",
    help="path to the input file containing the glmm model"),
  make_option("--varianceRatioFile", type="character",default="",
    help="path to the input file containing the variance ratio"),
  make_option("--SAIGEOutputFile", type="character", default="",
    help="path to the output file containing the SAIGE test results"),
  make_option("--numLinesOutput", type="numeric",default=10000,
    help="output results for every n markers [default=10000]"),
  make_option("--IsOutputAFinCaseCtrl", type="logical",default=FALSE,
    help="whether to output allele frequency in cases and controls for dichotomous traits [default=FALSE]"),
  make_option("--LOCO", type="logical", default=FALSE,
    help="Whether to apply the leave-one-chromosome-out option. This option has not been extensively tested.")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#try(if(length(which(opt == "")) > 0) stop("Missing arguments"))

dosageFilecolnames <- strsplit(opt$dosageFilecolnamesSkip,",")[[1]]

SPAGMMATtest(dosageFile=opt$dosageFile,
             dosageFileNrowSkip=opt$dosageFileNrowSkip,
             dosageFileNcolSkip=opt$dosageFileNcolSkip,
             dosageFilecolnamesSkip=dosageFilecolnames,
	     vcfFile=opt$vcfFile,
             vcfFileIndex=opt$vcfFileIndex,
	     vcfField=opt$vcfField,
	     bgenFile=opt$bgenFile,
	     bgenFileIndex=opt$bgenFileIndex,
	     savFile=opt$savFile,
	     savFileIndex=opt$savFileIndex,		
             chrom=opt$chrom,
             start=opt$start,
             end=opt$end,
	     sampleFile=opt$sampleFile,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,
             minMAF = opt$minMAF,
	     minMAC = opt$minMAC,
             numLinesOutput = opt$numLinesOutput,
	     IsOutputAFinCaseCtrl = opt$IsOutputAFinCaseCtrl,
	     LOCO = opt$LOCO
)



