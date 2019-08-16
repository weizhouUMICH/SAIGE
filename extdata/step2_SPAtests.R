#!/usr/bin/env Rscript

options(stringsAsFactors=F, digits=3, warn=1)
library(SAIGE)
require(optparse)

print(sessionInfo())

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
  make_option("--idstoExcludeFile", type="character",default="",
    help="Path to the file containing variant ids to be excluded from the bgen file. The file does not have a header and each line is for a marker ID."),
  make_option("--idstoIncludeFile", type="character",default="",
    help="Path to the file containing variant ids to be included from the bgen file. The file does not have a header and each line is for a marker ID."),
  make_option("--rangestoExcludeFile", type="character",default="",
    help="Path to the file containing genome regions to be excluded from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header."),
  make_option("--rangestoIncludeFile", type="character",default="",
    help="ath to the file containing genome regions to be included from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header."),
  make_option("--chrom", type="character",default="",
    help="chromosome in vcf to be tested. chrom must be specified for vcf/sav"),
  make_option("--start", type="numeric",default=1,
    help="start genome position in the vcf to be tested"),
  make_option("--end", type="numeric",default=250000000,
    help="end genome position in the vcf to be tested. If not specified, the entire vcf will be tested"),
  make_option("--IsDropMissingDosages", type="logical", default=FALSE,
    help="Whether to drop the samples with missing dosages for testing. If FALSE, the missing dosages with be mean imputed, otherwise, they will be removed before testing. This option only works for bgen, vcf, and sav input."),
  make_option("--minMAF", type="numeric", default=0,
    help="minimum minor allele frequency for markers to be tested"),
  make_option("--minMAC", type="numeric", default=0.5,
    help="minimum minor allele count for markers to be tested. Note final threshold will be the greater one between minMAF and minMAC"),
  make_option("--maxMAFforGroupTest", type="numeric", default=0.5,
    help="max MAF for markers tested in group test"),
   make_option("--minInfo", type="numeric", default=0,
    help="minimum Info for markers to be tested"),
  make_option("--sampleFile", type="character",default="",
    help="Path to the file that contains one column for IDs of samples in the dosage, vcf, sav, or bgen file with NO header"),
  make_option("--GMMATmodelFile", type="character",default="",
    help="path to the input file containing the glmm model(.rda), which is output from step 1"),
  make_option("--varianceRatioFile", type="character",default="",
    help="path to the input file containing the variance ratio, which is output from step 1"),
  make_option("--SAIGEOutputFile", type="character", default="",
    help="path to the output file containing the association test results"),
  make_option("--numLinesOutput", type="numeric",default=10000,
    help="output results for every n markers [default=10000]"),
  make_option("--IsSparse", type="logical",default=TRUE,
    help="Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for binary traits [default=TRUE]."), 
  make_option("--IsOutputAFinCaseCtrl", type="logical",default=FALSE,
    help="whether to output allele frequency in cases and controls for dichotomous traits [default=FALSE]"),
    make_option("--IsOutputNinCaseCtrl", type="logical",default=FALSE,
    help="whether to output sample sizes in cases and controls for dichotomous traits [default=FALSE]"),
  make_option("--LOCO", type="logical", default=FALSE,
    help="Whether to apply the leave-one-chromosome-out option. This option has not been extensively tested."),
  make_option("--condition", type="character",default="",
    help="For conditional analysis. Genetic marker ids (chr:pos_ref/alt if sav/vcf dosage input , marker id if bgen input) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A, Note that currently conditional analysis is only for bgen,vcf,sav input."),
  make_option("--groupFile", type="character", default="",
    help="Path to the file containing the group information for gene-based tests. Each line is for one gene/set of variants. The first element is for gene/set name. The rest of the line is for variant ids included in this gene/set. For vcf/sav, the genetic marker ids are in the format chr:pos_ref/alt. For bgen, the genetic marker ids should match the ids in the bgen file. Each element in the line is seperated by tab."),
  make_option("--sparseSigmaFile", type="character", default="",
    help="path to the output file containing the sparse Sigma output by step 1. from step 1. The suffix of this file is .mtx"),
  make_option("--kernel", type="character", default="linear.weighted",
    help="More options can be seen in the SKAT library"),
  make_option("--method", type="character",default="optimal.adj",
    help="method for gene-based test p-values. More options can be seen in the SKAT library"),
  make_option("--weights.beta", type="character", default="1,25",
    help="More options can be seen in the SKAT librar"),
  make_option("--r.corr", type="character", default=0,
    help="More options can be seen in the SKAT library"),
  make_option("--IsSingleVarinGroupTest",type="logical", default=FALSE,
    help="Whether to perform single-variant assoc tests for genetic markers included in the gene-based tests. By default, FALSE"),
  make_option("--cateVarRatioMinMACVecExclude",type="character", default="0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5",
    help="vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5']"),
  make_option("--cateVarRatioMaxMACVecInclude",type="character", default="1.5,2.5,3.5,4.5,5.5,10.5,20.5",
    help="vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='1.5,2.5,3.5,4.5,5.5,10.5,20.5']"),
  make_option("--singleGClambda",type="numeric", default=1,
    help="GC lambda values that can be used to adjust the gene-based tests results. This value is usually estimated based on the single-variant assoc test results. [default=1]"),
  make_option("--IsOutputPvalueNAinGroupTestforBinary", type="logical",default=FALSE,
    help="whether to output p value if not account for case-control imbalance when performing group test (only for binary traits). [default=FALSE]")	
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#try(if(length(which(opt == "")) > 0) stop("Missing arguments"))
set.seed(1)
dosageFilecolnames <- strsplit(opt$dosageFilecolnamesSkip,",")[[1]]
cateVarRatioMinMACVecExclude <- as.numeric(strsplit(opt$cateVarRatioMinMACVecExclude,",")[[1]])
cateVarRatioMaxMACVecInclude <- as.numeric(strsplit(opt$cateVarRatioMaxMACVecInclude,",")[[1]])
weights.beta <- as.numeric(strsplit(opt$weights.beta,",")[[1]])
print(cateVarRatioMinMACVecExclude)
print(cateVarRatioMaxMACVecInclude)
print(weights.beta)

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
	     sampleFile=opt$sampleFile,
	     idstoExcludeFile=opt$idstoExcludeFile,
	     idstoIncludeFile=opt$idstoIncludeFile,
	     rangestoExcludeFile=opt$rangestoExcludeFile,
	     rangestoIncludeFile=opt$rangestoIncludeFile,			
             chrom=opt$chrom,
             start=opt$start,
             end=opt$end,
	     IsDropMissingDosages=opt$IsDropMissingDosages,	
             minMAF = opt$minMAF,
	     minMAC = opt$minMAC,
	     maxMAFforGroupTest = opt$maxMAFforGroupTest,
	     minInfo = opt$minInfo,
	     GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,	
             numLinesOutput = opt$numLinesOutput,
	     IsSparse=opt$IsSparse,
	     IsOutputAFinCaseCtrl = opt$IsOutputAFinCaseCtrl,
		IsOutputNinCaseCtrl = opt$IsOutputNinCaseCtrl,
	     LOCO = opt$LOCO,
	     condition = opt$condition,
	     sparseSigmaFile=opt$sparseSigmaFile,
	     groupFile = opt$groupFile,
	     kernel=opt$kernel,
	     method=opt$method,
	     weights.beta=weights.beta,
	     r.corr=opt$r.corr,
	     IsSingleVarinGroupTest=opt$IsSingleVarinGroupTest,
             cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude,
             cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude,
	     singleGClambda=opt$singleGClambda,
		IsOutputPvalueNAinGroupTestforBinary=opt$IsOutputPvalueNAinGroupTestforBinary
)



