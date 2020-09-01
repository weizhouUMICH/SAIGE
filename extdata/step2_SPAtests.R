#!/usr/bin/env Rscript

#options(stringsAsFactors=F, scipen = 999)
options(stringsAsFactors=F)
library(SAIGE)
#library(SAIGE, lib.loc="../../install_dir/0.38")
#library(SAIGE, lib.loc="../../install_dir/0.36.3.3")
#library(SAIGE, lib.loc="../../install_dir/0.37")
print(sessionInfo())


library(optparse)
library(data.table)
library(methods)

option_list <- list(
  make_option("--vcfFile", type="character",default="",
    help="Path to vcf file."),
  make_option("--vcfFileIndex", type="character",default="",
    help="Path to vcf index file. Indexed by tabix. Path to index for vcf file by tabix, .tbi file by tabix -p vcf file.vcf.gz"),
  make_option("--vcfField", type="character",default="DS",
    help="DS or GT, [default=DS]"),
  make_option("--bgenFile", type="character",default="",
    help="Path to bgen file. Path to bgen file. Currently version 1.2 with 8 bit compression is supported"),
  make_option("--bgenFileIndex", type="character",default="",
    help="Path to the .bgi file (index of the bgen file)"),
  make_option("--savFile", type="character",default="",
    help="Path to the sav file."),
  make_option("--savFileIndex", type="character",default="",
    help="Path to the .s1r file (index of the sav file)."),
  make_option("--idstoExcludeFile", type="character",default="",
    help="Path to a file containing variant ids to be excluded from the bgen file. The file does not have a header and each line is for a marker ID."),
  make_option("--idstoIncludeFile", type="character",default="",
    help="Path to a file containing variant ids to be included from the bgen file. The file does not have a header and each line is for a marker ID."),
  make_option("--rangestoExcludeFile", type="character",default="",
    help="Path to a file containing genome regions to be excluded from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header."),
  make_option("--rangestoIncludeFile", type="character",default="",
    help="Path to a file containing genome regions to be included from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header."),
  make_option("--chrom", type="character",default="0",
    help="string chromosome in vcf to be tested. The string needs to exactly match the chromosome string in the vcf/sav file. For example, '1' does not match 'chr1'. If not specified, all markers in the vcf will be tested. If LOCO is specified, providing chrom will save computation cost."),
  make_option("--start", type="numeric",default=1,
    help="start genome position in the vcf to be tested [default=1]"),
  make_option("--end", type="numeric",default=250000000,
    help="end genome position in the vcf to be tested. If not specified, the whole genome will be tested [default=250000000]"),
  make_option("--IsDropMissingDosages", type="logical", default=FALSE,
    help="Whether to drop samples with missing dosages. If FALSE, the missing dosages with be mean imputed, otherwise, they will be removed before testing. This option only works for bgen, vcf, and sav input."),
  make_option("--minMAF", type="numeric", default=0,
    help="Minimum minor allele frequency for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0]."),
  make_option("--minMAC", type="numeric", default=0,
    help="Minimum minor allele count for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0.5]."),
  make_option("--maxMAFforGroupTest", type="numeric", default=0.5,
    help="Max MAF for markers tested in group test [default=0.5]"),
   make_option("--minInfo", type="numeric", default=0,
    help="Minimum Info for markers to be tested [default=0]"),
  make_option("--sampleFile", type="character",default="",
    help="Path to the file that contains one column for IDs of samples in the dosage file. For version >= 0.38, this file is only needed for bgen files. "),
  make_option("--GMMATmodelFile", type="character",default="",
    help="Path to the input file containing the glmm model, which is output from previous step. Will be used by load()"),
  make_option("--varianceRatioFile", type="character",default="",
    help="Path to the input file containing the variance ratio, which is output from the previous step"),
  make_option("--SAIGEOutputFile", type="character", default="",
    help="Path to the output file containing assoc test results"),
  make_option("--numLinesOutput", type="numeric",default=10000,
    help="Number of  markers to be output each time [default=10000]"),
  make_option("--IsSparse", type="logical",default=TRUE,
    help="Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for binary traits [default=TRUE]."),
  make_option("--SPAcutoff", type="numeric", default=2,
    help=" If the test statistic lies within the standard deviation cutoff of the
mean, p-value based on traditional score test is returned. Default value is 2."), 
  make_option("--IsOutputAFinCaseCtrl", type="logical",default=FALSE,
    help="whether to output allele frequency in cases and controls for dichotomous traits [default=FALSE]"),
  make_option("--IsOutputNinCaseCtrl", type="logical",default=FALSE,
    help="Whether to output sample sizes in cases and controls for dichotomous traits [default=FALSE]"),
  make_option("--IsOutputHetHomCountsinCaseCtrl", type="logical",default=FALSE,
    help="Whether to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns homN_Allele2_cases, hetN_Allele2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls will be output [default=FALSE]"),
  make_option("--LOCO", type="logical", default=FALSE,
    help="Whether to apply the leave-one-chromosome-out option. This option has not been extensively tested."),
  make_option("--condition", type="character",default="",
    help="For conditional analysis. Genetic marker ids (chr:pos_ref/alt if sav/vcf dosage input, marker id if bgen input) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A, Note that currently conditional analysis is only for bgen,vcf,sav input."),
  make_option("--sparseSigmaFile", type="character", default="",
    help="Path to the file containing the sparse Sigma output by step 1. The suffix of this file is .mtx"),
  make_option("--groupFile", type="character", default="",
    help="Path to the file containing the group information for gene-based tests. Each line is for one gene/set of variants. The first element is for gene/set name. The rest of the line is for variant ids included in this gene/set. For vcf/sav, the genetic marker ids are in the format chr:pos_ref/alt. For bgen, the genetic marker ids should match the ids in the bgen file. Each element in the line is seperated by tab."),
  make_option("--kernel", type="character", default="linear.weighted",
    help="More options can be seen in the SKAT library"),
  make_option("--method", type="character",default="optimal.adj",
    help="Method for gene-based test p-values. Methods other than optimal.adj have not been extensively tested. More options can be seen in the SKAT library"),
  make_option("--weights.beta.rare", type="character", default="1,25",
    help="parameters for the beta distribution to weight genetic markers with MAF <= weightMAFcutoff in gene-based tests. More options can be seen in the SKAT library"),
  make_option("--weights.beta.common", type="character", default="1,25",
    help="parameters for the beta distribution to weight genetic markers with MAF > weightMAFcutoff in gene-based tests. More options can be seen in the SKAT library. NOTE: this argument is not fully developed. currently, weights.beta.common is euqal to weights.beta.rare"),
  make_option("--weightMAFcutoff", type="numeric", default="0.01",
    help="See document above for weights.beta.rare and weights.beta.common"),
  make_option("--r.corr", type="character", default=0,
    help="More options can be seen in the SKAT library"),
  make_option("--IsSingleVarinGroupTest",type="logical", default=FALSE,
    help="Whether to perform single-variant assoc tests for genetic markers included in the gene-based tests. By default, FALSE"),
  make_option("--cateVarRatioMinMACVecExclude",type="character", default="0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5",
    help="vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5']"),
  make_option("--cateVarRatioMaxMACVecInclude",type="character", default="1.5,2.5,3.5,4.5,5.5,10.5,20.5",
    help="vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='1.5,2.5,3.5,4.5,5.5,10.5,20.5']"),
  make_option("--dosageZerodCutoff",type="numeric", default=0.2,
    help="In gene- or region-based tests, for each variants with MAC <= 10, dosages <= dosageZerodCutoff with be set to 0. [default=0.2]"),
  make_option("--IsOutputPvalueNAinGroupTestforBinary", type="logical",default=FALSE,
    help="Whether to output p value if not account for case-control imbalance when performing group test (only for binary traits). [default=FALSE]"),
  make_option("--IsAccountforCasecontrolImbalanceinGroupTest", type="logical",default=TRUE,
    help="Whether to account for unbalanced case-control ratios for binary tratis in gene- or region-based tests. [default=TRUE]"),
  make_option("--weightsIncludeinGroupFile", type="logical",default=FALSE,
    help="Whether to specify customized weight for makers in gene- or region-based tests. If TRUE, weights are included in the group file. For vcf/sav, the genetic marker ids and weights are in the format chr:pos_ref/alt;weight. For bgen, the genetic marker ids should match the ids in the bgen filE, e.g. SNPID;weight. Each element in the line is seperated by tab. [default=FALSE]"
),
  make_option("--weights_for_G2_cond",type="character", default=NULL, 
    help="vector of float. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, delimited by comma. e.g. '1,2,3"),
  make_option("--IsOutputBETASEinBurdenTest", type="logical",default=FALSE,
    help="Whether to output effect sizes for burden tests. [default=FALSE]"),
  make_option("--IsOutputlogPforSingle", type="logical",default=FALSE,
    help="Whether to output -log10 p-values instead of p-values. [default=FALSE]"),
  make_option("--analysisType", type="character", default='additive',
    help="'additive' (default), 'recessive' or 'dominant'")
)


parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


convertoNumeric = function(x,stringOutput){
	y= tryCatch(expr = as.numeric(x),warning = function(w) {return(NULL)})
	if(is.null(y)){
		stop(stringOutput, " is not numeric\n")
	}else{
		cat(stringOutput, " is ", y, "\n")
	}
	return(y)	
}


#weights.beta.rare <- as.numeric(strsplit(opt$weights.beta.rare,",")[[1]])
weights.beta.rare <- convertoNumeric(x=strsplit(opt$weights.beta.rare,",")[[1]], "weights.beta.rare")
weights.beta.common <- convertoNumeric(x=strsplit(opt$weights.beta.common,",")[[1]], "weights.beta.common")
if(sum(weights.beta.common!=weights.beta.rare) > 0){stop("weights.beta.common option is not functioning, so weights.beta.common needs to be equal to weights.beta.rare")}

cateVarRatioMinMACVecExclude <- convertoNumeric(x=strsplit(opt$cateVarRatioMinMACVecExclude,",")[[1]], "cateVarRatioMinMACVecExclude")
cateVarRatioMaxMACVecInclude <- convertoNumeric(x=strsplit(opt$cateVarRatioMaxMACVecInclude,",")[[1]], "cateVarRatioMaxMACVecInclude")
if(is.null(opt$weights_for_G2_cond)){
	weights_for_G2_cond=NULL
}else{
	weights_for_G2_cond <- convertoNumeric(x=strsplit(opt$weights_for_G2_cond,",")[[1]], "weights_for_G2_cond")
}


#try(if(length(which(opt == "")) > 0) stop("Missing arguments"))

SPAGMMATtest(vcfFile=opt$vcfFile,
             vcfFileIndex=opt$vcfFileIndex,
             vcfField=opt$vcfField,
             bgenFile=opt$bgenFile,
             bgenFileIndex=opt$bgenFileIndex,
             savFile=opt$savFile,
             savFileIndex=opt$savFileIndex,
	     idstoExcludeFile=opt$idstoExcludeFile,
	     idstoIncludeFile=opt$idstoIncludeFile,
             rangestoExcludeFile=opt$rangestoExcludeFile,
             rangestoIncludeFile=opt$rangestoIncludeFile,  	  
             chrom=opt$chrom,
             start=opt$start,
             end=opt$end,
             IsDropMissingDosages=opt$IsDropMissingDosages,
             sampleFile=opt$sampleFile,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,
             minMAF = opt$minMAF,
             minMAC = opt$minMAC,
             numLinesOutput = opt$numLinesOutput,
             IsOutputAFinCaseCtrl = opt$IsOutputAFinCaseCtrl,
	     IsOutputNinCaseCtrl = opt$IsOutputNinCaseCtrl,
             condition = opt$condition,
	     maxMAFforGroupTest = opt$maxMAFforGroupTest,
	     groupFile = opt$groupFile,
	     sparseSigmaFile = opt$sparseSigmaFile,
	     minInfo=opt$minInfo,
	     cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
             cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
	     IsSingleVarinGroupTest = opt$IsSingleVarinGroupTest,
	     dosageZerodCutoff=opt$dosageZerodCutoff,
	     IsOutputPvalueNAinGroupTestforBinary=opt$IsOutputPvalueNAinGroupTestforBinary,
	     IsAccountforCasecontrolImbalanceinGroupTest=opt$IsAccountforCasecontrolImbalanceinGroupTest,
	     method=opt$method,
	     kernel=opt$kernel,
	     weights.beta.rare=weights.beta.rare,
	     weights.beta.common=weights.beta.common,
	     weightMAFcutoff=opt$weightMAFcutoff,
	     r.corr=opt$r.corr,
	     weightsIncludeinGroupFile=opt$weightsIncludeinGroupFile,
	     weights_for_G2_cond=weights_for_G2_cond,
		IsOutputBETASEinBurdenTest=opt$IsOutputBETASEinBurdenTest,
	SPAcutoff=opt$SPAcutoff,
	     LOCO=opt$LOCO,
	     IsOutputlogPforSingle=opt$IsOutputlogPforSingle,
	     analysisType = opt$analysisType
)
