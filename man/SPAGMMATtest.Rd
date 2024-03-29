% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAIGE_SPATest.R
\name{SPAGMMATtest}
\alias{SPAGMMATtest}
\title{Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.}
\usage{
SPAGMMATtest(
  bgenFile = "",
  bgenFileIndex = "",
  vcfFile = "",
  vcfFileIndex = "",
  vcfField = "DS",
  savFile = "",
  savFileIndex = "",
  sampleFile = "",
  idstoExcludeFile = "",
  idstoIncludeFile = "",
  rangestoExcludeFile = "",
  rangestoIncludeFile = "",
  chrom = "",
  start = 1,
  end = 2.5e+08,
  IsDropMissingDosages = FALSE,
  minMAC = 0.5,
  minMAF = 0,
  maxMAFforGroupTest = 0.5,
  minInfo = 0,
  GMMATmodelFile = "",
  varianceRatioFile = "",
  SPAcutoff = 2,
  SAIGEOutputFile = "",
  numLinesOutput = 10000,
  IsSparse = TRUE,
  IsOutputAFinCaseCtrl = FALSE,
  IsOutputHetHomCountsinCaseCtrl = FALSE,
  IsOutputNinCaseCtrl = FALSE,
  IsOutputlogPforSingle = FALSE,
  LOCO = TRUE,
  condition = "",
  sparseSigmaFile = "",
  groupFile = "",
  kernel = "linear.weighted",
  method = "optimal.adj",
  weights.beta.rare = c(1, 25),
  weights.beta.common = c(1, 25),
  weightMAFcutoff = 0.01,
  weightsIncludeinGroupFile = FALSE,
  weights_for_G2_cond = NULL,
  r.corr = 0,
  IsSingleVarinGroupTest = TRUE,
  cateVarRatioMinMACVecExclude = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
  cateVarRatioMaxMACVecInclude = c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
  dosageZerodCutoff = 0.2,
  IsOutputPvalueNAinGroupTestforBinary = FALSE,
  IsAccountforCasecontrolImbalanceinGroupTest = TRUE,
  IsOutputBETASEinBurdenTest = FALSE,
  IsOutputMAFinCaseCtrlinGroupTest = FALSE,
  X_PARregion = "60001-2699520,154931044-155270560",
  is_rewrite_XnonPAR_forMales = FALSE,
  sampleFile_male = "",
  method_to_CollapseUltraRare = "absence_or_presence",
  MACCutoff_to_CollapseUltraRare = 10,
  DosageCutoff_for_UltraRarePresence = 0.5
)
}
\arguments{
\item{bgenFile}{character. Path to bgen file. Currently version 1.2 with 8 bit compression is supported}

\item{bgenFileIndex}{character. Path to the .bgi file (index of the bgen file)}

\item{vcfFile}{character. Path to vcf file}

\item{vcfFileIndex}{character. Path to index for vcf file by tabix, ".tbi" by "tabix -p vcf file.vcf.gz"}

\item{vcfField}{character. genotype field in vcf file to use. "DS" for dosages or "GT" for genotypes. By default, "DS".}

\item{savFile}{character. Path to sav file}

\item{savFileIndex}{character. Path to index for sav file .s1r}

\item{sampleFile}{character. Path to the file that contains one column for IDs of samples in the bgen file with NO header}

\item{idstoExcludeFile}{character. Path to the file containing variant ids to be excluded from the bgen file. The file does not have a header and each line is for a marker ID.}

\item{idstoIncludeFile}{character. Path to the file containing variant ids to be included from the bgen file. The file does not have a header and each line is for a marker ID.}

\item{rangestoExcludeFile}{character. Path to the file containing genome regions to be excluded from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header}

\item{rangestoIncludeFile}{character. Path to the file containing genome regions to be included from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header}

\item{chrom}{character. string for the chromosome to include from vcf file. Required for vcf file. Note: the string needs to exactly match the chromosome string in the vcf/sav file. For example, "1" does not match "chr1". If LOCO is specified, providing chrom will save computation cost}

\item{start}{numeric. start genome position to include from vcf file. By default, 1}

\item{end}{numeric. end genome position to include from vcf file. By default, 250000000}

\item{IsDropMissingDosages}{logical. whether to drop missing dosages (TRUE) or to mean impute missing dosages (FALSE). By default, FALSE. This option only works for bgen, vcf, and sav input.}

\item{minMAC}{numeric. Minimum minor allele count of markers to test. By default, 0.5. The higher threshold between minMAC and minMAF will be used}

\item{minMAF}{numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used}

\item{maxMAFforGroupTest}{numeric. Maximum minor allele frequency of markers to test in group test. By default 0.5.}

\item{minInfo}{numeric. Minimum imputation info of markers to test. By default, 0. This option only works for bgen, vcf, and sav input}

\item{GMMATmodelFile}{character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()}

\item{varianceRatioFile}{character. Path to the input file containing the variance ratio, which is output from the previous step}

\item{SPAcutoff}{by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)}

\item{SAIGEOutputFile}{character. Path to the output file containing assoc test results}

\item{numLinesOutput}{numeric. Number of  markers to be output each time. By default, 10000}

\item{IsSparse}{logical. Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for dichotomous traits. By default, TRUE}

\item{IsOutputAFinCaseCtrl}{logical. Whether to output allele frequency in cases and controls. By default, FALSE}

\item{IsOutputHetHomCountsinCaseCtrl}{logical. Whether to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns "homN_Allele2_cases", "hetN_Allele2_cases", "homN_Allele2_ctrls", "hetN_Allele2_ctrls" will be output.}

\item{IsOutputNinCaseCtrl}{logical. Whether to output sample sizes in cases and controls. By default, FALSE}

\item{IsOutputlogPforSingle}{logical. Whether to output log(Pvalue) for single-variant assoc tests. By default, FALSE. If TRUE, the log(Pvalue) instead of original P values will be output}

\item{LOCO}{logical. Whether to apply the leave-one-chromosome-out option. By default, TRUE}

\item{condition}{character. For conditional analysis. Genetic marker ids (chr:pos_ref/alt if sav/vcf dosage input , marker id if bgen input) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A, Note that currently conditional analysis is only for bgen,vcf,sav input.}

\item{sparseSigmaFile}{character. Path to the file containing the sparseSigma from step 1. The suffix of this file is ".mtx".}

\item{groupFile}{character. Path to the file containing the group information for gene-based tests. Each line is for one gene/set of variants. The first element is for gene/set name. The rest of the line is for variant ids included in this gene/set. For vcf/sav, the genetic marker ids are in the format chr:pos_ref/alt. For bgen, the genetic marker ids should match the ids in the bgen file. Each element in the line is seperated by tab.}

\item{kernel}{character. For gene-based test. By default, "linear.weighted". More options can be seen in the SKAT library}

\item{method}{character. method for gene-based test p-values. By default, "optimal.adj". More options can be seen in the SKAT library}

\item{weights.beta.rare}{vector of numeric. parameters for the beta distribution to weight genetic markers with MAF <= weightMAFcutoff in gene-based tests.By default, "c(1,25)". More options can be seen in the SKAT library}

\item{weights.beta.common}{vector of numeric. parameters for the beta distribution to weight genetic markers with MAF > weightMAFcutoff in gene-based tests.By default, "c(1,25)". More options can be seen in the SKAT library. NOTE: this argument is not fully developed. currently, weights.beta.common is euqal to weights.beta.rare}

\item{weightMAFcutoff}{numeric. Between 0 and 0.5. See document above for weights.beta.rare and weights.beta.common. By default, 0.01}

\item{weightsIncludeinGroupFile}{logical. Whether to specify customized weight for makers in gene- or region-based tests. If TRUE, weights are included in the group file. For vcf/sav, the genetic marker ids and weights are in the format chr:pos_ref/alt;weight. For bgen, the genetic marker ids should match the ids in the bgen filE, e.g. SNPID;weight. Each element in the line is seperated by tab. By default, FALSE}

\item{weights_for_G2_cond}{vector of float. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, delimited by comma. By default, "c(1,2)"}

\item{r.corr}{numeric. bewteen 0 and 1. parameters for gene-based tests.  By default, 0.  More options can be seen in the SKAT library}

\item{IsSingleVarinGroupTest}{logical. Whether to perform single-variant assoc tests for genetic markers included in the gene-based tests. By default, FALSE}

\item{cateVarRatioMinMACVecExclude}{vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used}

\item{cateVarRatioMaxMACVecInclude}{vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(1.5,2.5,3.5,4.5,5.5,10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used}

\item{dosageZerodCutoff}{numeric. In gene- or region-based tests, for each variants with MAC <= 10, dosages <= dosageZerodCutoff with be set to 0. By default, 0.2.}

\item{IsOutputPvalueNAinGroupTestforBinary}{logical. In gene- or region-based tests for binary traits. if IsOutputPvalueNAinGroupTestforBinary is TRUE, p-values without accounting for case-control imbalance will be output. By default, FALSE}

\item{IsAccountforCasecontrolImbalanceinGroupTest}{logical. In gene- or region-based tests for binary traits. If IsAccountforCasecontrolImbalanceinGroupTest is TRUE, p-values after accounting for case-control imbalance will be output. By default, TRUE}

\item{IsOutputBETASEinBurdenTest}{logical. Output effect size (BETA and SE) for burden tests. By default, FALSE}

\item{IsOutputMAFinCaseCtrlinGroupTest}{logical. Whether to output minor allele frequency in cases and controls in set-based tests By default, FALSE}

\item{X_PARregion}{character. ranges of (pseudoautosomal) PAR region on chromosome X, which are seperated by comma and in the format start:end. By default: '60001-2699520,154931044-155260560' in the UCSC build hg19. For males, there are two X alleles in the PAR region, so PAR regions are treated the same as autosomes. In the NON-PAR regions (outside the specified PAR regions on chromosome X), for males, there is only one X allele. If is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in the NON-PAR regions on chromosome X will be multiplied by 2.}

\item{is_rewrite_XnonPAR_forMales}{logical. Whether to rewrite gentoypes or dosages of variants in the NON-PAR regions on chromosome X for males (multiply by 2). By default, FALSE. Note, only use is_rewrite_XnonPAR_forMales=TRUE when the specified VCF or Bgen file only has variants on chromosome X. When is_rewrite_XnonPAR_forMales=TRUE, the program does not check the chromosome value by assuming all variants are on chromosome X}

\item{sampleFile_male}{character. Path to the file containing one column for IDs of MALE samples in the bgen or vcf file with NO header. Order does not matter}

\item{method_to_CollapseUltraRare}{character. Method to collpase the ultra rare variants in the set-based association tests. This argument can be 'absence_or_presence', 'sum_geno', or ''. absence_or_presence:  For the resulted collpased marker, any individual having DosageCutoff_for_UltraRarePresence <= dosage < 1+DosageCutoff_for_UltraRarePresence for any ultra rare variant has 1 in the genotype vector, having dosage >= 1+DosageCutoff_for_UltraRarePresence for any ultra rare variant has 2 in the genotype vector, otherwise 0. sum_geno: Ultra rare variants with MAC <=  MACCutoff_to_CollapseUltraRare will be collpased for set-based tests in the 'sum_geno' way and the resulted collpased marker's genotype equals weighted sum of the genotypes of all ultra rare variants. NOTE: this option sum_geno currently is NOT active. By default, "absence_or_presence".}

\item{MACCutoff_to_CollapseUltraRare}{numeric. MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10.}

\item{DosageCutoff_for_UltraRarePresence}{numeric. Dosage cutoff to determine whether the ultra rare variants are absent or present in the samples. Dosage >= DosageCutoff_for_UltraRarePresence indicates the varaint in present in the sample. 0< DosageCutoff_for_UltraRarePresence <= 2. By default, 0.5.}
}
\value{
SAIGEOutputFile
}
\description{
Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
}
