##################
# For install
install.packages("devtools")
library(devtools)

#install MetaSKAT
devtools::install_github("leeshawn/MetaSKAT") 


install.packages("RcppArmadillo")
install.packages("RcppEigen")
install.packages("BH")

conda install -c anaconda cmake

####################
#
bgenFile = "";
bgenFileIndex = "";
vcfFile = "./input/genotype_10markers.vcf.gz";
vcfFileIndex = "./input/genotype_10markers.vcf.gz.tbi";
vcfField = "GT";
savFile = "";
savFileIndex = "";
sampleFile = "./input/samplelist.txt";
idstoExcludeFile = "";
idstoIncludeFile = "";
rangestoExcludeFile = "";
rangestoIncludeFile = "";
chrom = "1";
start = 1;
end = 250000000;
minMAC = 0.5;
minMAF = 0;
maxMAFforGroupTest = 0.01;
minInfo = 0;
GMMATmodelFile = "./output/example_quantitative.rda";
varianceRatioFile = "./output/example_quantitative_cate.varianceRatio.txt";
SPAcutoff=2;
SAIGEOutputFile = "../example_quantitative.SAIGE.gene.txt";
numLinesOutput = 1;
IsSparse=TRUE;
LOCO=FALSE;
sparseSigmaFile="./output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx";
groupFile="../group_multiSets.txt";
kernel="linear.weighted";
method="optimal.adj";
weights.beta.rare = c(1,25);
weightMAFcutoff = 0.01;
weightsIncludeinGroupFile=FALSE;
r.corr=0;
IsSingleVarinGroupTest = TRUE;
cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5);
cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5);
dosageZerodCutoff = 0.2;
method_to_CollapseUltraRare="absence_or_presence";
MACCutoff_to_CollapseUltraRare = 10;
DosageCutoff_for_UltraRarePresence = 0.5;
function_group_test =c("lof", "missense");
MAF_cutoff=c(0.001,0.01);





Rscript createSparseGRM.R	\
--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
--nThreads=4  \
--outputPrefix=./output/sparseGRM	\
--numRandomMarkerforSparseKin=2000	\
--relatednessCutoff=0.125

Rscript step1_fitNULLGLMM.R     \
--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
--phenoCol=y_quantitative \
--covarColList=x1,x2 \
--sampleIDColinphenoFile=IID \
--traitType=quantitative       \
--invNormalize=TRUE     \
--outputPrefix=./output/example_quantitative \
--outputPrefix_varRatio=./output/example_quantitative_cate	\
--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx   \
--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--nThreads=4 \
--LOCO=FALSE	\
--skipModelFitting=FALSE \
--IsSparseKin=TRUE      \
--isCateVarianceRatio=TRUE	

Rscript step2_SPAtests.R \
--vcfFile=./input/genotype_10markers.vcf.gz \
--vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
--vcfField=GT \
--chrom=1 \
--minMAF=0 \
--minMAC=0.5 \
--maxMAFforGroupTest=0.01       \
--sampleFile=./input/samplelist.txt \
--GMMATmodelFile=./output/example_quantitative.rda \
--varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
--SAIGEOutputFile=./output/example_quantitative.SAIGE.gene.txt \
--numLinesOutput=1 \
--groupFile=./input/groupFile_geneBasedtest_simulation.txt    \
--sparseSigmaFile=./output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx       \
--IsSingleVarinGroupTest=TRUE \
--LOCO=FALSE	\
