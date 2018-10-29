
#binary
#For single-variant association tests. 
#Not use sparse GRM and not use categorical variance ratios#
#randomly selected markers with MAC >= 20 are used to estimate the variance ratio
#step 1: fit the NULL glmm 
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=4 \
        --LOCO=FALSE

#step 2: perform the single-variant association tests
Rscript step2_SPAtests.R \
	--vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_1/2


#For gene-based test

#step 1: fit the NULL glmm
#step 1 model result from the single-variant assoc test can be re-used, except that for gene-based tests, variance ratios for multiple MAC categories and a sparse GRM need to be used. If IsSparseKin=TRUE and no sparseSigmaFile and sparseSigmaSampleIDFile are specified, a sparse GRM will be created based on the relatednessCutoff. sparseSigmaFile and sparseSigmaSampleIDFile can be used to specify a pre-calcuated sparse GRM and the sample ids for the sparse GRM. Tested samples would be a subset of samples in the pre-calcuated GRM. 

#To activate the variance ratio estimation based multiple MAC categories, --isCateVarianceRatio=TRUE
#cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude are used to specify the MAC categories
#by default --cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5
#--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,20.5
#corresponding to
#0.5 < MAC <=  1.5
#1.5 < MAC <=  2.5
#2.5 < MAC <=  3.5
#3.5 < MAC <=  4.5
#4.5 < MAC <=  5.5
#5.5 < MAC <=  10.5
#10.5 < MAC <=  20.5
#20.5 < MAC


#with no pre-calcuted sparse GRM
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=4 \
        --LOCO=FALSE	\
        --skipModelFitting=TRUE	\
	--IsSparseKin=TRUE	\
	--relatednessCutoff=0.125	\
	--isCateVarianceRatio=TRUE	

#with pre-calcauted sparse GRM, ./output/example_binary.varianceRatio.txt.pre-cal.sparseGRM.mtx and the corresponding sample ids ./output/example_binary.varianceRatio.txt.sparseGRM.mtx.sample
#if SAIGE was used to generate the pre-cal sparse GRM,  the sample ids for the pre-cal sparse GRM  should be the same as the sampleID in model file(.rda)
#the R code below can be used to extract and write the sample ids for the pre-cal sparse GRM
#
#load("mod.rda")
#write.table(modglmm$sampleID, "./output/example_binary.varianceRatio.txt.sparseGRM.mtx.sample", quote=F, col.names=F, row.names=F)
#
#The following step 1 job will take the pre-cal sparse GRM as an input and output another sparse GRM specifically for the tested phenotype
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=4 \
        --LOCO=FALSE    \
        --skipModelFitting=TRUE \
        --IsSparseKin=TRUE      \
	--sparseSigmaFile=./output/example_binary.varianceRatio.txt.pre-cal.sparseGRM.mtx	\
	--sparseSigmaSampleIDFile=./output/example_binary.varianceRatio.txt.sparseGRM.mtx.sample	\
        --isCateVarianceRatio=TRUE


#Perform gene-based/region-based tests according to the group file specified in groupFile
#IsSingleVarinGroupTest=TRUE is to perform single-variant assoc tests as well for markers included in the gene-based tests
#only vcf, sav, and bgen dosage file formats can be used for gene-based tests
#to perform gene-based tests, --groupFile is used to specify a group file 
#Each line is for one gene/set of
#          variants. The first element is for gene/set name. The rest of
#          the line is for variant ids included in this gene/set. For
#          vcf/sav, the genetic marker ids are in the format
#          chr:pos_ref/alt. For begen, the genetic marker ids should
#          match the ids in the bgen file. Each element in the line is
#          seperated by tab.


Rscript step2_SPAtests.R \
	--vcfFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav \
        --vcfFileIndex=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav.s1r \
        --vcfField=DS \
        --chrom=chr1 \
        --minMAF=0 \
        --minMAC=0.5 \
	--maxMAFforGroupTest=0.01	\
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.gene.txt \
        --numLinesOutput=1 \
	--groupFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.grp.file	\
	--sparseSigmaFile=./output/example_binary.varianceRatio.txt.sparseGRM.mtx	\
        --IsOutputAFinCaseCtrl=TRUE	\
	--IsSingleVarinGroupTest=TRUE		





#quantitative
#step 1: fit the NULL glmm. For quantitative traits, if not normally distributed, inverse normalization needs to be specified to be TRUE --invNormalize=TRUE
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
	--invNormalize=TRUE	\
        --outputPrefix=./output/example_quantitative \
        --nThreads=4 \
        --LOCO=FALSE


Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
	--GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_1/2  #conditional analysis can be performed if a conditioning genetic marker is specified (chr:pos_ref/alt)



##Gene-based tests


