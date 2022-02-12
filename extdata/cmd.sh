#check the help info for step 1
#Rscript step1_fitNULLGLMM.R --help

#For Binary traits

#For single-variant association tests.
#randomly selected markers with MAC >= 20 in the plink file are used to estimate the variance ratio
#--minMAFforGRM can be used for specify the minimum MAF of markers in the plink file used for GRM, by default=0.01
#step 1: fit the NULL glmm

# Step 1: Fit the null model using a full GRM

#specify categorical covariates using --qCovarColList
Rscript step1_fitNULLGLMM.R     \
	--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr	\
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10 	\
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=2 

#set all covariates as offset --isCovariateOffset=TRUE. The step 1 will be faster
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=30	\
	--isCovariateOffset=TRUE	\
        --IsOverwriteVarianceRatioFile=TRUE	

#Fit the null model using a sparse GRM
Rscript step1_fitNULLGLMM.R     \
	--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
	--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx	\
	--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt	\
	--useSparseGRMtoFitNULL=TRUE	\
	--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=2    \
        --isCovariateOffset=TRUE	

#Fit the null model using a sparse GRM and do not estiamte variance ratios
Rscript step1_fitNULLGLMM.R     \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=2    \
        --isCovariateOffset=TRUE	\
	--skipVarianceRatioEstimation=TRUE

#Estimate categorical variance ratios (--isCateVarianceRatio=TRUE)
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=64    \
        --isCovariateOffset=TRUE	\
	--isCateVarianceRatio=TRUE	\
	--useSparseGRMforVarRatio=TRUE	\
	--IsOverwriteVarianceRatioFile=TRUE

Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
	--invNormalize=TRUE	\
        --traitType=quantitative        \
        --outputPrefix=./output/example_quantitative \
        --nThreads=64    \
        --isCovariateOffset=TRUE    \
        --isCateVarianceRatio=FALSE      \
        --useSparseGRMforVarRatio=TRUE  \
        --IsOverwriteVarianceRatioFile=TRUE


########Step 2
#single-variant assoc test
Rscript step2_SPAtests.R	\
	--bgenFile=./input/genotype_100markers.bgen    \
	--bgenFileIndex=./input/genotype_100markers.bgen.bgi \
        --SAIGEOutputFile=./genotype_100markers_single_out.txt \
        --chrom=1 \
        --LOCO=TRUE    \
        --AlleleOrder=ref-first \
        --minMAF=0 \
        --minMAC=1 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt   \
        --numLinesOutput=10	\
	--condition=1:13_A/C,1:79_A/C

Rscript step2_SPAtests.R	\
	--bedFile=./input/genotype_100markers.bed	\
        --bimFile=./input/genotype_100markers.bim	\
        --famFile=./input/genotype_100markers.fam	\
        --SAIGEOutputFile=./output/new_single_plink.txt	\
        --chrom=1	\
        --LOCO=TRUE	\
	--minMAF=0 \
        --minMAC=1 \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt   \
        --AlleleOrder=alt-first	\
        --is_output_moreDetails=TRUE	\
	--condition=1:13_A/C,1:79_A/C

Rscript step2_SPAtests.R        \
	--vcfFile=./input/genotype_100markers.vcf.gz	\
	--vcfFileIndex=./input/genotype_100markers.vcf.gz.csi     \
	--vcfField=GT	\
	--SAIGEOutputFile=./output/new_single_vcf.txt \
	--chrom=1       \
        --LOCO=TRUE     \
	--minMAF=0 \
        --minMAC=1 \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt   \
        --is_output_moreDetails=TRUE	\
	--condition=1:13_A/C,1:79_A/C


##group tests
##fit the null model using the full GRM (--useSparseGRMtoFitNULL=FALSE) and estimate the categorical variance ratio (--isCateVarianceRatio=TRUE) using the full GRM and sparse GRM (--useSparseGRMforVarRatio=TRUE)
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_fullGRM \
	--outputPrefix_varRatio=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio	\
        --nThreads=2    \
        --isCovariateOffset=TRUE    \
        --isCateVarianceRatio=TRUE      \
        --useSparseGRMforVarRatio=TRUE



Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative        \
        --outputPrefix=./output/example_quantitative_fullGRM \
        --outputPrefix_varRatio=./output/example_quantitative_fullGRM_sparseGRM_categorical_varRatio  \
        --nThreads=2    \
        --isCovariateOffset=TRUE    \
        --isCateVarianceRatio=TRUE      \
        --useSparseGRMforVarRatio=TRUE

##conduct group test. using --function_group_test to list different annotations and --maxMAFforGroupTest for different max MAF cutoffs.
##By default, SKAT-O, SKAT, and BURDEN tests are performed 
Rscript step2_SPAtests.R        \
        --bgenFile=./input/genotype_100markers.bgen    \
        --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
        --SAIGEOutputFile=./genotype_100markers_bgen_groupTest_out.txt \
        --chrom=1 \
        --LOCO=TRUE    \
        --AlleleOrder=ref-first \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary_fullGRM.rda \
        --varianceRatioFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt   \
        --numLinesOutput=10     \
	--groupFile=./input/group_new_chrposa1a2.txt	\
	--sparseSigmaFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx	\
	--function_group_test="lof,missense;lof,missense;lof;synonymous"	\
	--maxMAFforGroupTest=0.0001,0.001,0.01

	#--groupFile=./input/group_new_snpid.txt	\
##if --r.corr=1 , only BURDEN test is performed
Rscript step2_SPAtests.R        \
        --bgenFile=./input/genotype_100markers.bgen    \
        --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
        --SAIGEOutputFile=./output/genotype_100markers_bgen_groupTest_onlyBURDEN.out.txt \
        --chrom=1 \
        --LOCO=TRUE    \
        --AlleleOrder=ref-first \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary_fullGRM.rda \
        --varianceRatioFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt   \
        --numLinesOutput=10     \
        --groupFile=./input/group_new_chrposa1a2.txt \
        --sparseSigmaFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx       \
        --function_group_test="lof,missense;lof,missense;lof;synonymous"        \
        --maxMAFforGroupTest=0.0001,0.001,0.01	\
	--r.corr=1

##use --condition= to perform conditioning analysis
Rscript step2_SPAtests.R        \
        --bgenFile=./input/genotype_100markers.bgen    \
        --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
        --SAIGEOutputFile=./genotype_100markers_bgen_groupTest_conditional.out.txt \
        --chrom=1 \
        --LOCO=TRUE    \
        --AlleleOrder=ref-first \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary_fullGRM.rda \
        --varianceRatioFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt   \
        --numLinesOutput=10     \
        --groupFile=./input/group_new_chrposa1a2.txt \
        --sparseSigmaFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx       \
        --function_group_test="lof,missense;lof,missense;lof;synonymous"        \
        --maxMAFforGroupTest=0.0001,0.001,0.01	\
	--condition=1:13_A/C,1:79_A/C



##plink file
Rscript step2_SPAtests.R        \
	--bedFile=./input/genotype_100markers.bed       \
        --bimFile=./input/genotype_100markers.bim       \
        --famFile=./input/genotype_100markers.fam       \
        --SAIGEOutputFile=./output/genotype_100markers_plink_groupTest_out.txt \
        --chrom=1 \
        --LOCO=TRUE    \
        --AlleleOrder=alt-first \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary_fullGRM.rda \
        --varianceRatioFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt   \
        --numLinesOutput=10     \
        --groupFile=./input/group_new_snpid.txt \
        --sparseSigmaFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx       \
        --function_group_test="lof,missense;lof,missense;lof;synonymous"        \
        --maxMAFforGroupTest=0.0001,0.001,0.01


##vcf file
Rscript step2_SPAtests.R        \
	--vcfFile=./input/genotype_100markers.vcf.gz    \
        --vcfFileIndex=./input/genotype_100markers.vcf.gz.csi     \
        --vcfField=GT   \
        --SAIGEOutputFile=./output/genotype_100markers_vcf_groupTest_out.txt \
        --chrom=1 \
        --LOCO=TRUE    \
        --AlleleOrder=alt-first \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary_fullGRM.rda \
        --varianceRatioFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt   \
        --numLinesOutput=10     \
        --groupFile=./input/group_new_chrposa1a2.txt \
        --sparseSigmaFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx       \
        --function_group_test="lof,missense;lof,missense;lof;synonymous"        \
        --maxMAFforGroupTest=0.0001,0.001,0.01

\
        
Rscript step2_SPAtests.R        \
        --vcfFile=./input/genotype_10markers.missingness.vcf.gz	\
	--vcfFileIndex=./input/genotype_10markers.missingness.vcf.gz.csi	\
        --vcfField=GT   \
        --SAIGEOutputFile=./output/genotype_100markers_vcf_groupTest_out.txt \
        --chrom=1 \
        --LOCO=TRUE    \
        --AlleleOrder=alt-first \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary_fullGRM.rda \
        --varianceRatioFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt   \
        --numLinesOutput=10     \
        --groupFile=./input/group_new_chrposa1a2.txt \
        --sparseSigmaFile=./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx       \
        --function_group_test="lof,missense;lof,missense;lof;synonymous"        \
        --maxMAFforGroupTest=0.0001,0.001,0.01
