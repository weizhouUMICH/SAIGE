#check the help info for step 1
#Rscript step1_fitNULLGLMM.R --help

#For Binary traits

#For single-variant association tests.
#Not use sparse GRM and not use categorical variance ratios#
#randomly selected markers with MAC >= 20 in the plink file are used to estimate the variance ratio
#--minMAFforGRM can be used for specify the minimum MAF of markers in the plink file used for GRM, by default=0.01
#step 1: fit the NULL glmm

#specify categorical variable using --qCovarColList
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

#set all covariates as offset --noEstFixedEff=TRUE. The step 1 will be faster
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2,a9,a10 \
        --qCovarColList=a9,a10  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=2	\
	--noEstFixedEff=TRUE 


#Fit the null model using a sparse GRM
Rscript step1_fitNULLGLMM.R     \
	--sparseGRMFile


