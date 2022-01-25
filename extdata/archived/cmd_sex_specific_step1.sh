Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_includenonAutoforvarRatio \
        --nThreads=4    \
        --LOCO=FALSE    \
        --minMAFforGRM=0.01	\
	--FemaleOnly=TRUE	\
	--sexCol=a1	\
	--FemaleCode=1


Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_includenonAutoforvarRatio \
        --nThreads=4    \
        --LOCO=FALSE    \
        --minMAFforGRM=0.01	\
	--MaleOnly=TRUE	\
	--sexCol=a1	\
	--MaleCode=2


Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_includenonAutoforvarRatio \
        --nThreads=4    \
        --LOCO=FALSE    \
        --minMAFforGRM=0.01	\
	--FemaleOnly=TRUE	\
	--FemaleCode=1
