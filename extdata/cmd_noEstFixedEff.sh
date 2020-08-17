Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_includenonAutoforvarRatio_noEstFixedEff \
        --nThreads=4    \
        --LOCO=TRUE    \
        --minMAFforGRM=0.01	\
	--noEstFixedEff=TRUE > tst.log_noEstFixedEff

Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_includenonAutoforvarRatio_noEstFixedEff_FALSE \
        --nThreads=4    \
        --LOCO=FALSE    \
        --minMAFforGRM=0.01	\
	--noEstFixedEff=FALSE > tst.log_withEstFixedEff


#quantitative traits
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
        --invNormalize=TRUE     \
        --outputPrefix=./output/example_quantitative_noEstFixedEff \
        --nThreads=4 \
        --LOCO=FALSE    \
        --tauInit=1,0	\
	--noEstFixedEff=TRUE > tst.log_noEstFixedEff_quant.txt

#quantitative traits
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
        --invNormalize=TRUE     \
        --outputPrefix=./output/example_quantitative_noEstFixedEff_FALSE \
        --nThreads=4 \
        --LOCO=FALSE    \
        --tauInit=1,0	\
	--noEstFixedEff=FALSE > tst.log_noEstFixedEff_FALSE_quant.txt


