Rscript step1_fitNULLGLMM_0.41.R             --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly         --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt         --phenoCol=y_binary         --covarColList=x1,x2         --sampleIDColinphenoFile=IID         --traitType=binary                --outputPrefix=./output/example_binary_includenonAutoforvarRatio_0.41         --nThreads=4            --LOCO=FALSE            --minMAFforGRM=0.01


Rscript step2_SPAtests_0.42.2.R        \
        --bgenFile=./input/genotype_10markers.missingness.bgen  \
        --bgenFileIndex=./input/genotype_10markers.missingness.bgen.bgi \
        --SAIGEOutputFile=./output/example_binary.SAIGE.bgen.genotype.dropmissing_0.42.2.txt \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/samplefileforbgen_1000samples.txt \
        --GMMATmodelFile=./output/example_binary_includenonAutoforvarRatio_0.41.rda \
        --varianceRatioFile=./output/example_binary_includenonAutoforvarRatio_0.41.varianceRatio.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=FALSE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE	\
	--LOCO=FALSE

