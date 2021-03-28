##0.44.3: does not use the bgen C++ library for reading the bgen files in Step 2 any more. re-code the score tests and SPA using C++ to speed up Step 2. This version can only take BGEN files as input for Step 2 for Single-variant association tests. The function SPAGMMATtest_new should be used for Step 2. LOCO=FALSE



#step 0: create a sparse GRM for a data set. This sparse GRM only needs to be created once for each data set, e.g. a biobank,  and can be used for all different phenotypes as long as all tested samples are in the sparse GRM.
Rscript createSparseGRM.R       \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --nThreads=4  \
        --outputPrefix=./output/sparseGRM       \
        --numRandomMarkerforSparseKin=1000      \
        --relatednessCutoff=0.125

##This version can also use sparse GRM to fit the null model in Step 1 by set useSparseGRMtoFitNULL=TRUE, --nThreads=1
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
	--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx    \
        --sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  \
        --nThreads=1	\
        --LOCO=FALSE    \
        --relatednessCutoff=0.125       \
	--useSparseGRMtoFitNULL=TRUE  


#run Step 2, only bgen files can be used and samples with missing doages with be mean imputed
Rscript step2_SPAtests.R        \
        --bgenFile=./input/genotype_100markers.bgen     \
        --bgenFileIndex=./input/genotype_100markers.bgen.bgi    \
        --SAIGEOutputFile=./output/example_binary.SAIGE.bgen.genotype.txt \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --numLinesOutput=10 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE	\
	--LOCO=FALSE
