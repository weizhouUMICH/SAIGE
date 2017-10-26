##step 1: fit the NULL GLMM
Rscript step1_fitNULLGLMM.R \
	--plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
	--phenoFile=./input/pheno_1000samples.txt \
	--phenoCol=y \
	--covarColList=X2,X3 \
	--sampleIDColinphenoFile=IID \
	--traitType=binary \
	--outputPrefix=./output/example \
	--nThreads=4

##step 2: perfrom score test with SPA applied for each marker
Rscript step2_SPATests.R \
	--dosageFile=./input/dosage_10markers.txt \
	--dosageFileNrowSkip=1 \
	--dosageFileNcolSkip=6 \
	--dosageFilecolnamesSkip=CHR,SNP,CM,POS,EFFECT_ALLELE,ALT_ALLELE \
	--minMAF=0.0001 \
	--sampleFile=./input/sampleIDindosage.txt \
	--phenoFile=./input/pheno_1000samples.txt \
	--phenoCol=y \
        --covarColList=X2,X3 \
        --sampleIDColinphenoFile=IID \
	--GMMATmodelFile=./output/example.rda \
	--varianceRatioFile=./output/example.varianceRatio.txt \
	--SPAGMMAToutputFile=./output/example.SAIGE.txt
