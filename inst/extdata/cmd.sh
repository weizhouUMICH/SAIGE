
#step 1: fit the NULL GLMM
Rscript step1_fitNULLGLMM.R \
	--plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
	--phenoFile=./input/pheno_1000samples.txt \
	--phenoCol=y \
	--covarColList=X2,X3 \
	--sampleIDColinphenoFile=IID \
	--traitType=binary \
	--outputPrefix=./output/example \
	--nThreads=4

##step 2: perfrom score test with SPA applied for each marker with the plain text dosage file 
Rscript step2_SPATests_dosage_plain.R \
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
	--SAIGEOutputFile=./output/example.SAIGE.txt \
	--numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE

#step 2: perfrom score test with SPA applied for each marker with the plain text vcf file
Rscript step2_SPATests_dosage_vcf.R \
	--vcfFile=./input/dosage_10markers.vcf.gz \
	--vcfFileIndex=./input/dosage_10markers.vcf.gz.tbi \
	--vcfField=DS \
	--chrom=1 \
	--minMAF=0.0001 \
	--minMAC=1 \
	--sampleFile=./input/sampleIDindosage.txt \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=X2,X3 \
        --sampleIDColinphenoFile=IID \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example.SAIGE.vcf.dosage.txt \
	--numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE	

Rscript step2_SPATests_dosage_vcf.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=X2,X3 \
        --sampleIDColinphenoFile=IID \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example.SAIGE.vcf.genotype.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE

