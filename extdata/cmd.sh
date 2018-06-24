#step 1: fit the NULL GLMM
Rscript step1_fitNULLGLMM.R \
	--plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
	--phenoFile=./input/pheno_1000samples.txt \
	--phenoCol=y \
	--covarColList=x1,x2 \
	--sampleIDColinphenoFile=IID \
	--traitType=binary \
	--outputPrefix=./output/example \
	--nThreads=4 \
	--LOCO=TRUE

##step 2: perfrom score test with SPA applied for each marker
######plain text dosage file
 Rscript step2_SPAtests.R \
	--dosageFile=./input/dosage_10markers.txt \
	--dosageFileNrowSkip=1 \
	--dosageFileNcolSkip=6 \
	--dosageFilecolnamesSkip=CHR,SNP,CM,POS,EFFECT_ALLELE,ALT_ALLELE \
	--minMAF=0.0001 \
	--sampleFile=./input/sampleIDindosage.txt \
	--GMMATmodelFile=./output/example.rda \
	--varianceRatioFile=./output/example.varianceRatio.txt \
	--SAIGEOutputFile=./output/example.plainDosage.SAIGE.txt \
	--numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE


#######vcf file (hard call genotypes) 
Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example.SAIGE.vcf.genotype.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE


########vcf file (dosages)
Rscript step2_SPAtests.R \
        --vcfFile=./input/dosage_10markers.vcf.gz \
        --vcfFileIndex=./input/dosage_10markers.vcf.gz.tbi \
        --vcfField=DS \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example.SAIGE.vcf.dosage.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE


########sav file
Rscript step2_SPAtests.R	\
	--savFile=./input/dosage_10markers.sav	\
	--savFileIndex=./input/dosage_10markers.sav.s1r	\
	--minMAF=0.0001 \
        --minMAC=1 \
	--vcfField=DS \
	--chrom=1 \
        --sampleFile=./input/samplefileforbgen_10000samples.txt \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example.SAIGE.sav.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE

##########bgen file
Rscript step2_SPAtests.R \
	--bgenFile=./input/genotype_100markers.bgen \
	--bgenFileIndex=./input/genotype_100markers.bgen.bgi \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/samplefileforbgen_10000samples.txt \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example.SAIGE.bgen.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE
