#binary
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binarytest \
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
        --GMMATmodelFile=./output/example_binarytest.rda \
        --varianceRatioFile=./output/example_binarytest.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binarytest.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_1/2

Rscript step1_fitNULLGLMM_0.29.4.R.3.5.1.R      \
        --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binarytest_0.29.4.R.3.5.1 \
        --nThreads=4 \
        --LOCO=FALSE

Rscript step2_SPAtests_0.29.4.R.3.5.1.R \
        --vcfFile=./input/dosage_10markers.vcf.gz \
        --vcfFileIndex=./input/dosage_10markers.vcf.gz.tbi \
        --vcfField=DS \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example_binarytest_0.29.4.R.3.5.1.rda \
        --varianceRatioFile=./output/example_binarytest_0.29.4.R.3.5.1.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binarytest_0.29.4.R.3.5.1_cond.SAIGE.vcf.dosage.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \

#quantitative
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
	--invNormalize=TRUE	\
        --outputPrefix=./output/example_quantitativetest \
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
	--GMMATmodelFile=./output/example_quantitativetest_0.29.4.R.3.5.1.rda \
        --varianceRatioFile=./output/example_quantitativetest_0.29.4.R.3.5.1.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitativetest.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_1/2

        #--GMMATmodelFile=./output/example_quantitativetest.rda \
        #--varianceRatioFile=./output/example_quantitativetest.varianceRatio.txt \
        #--SAIGEOutputFile=./output/example_quantitativetest.SAIGE.vcf.genotype.txt_cond \


Rscript step1_fitNULLGLMM_0.29.4.R.3.5.1.R      \
	--plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
        --invNormalize=TRUE     \
        --outputPrefix=./output/example_quantitativetest_0.29.4.R.3.5.1 \
        --nThreads=4 	\
        --LOCO=FALSE

Rscript step2_SPAtests_0.29.4.R.3.5.1.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example_quantitativetest_0.29.4.R.3.5.1.rda \
        --varianceRatioFile=./output/example_quantitativetest_0.29.4.R.3.5.1.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitativetest_0.29.4.R.3.5.1.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_1/2
