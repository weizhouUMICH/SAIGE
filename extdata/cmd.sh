#step 1: fit the NULL GLMM


#For single-variant association tests. Not use sparse GRM and not use categorical variance ratios
Rscript step1_fitNULLGLMM.R \
	--plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
	--phenoFile=./input/pheno_1000samples.txt \
	--phenoCol=y \
	--covarColList=x1,x2 \
	--sampleIDColinphenoFile=IID \
	--traitType=quantitative	\
	--invNormalize=TRUE	\
	--outputPrefix=./output/exampletest \
	--nThreads=4 \
	--LOCO=FALSE

Rscript step1_fitNULLGLMM.R \
        --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binarytest \
        --nThreads=4 \
        --LOCO=FALSE


Rscript step1_fitNULLGLMM.R \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly	\
	--phenoFile=./input/Prev_0.1_nfam_1000.pheno_seed_31_tau_1_pheno.txt	\
	--phenoCol=y \
	--covarColList=x1,x2 \
	--sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binarytest \
        --nThreads=4 \
        --LOCO=FALSE
	



Rscript step2_SPAtests.R \
        --vcfFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav \
        --vcfFileIndex=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav.s1r \
        --vcfField=DS \
        --chrom=chr1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=/net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation_08_2018/jobs/SAIGE_SKATO/step2/jobs/samplelist.txt \
        --GMMATmodelFile=./output/example_binarytest.rda \
        --varianceRatioFile=./output/example_binarytest.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binarytest.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
	--condition=chr1:5_A/C



        --condition=1:4_1/2

seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav

#conditional analysis in step 1
Rscript step1_fitNULLGLMM.R \
        --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt_withdosages.txt \
        --phenoCol=y \
        --covarColList=x1,x2,a1 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=./output/example_cond \
        --nThreads=4 \
        --LOCO=FALSE

#conditional analysis in step 2
#        --vcfFile=./input/genotype_10markers.vcf.gz \
 #       --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
#       --savFile=/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE/extdata/input/dosage_10markers.sav      \
#        --savFileIndex=./input/dosage_10markers.sav.s1r \


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

Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/exampletest.rda \
        --varianceRatioFile=./output/exampletest.varianceRatio.txt \
        --SAIGEOutputFile=./output/exampletest.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_1/2

#old package
#binary 
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

Rscript step1_fitNULLGLMM.R	\
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






        --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative\
        --outputPrefix=./output/exampletest_0.29.4.R.3.5.1 \
        --nThreads=4 \
        --LOCO=FALSE




Rscript step1_fitNULLGLMM_0.29.4.R.3.5.1.R	\
	--plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt \
        --phenoCol=y \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative\
        --outputPrefix=./output/exampletest_0.29.4.R.3.5.1 \
        --nThreads=4 \
        --LOCO=FALSE

#old package, conditional analysis in step 1
Rscript step1_fitNULLGLMM_0.29.4.R.3.5.1.R      \
        --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
        --phenoFile=./input/pheno_1000samples.txt_withdosages.txt \
        --phenoCol=y \
        --covarColList=x1,x2,a1 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary\
        --outputPrefix=./output/exampletest_0.29.4.R.3.5.1_cond \
        --nThreads=4 \
        --LOCO=FALSE

Rscript step2_SPAtests_0.29.4.R.3.5.1.R	\
	--vcfFile=./input/dosage_10markers.vcf.gz \
        --vcfFileIndex=./input/dosage_10markers.vcf.gz.tbi \
        --vcfField=DS \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_cond2.SAIGE.vcf.dosage.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_A/C	

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
        --SAIGEOutputFile=./output/example_cond2.SAIGE.vcf.dosage.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_A/C




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
        --IsOutputAFinCaseCtrl=TRUE	\
	--condition=1:4_1/2


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
        --IsOutputAFinCaseCtrl=TRUE	\
	--condition=1:9_1/2,1:11_1/2


########sav file
	--savFile=./input/dosage_10markers.sav	\
Rscript step2_SPAtests.R	\
	--savFile=/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE/extdata/input/dosage_10markers.sav	\
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
