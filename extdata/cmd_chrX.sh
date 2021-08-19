Rscript step2_SPAtests.R \
	--vcfFile=./input/genotype_10markers_chrX.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers_chrX.vcf.gz.tbi       \
        --vcfField=GT \
        --chrom=chrX \
        --sampleFile_male=./input/sampleid_males.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.vcf.genotype.chrX.gene.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --is_rewrite_XnonPAR_forMales=TRUE      \
        --X_PARregion=1-9,12-15	\
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest_simulation.chrX.txt    \
        --sparseSigmaFile=./output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsSingleVarinGroupTest=TRUE



Rscript step2_SPAtests.R        \
        --vcfFile=./input/genotype_10markers.missingness_chrX.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.missingness_chrX.vcf.gz.tbi       \
        --vcfField=GT \
        --chrom=chrX \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile_male=./input/sampleid_males.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.vcf.genotype.missingness_chrX.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --is_rewrite_XnonPAR_forMales=TRUE      \
        --X_PARregion=1-9,12-15	\
	--IsDropMissingDosages=TRUE


Rscript step2_SPAtests.R        \
        --vcfFile=./input/genotype_10markers.missingness_chrX.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.missingness_chrX.vcf.gz.tbi       \
        --vcfField=GT \
        --chrom=chrX \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile_male=./input/sampleid_males.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.vcf.genotype.missingness_noDrop_chrX.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --is_rewrite_XnonPAR_forMales=TRUE      \
        --X_PARregion=12-15	\
	--IsDropMissingDosages=FALSE
