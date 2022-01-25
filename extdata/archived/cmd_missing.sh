Rscript step2_SPAtests.R        \
        --vcfFile=./input/genotype_10markers_triallelic.missing.withVarianceallMissing.vcf.gz \
	--vcfFileIndex=./input/genotype_10markers_triallelic.missing.withVarianceallMissing.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=chr1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.vcf.genotype.dropmissing.withVarianceallMissing.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=TRUE



Rscript step2_SPAtests.R \
	--vcfFile=./input/genotype_10markers_triallelic.missing.withVarianceallMissing.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers_triallelic.missing.withVarianceallMissing.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=chr1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.gene.triallelic.withVarianceallMissing.txt \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest_simulation.triallelic.txt   \
        --sparseSigmaFile=./output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsSingleVarinGroupTest=TRUE
