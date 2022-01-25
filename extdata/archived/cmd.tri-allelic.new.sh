Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers_triallelic.vcf.gz	\
        --vcfFileIndex=./input/genotype_10markers_triallelic.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.gene.triallelic.newtemp.txt \
        --numLinesOutput=1 \
	--groupFile=./input/groupFile_geneBasedtest_simulation.triallelic.txt	\
        --sparseSigmaFile=./output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsSingleVarinGroupTest=TRUE



Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers_triallelic.missing.vcf.gz  \
        --vcfFileIndex=./input/genotype_10markers_triallelic.missing.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.gene.triallelic.missing.newtemp.txt \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest_simulation.triallelic.txt   \
        --sparseSigmaFile=./output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsSingleVarinGroupTest=TRUE	\
	--IsDropMissingDosages=TRUE




	--condition=chr1:4_A/C



Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers.vcf.gz	\
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.gene.triallelic_control.newtemp.txt \
        --numLinesOutput=1 \
	--groupFile=./input/groupFile_geneBasedtest_simulation.triallelic_control.txt	\
        --sparseSigmaFile=./output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsSingleVarinGroupTest=TRUE
