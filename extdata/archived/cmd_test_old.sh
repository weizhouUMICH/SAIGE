Rscript step2_SPAtests.R        \
        --bgenFile=./input/genotype_100markers.bgen     \
        --bgenFileIndex=./input/genotype_100markers.bgen.bgi    \
        --SAIGEOutputFile=./test_new/old_bgen.txt	\
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=TRUE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE	\
	--LOCO=FALSE


rm ./test_new/plink.txt*
Rscript step2_SPAtests_test_PLINK.R &> step2_SPAtests_test_PLINK.log
rm ./test_new/bgen.txt*
Rscript step2_SPAtests_test.R &> step2_SPAtests_test.log

