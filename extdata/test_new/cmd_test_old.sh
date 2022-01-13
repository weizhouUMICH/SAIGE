##single-variant assoc tests

Rscript step2_SPAtests_old.R        \
        --bgenFile=../input/genotype_100markers.bgen     \
        --bgenFileIndex=../input/genotype_100markers.bgen.bgi    \
        --SAIGEOutputFile=./old_single_bgen.txt	\
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=../input/samplelist.txt \
        --GMMATmodelFile=../output/example_binary.rda \
        --varianceRatioFile=../output/example_binary.varianceRatio.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=TRUE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE	\
	--LOCO=FALSE &> step2_SPAtests_old_single.log 

rm ./new_single_bgen*
Rscript step2_SPAtests_new_BGEN_single.R &> step2_SPAtests_new_BGEN_single.log

rm ./new_single_plink*
Rscript step2_SPAtests_new_PLINK_single.R
diff ./new_single_plink.txt ./new_single_bgen.txt


##single read data
step1Prefix=/net/csgspare3/snowwhite.archive/zczhao/SAIGE-GENE-UPDATE/realdata/step1/output/UKB_WES_200k_X427.2
trait=X427.2
chrom=2
inpath=~/projects/Dec2021/SAIGE/extdata
#Rscript ${inpath}/test_new/step2_SPAtests_old.R        \
Rscript ${inpath}/test_new/step2_SPAtests_currentgithub.R       \
	--bgenFile=/net/csgspare2/spare1/wenjianb/Data/UKBB_WES_200k/ukb23155_c${chrom}_b0_v1_ref_first.bgen    \
        --bgenFileIndex=/net/csgspare2/spare1/wenjianb/Data/UKBB_WES_200k/ukb23155_c${chrom}_b0_v1_ref_first.bgen.bgi   \
        --sampleFile=/net/csgspare2/spare1/wenjianb/Data/UKBB_WES_200k/ukb23155_b0_v1.SAIGE.sample      \
        --GMMATmodelFile=${step1Prefix}.rda     \
        --varianceRatioFile=${step1Prefix}.varianceRatio.txt    \
	--SAIGEOutputFile=./UKB_WES_200k_${trait}_chr${chrom}.OLD_out.txt    \
        --chrom=2 \
        --minMAF=0 \
        --minMAC=0.5 \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=FALSE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE   \
        --LOCO=FALSE &> step2_SPAtests_old_single.log

rm ./new_single_bgen.txt.UKBB.single*
Rscript step2_SPAtests_new_BGEN_UKBB_single.R &> step2_SPAtests_new_BGEN_UKBB_single.log

plink2 --bgen /net/csgspare2/spare1/wenjianb/Data/UKBB_WES_200k/ukb23155_c2_b0_v1_ref_first.bgen ref-first	\
	--export vcf	\
	--out ukb23155_c2_b0_v1_ref_first









##group test
#make group file: from new to old
convertGroupFile_new_to_old.R

MAFarray=( 0.01 0.001 0.05 0.1 )
for maf in "${MAFarray[@]}"
do
   Rscript step2_SPAtests_old.R        \
        --bgenFile=../input/genotype_100markers.bgen     \
        --bgenFileIndex=../input/genotype_100markers.bgen.bgi    \
        --SAIGEOutputFile=./old_group_bgen_${maf}.txt \
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=../input/samplelist.txt \
	--groupFile=./group_old.txt	\
        --GMMATmodelFile=../output/example_binary.rda \
        --varianceRatioFile=../output/example_binary_cate_v2.varianceRatio.txt	\
        --numLinesOutput=2 \
	--maxMAFforGroupTest=${maf}	\
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=TRUE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE   \
        --LOCO=FALSE	\
	--IsSingleVarinGroupTest=TRUE
done

outfile=old_group_bgen.txt
rm ${outfile}
maf=0.01
head -n 1 old_group_bgen_${maf}.txt > ${outfile}
MAFarray=( 0.01 0.001 0.05 0.1 )
for maf in "${MAFarray[@]}"
do
	tail -n +2 old_group_bgen_${maf}.txt | awk -v maf2=$maf '{print maf2"_"$0}' >> ${outfile}
done

##new
#Rscript step2_SPAtests_test_BGEN_group.R

###
grep GENE1_lof group_old.txt > group_old_GENE1_lof.txt
maf=0.01
   Rscript step2_SPAtests_old.R        \
        --bgenFile=../input/genotype_100markers.bgen     \
        --bgenFileIndex=../input/genotype_100markers.bgen.bgi    \
        --SAIGEOutputFile=./old_group_bgen_${maf}_GENE1_lof.txt \
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=../input/samplelist.txt \
        --groupFile=./group_old_GENE1_lof.txt     \
        --GMMATmodelFile=../output/example_binary.rda \
	--varianceRatioFile=../output/example_binary_cate_v2.varianceRatio.txt	\
        --numLinesOutput=2 \
        --maxMAFforGroupTest=${maf}     \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=TRUE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE   \
        --LOCO=FALSE    \
	--IsSingleVarinGroupTest=TRUE	\
        --sparseSigmaFile=../output/example_binary_cate_v2.varianceRatio.txt_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseSigma.mtx &> step2_SPAtests_old.log

        #--IsSingleVarinGroupTest=TRUE &> step2_SPAtests_old.log
	#--varianceRatioFile=./example_binary_cate_v2.varianceRatio.txt.fake	\
        #--varianceRatioFile=../output/example_binary_cate_v2.varianceRatio.txt  \

grep GENE1 ./group_plink.txt > ./group_plink_GENE1.txt
rm ./new_group_bgen_GENE1*
Rscript step2_SPAtests_test_BGEN_group_test.R &> step2_SPAtests_test_BGEN_group_test.log


#########Conditional analysis


Rscript step2_SPAtests_old.R        \
        --bgenFile=../input/genotype_100markers.bgen     \
        --bgenFileIndex=../input/genotype_100markers.bgen.bgi    \
        --SAIGEOutputFile=./old_single_bgen_cond.txt \
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=../input/samplelist.txt \
        --GMMATmodelFile=../output/example_binary.rda \
        --varianceRatioFile=../output/example_binary.varianceRatio.txt.fake \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=TRUE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE   \
        --LOCO=FALSE	\
	--condition=rs18,rs33,rs50 &> step2_SPAtests_old_single_cond.log 

rm ./new_single_bgen_cond*
Rscript step2_SPAtests_new_BGEN_single_cond.R &> step2_SPAtests_new_BGEN_single_cond.log 

rm ./new_single_plink*
Rscript step2_SPAtests_new_PLINK_single.R
diff ./new_single_plink.txt ./new_single_BGEN.txt



##gene based
maf=0.01
   Rscript step2_SPAtests_old.R        \
        --bgenFile=../input/genotype_100markers.bgen     \
        --bgenFileIndex=../input/genotype_100markers.bgen.bgi    \
        --SAIGEOutputFile=./old_group_bgen_${maf}_GENE1_lof_cond.txt \
        --chrom=1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=../input/samplelist.txt \
        --groupFile=./group_old_GENE1_lof.txt     \
        --GMMATmodelFile=../output/example_binary.rda \
        --varianceRatioFile=../output/example_binary_cate_v2.varianceRatio.txt  \
        --numLinesOutput=2 \
        --maxMAFforGroupTest=${maf}     \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsDropMissingDosages=TRUE     \
        --IsOutputHetHomCountsinCaseCtrl=TRUE   \
        --LOCO=FALSE    \
        --IsSingleVarinGroupTest=TRUE   \
	--condition=rs18,rs33,rs50	\
        --sparseSigmaFile=../output/example_binary_cate_v2.varianceRatio.txt_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseSigma.mtx &> step2_SPAtests_old_group_bgen_cond.log

rm ./new_group_bgen_GENE1_cond*
Rscript step2_SPAtests_test_BGEN_group_test_cond.R &> step2_SPAtests_test_BGEN_group_test_cond.log



