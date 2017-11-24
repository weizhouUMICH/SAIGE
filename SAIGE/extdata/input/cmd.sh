trait=$2

Rscript step2_SPATests.R \	
	--vcfFile=realdata.vcf.gz \
	--vcfIndexFile=realdata.vcf.gz.tbi \ 
	--sampleFile=readdata.vcf.samplelist.txt \ 
        --phenoFile=/net/dumbo/home/zhowei/projects/UKBIOBANK/SAIGE/testInCloud10292017/step1/input/UKB500-phenome-v3-withX.txt \
        --phenoCol=${trait} \
	--traitType=binary \
        --GMMATmodelFile=/net/dumbo/home/zhowei/projects/UKBIOBANK/SPAGMMAT/step1/output/${trait}.rda\
        --varianceRatioFile=/net/dumbo/home/zhowei/projects/UKBIOBANK/SPAGMMAT/step1/output/${trait}.varianceRatio.txt \
        --SPAGMMAToutputFile=${trait}.realdata.vcf.SAIGE.txt
#Rscript step2_SPATests_oldsaige.R \
#at(nrow(sampleListinDosage), " sample IDs are found in sample file\n")
