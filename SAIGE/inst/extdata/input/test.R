#library(SAIGE, lib="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/tempInstall")
library(SAIGE)
vcfFile="realdata.vcf.gz"
vcfFileIndex="realdata.vcf.gz.tbi"
#sampleFile="readdata.vcf.samplelist.txt"
sampleFile="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/gz_quantitaive_08022017/UKbgen/ukbgen.sample"
GMMATmodelFile="/net/dumbo/home/zhowei/projects/UKBIOBANK/SPAGMMAT/step1/output/CAD.rda"
varianceRatioFile="/net/dumbo/home/zhowei/projects/UKBIOBANK/SPAGMMAT/step1/output/CAD.varianceRatio.txt"
SAIGEOutputFile="./X041"
phenoCol="CAD"
phenoFile="/net/dumbo/home/zhowei/projects/UKBIOBANK/pheno/pheno.CAD"
covarColList = c("Sex", "birthYear","PC1.wb","PC2.wb","PC3.wb","PC4.wb")
sampleIDColinphenoFile="IID"
centerVariables="birthYear"

dosageFile="dosage_10markers.txt"
sampleInDosage="sampleIDindosage.txt"
GMMATmodelFile2 =




SPAGMMATtest(dosageFile=dosageFile,
             dosageFileNrowSkip=1,
	     dosageFileNcolSkip=6,
             sampleFile=sampleInDosage,
             phenoFile=phenoFile,
             phenoCol=phenoCol,
             covarColList = c("Sex", "birthYear","PC1.wb","PC2.wb","PC3.wb","PC4.wb"),
             sampleIDColinphenoFile="IID",
             centerVariables="birthYear",
             GMMATmodelFile=GMMATmodelFile,
             varianceRatioFile=varianceRatioFile,
             SAIGEOutputFile=SAIGEOutputFile,
             minMAC = 10,
             numLinesOutput = 10
)






SPAGMMATtest(vcfFile=vcfFile,
             vcfFileIndex=vcfFileIndex,
	     chrom="chr1",
	     start=60000,	
	     end=100000,	
             sampleFile=sampleFile,
             phenoFile=phenoFile,
             phenoCol=phenoCol,
             covarColList = c("Sex", "birthYear","PC1.wb","PC2.wb","PC3.wb","PC4.wb"),
             sampleIDColinphenoFile="IID",
             centerVariables="birthYear",
             GMMATmodelFile=GMMATmodelFile,
             varianceRatioFile=varianceRatioFile,
             SAIGEOutputFile=SAIGEOutputFile,
	     minMAC = 10,
	     numLinesOutput = 10
)
