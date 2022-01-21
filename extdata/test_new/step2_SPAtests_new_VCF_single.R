#!/usr/bin/env Rscript

#options(stringsAsFactors=F, scipen = 999)
options(stringsAsFactors=F)
#library(SAIGE)
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.44.6.5")
#library(SAIGE, lib.loc="../../install_dir/0.43.1")
#library(SAIGE, lib.loc="../../install_dir/0.36.3.3")
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.44.2.b")
print(sessionInfo())

#library(SAIGE, lib.loc="~/projects/Dec2021/install")
library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_old_check/install")

library(optparse)
library(data.table)
library(methods)

SPAGMMATtest(vcfFile="../input/genotype_100markers.vcf.gz",
	     vcfFileIndex="../input/genotype_100markers.vcf.gz.csi",
	     vcfField="GT",
	     SAIGEOutputFile="./new_single_vcf.txt",
	     chrom="1",
	     GMMATmodelFile="../output/example_binary.rda",
	     varianceRatioFile="../output/example_binary.varianceRatio.txt",
	     LOCO=FALSE,
	     is_output_moreDetails=TRUE,
	     numLinesOutput = 100 
	     )

