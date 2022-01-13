#!/usr/bin/env Rscript

#options(stringsAsFactors=F, scipen = 999)
options(stringsAsFactors=F)
#library(SAIGE)
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.44.6.5")
#library(SAIGE, lib.loc="../../install_dir/0.43.1")
#library(SAIGE, lib.loc="../../install_dir/0.36.3.3")
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.44.2.b")
print(sessionInfo())

library(SAIGE, lib.loc="~/projects/Dec2021/install")


library(optparse)
library(data.table)
library(methods)


SPAGMMATtest(bedFile="/net/csgspare2/spare1/wenjianb/Data/UKBB_WES_200k/ukb23155_c2_b0_v1.bed",
	     bimFile="/net/csgspare2/spare1/wenjianb/Data/UKBB_WES_200k/ukb23155_c2_b0_v1.bim",
	     famFile="./ukb23155_c2_b0_v1_ref_first.fam", 
             SAIGEOutputFile="./new_single_plink.txt.UKBB.single",
             sampleFile="/net/csgspare2/spare1/wenjianb/Data/UKBB_WES_200k/ukb23155_b0_v1.SAIGE.sample",
	     chrom="2",
	     GMMATmodelFile="/net/csgspare3/snowwhite.archive/zczhao/SAIGE-GENE-UPDATE/realdata/step1/output/UKB_WES_200k_X427.2.rda",
	     varianceRatioFile="/net/csgspare3/snowwhite.archive/zczhao/SAIGE-GENE-UPDATE/realdata/step1/output/UKB_WES_200k_X427.2.varianceRatio.txt",
	     LOCO=FALSE,
	     is_output_moreDetails=TRUE,
	     max_markers_region = 10,
	     numLinesOutput = 1
	     )
	     #AlleleOrder="ref-first",

