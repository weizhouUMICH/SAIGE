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

SPAGMMATtest(bgenFile="./input/genotype_100markers.bgen",
	     bgenFileIndex = "./input/genotype_100markers.bgen.bgi",
	     SAIGEOutputFile="./test_new/bgen.txt",
	     chrom="3",
	     sampleFile="./input/samplelist.txt",
	     GMMATmodelFile="./output/example_binary.rda",
	     varianceRatioFile="./output/example_binary.varianceRatio.txt",
	     LOCO=FALSE,
	     AlleleOrder="ref-first",
	     is_output_moreDetails=TRUE
	     )

