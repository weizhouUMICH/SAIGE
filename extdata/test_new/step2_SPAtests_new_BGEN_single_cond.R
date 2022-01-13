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

SPAGMMATtest(bgenFile="../input/genotype_100markers.bgen",
             bgenFileIndex = "../input/genotype_100markers.bgen.bgi",
             SAIGEOutputFile="./new_single_bgen_cond.txt",
             sampleFile="../input/samplelist.txt",
	     chrom="1",
	     GMMATmodelFile="../output/example_binary.rda",
	     varianceRatioFile="../output/example_binary.varianceRatio.txt.fake",
	     LOCO=FALSE,
	     AlleleOrder="ref-first",
	     is_output_moreDetails=TRUE,
	     max_markers_region = 10,
	     condition=c("rs18,rs33,rs50")
	     )

