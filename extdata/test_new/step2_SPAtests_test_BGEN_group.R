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
             SAIGEOutputFile="./new_group_bgen.txt",
             sampleFile="../input/samplelist.txt",
	     chrom="1",
	     groupFile="./group_plink.txt",
	     maxMAFforGroupTest=c(0.01,0.001,0.05,0.1),
	     function_group_test =c("lof", "missense", "lof;missense"),
	     GMMATmodelFile="../output/example_binary.rda",
	     varianceRatioFile="../output/example_binary_cate_v2.varianceRatio.txt",
	     LOCO=FALSE,
	     AlleleOrder="ref-first",
	     is_output_moreDetails=TRUE,
	     max_markers_region = 10,
	     sparseSigmaFile="../output/example_binary_cate_v2.varianceRatio.txt_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseSigma.mtx"
	     )

