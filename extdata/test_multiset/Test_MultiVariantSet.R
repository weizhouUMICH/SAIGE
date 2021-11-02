library(SAIGE)
library(SKAT)
setwd("~/Github/SAIGE/extdata")
source("~/Github/SAIGE/R/SAIGE_GENE_MultiVariantSet_Func.R")
source("~/Github/SAIGE/R/SAIGE_GENE_MultiVariantSet_Group.R")
source("~/Github/SAIGE/R/SAIGE_GENE_MultiVariantSet_Test.R")
#source("/Users/lee7801/Dropbox/R_Packages/SAIGE/R/SAIGE_GENE_MultiVariantSet.R")


source("~/Github/SAIGE/extdata/test_multiset/Test_MultiVariantSet_Code.R")
re_out


setwd("~/Github/SAIGE/extdata")
source("~/Github/SAIGE/R/SAIGE_GENE_MultiVariantSet.R")
out = SAIGE_GENE_MultiVariantSets(vcfFile = "./input/genotype_10markers.vcf.gz",
                            vcfFileIndex = "./input/genotype_10markers.vcf.gz.tbi",
                            sampleFile = "./input/samplelist.txt",
                            GMMATmodelFile = "./output/example_quantitative.rda",
                            varianceRatioFile = "./output/example_quantitative_cate.varianceRatio.txt",
                            SAIGEOutputFile = "./test_multiset/example_quantitative.SAIGE.gene.txt",
                            sparseSigmaFile="./output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx",
                            groupFile="./test_multiset/group_multiSets.txt",
                            LOCO=FALSE,
                            chrom = "1",
                            start = 1,
                            end = 250000000,
                            function_group_test =c("lof", "missense"),
                            MAF_cutoff=c(0.001,0.01)
                            )
out

Gx = getGenoOfGene_vcf(marker_group_line, minInfo)