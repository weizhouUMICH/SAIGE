options(stringsAsFactors=F)

## load R libraries
#library(SAIGE)
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.35.2.mmSKAT.debugged.R-3.5.1.test2_subsetSparseSigma_speedup_test2")
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.35.3.2")
library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/0.35.3.3")


require(optparse) #install.packages("optparse")

## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile", type="character",default="",
    help="path to plink file for creating the genetic relationship matrix (GRM)"),
    make_option("--nThreads", type="integer", default=16,
    help="Number of threads (CPUs) to use"),
  make_option("--memoryChunk", type="numeric", default=2,
   help="The size (Gb) for each memory chunk. By default, 2"),
  make_option("--outputPrefix", type="character", default="~/",
    help="path and prefix to the output files [default='~/']"),
  make_option("--numRandomMarkerforSparseKin", type="integer", default=500,
    help="number of randomly selected markers (MAF >= 1%) to be used to identify related samples for sparse GRM [default=500]"),
  make_option("--relatednessCutoff", type="numeric", default=0.125,
    help="The threshold to treat two samples as unrelated if IsSparseKin is TRUE [default=0.125]"),
  make_option("--isDiagofKinSetAsOne", type="logical", default=FALSE,
    help="Whether to set the diagnal elements in GRM to be 1 [default='FALSE'].")
)



## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

set.seed(1)

createSparseGRM(plinkFile = opt$plinkFile,
                outputPrefix=opt$outputPrefix,
                numRandomMarkerforSparseKin = opt$numRandomMarkerforSparseKin,
                relatednessCutoff = opt$relatednessCutoff,
                memoryChunk = opt$memoryChunk,
                isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
                nThreads = opt$nThreads)
