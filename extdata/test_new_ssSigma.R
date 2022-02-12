> library(Matrix)
> ssigma_old = readMM("./output/example_binary_fullGRM_sparseGRM_categorical_varRatio.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx")
> library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_old_check/install_0.93")
> sparseGRMFile="output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
> sparseGRMSampleIDFile="output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
> load("example_binary_fullGRM.rda")
Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection
In addition: Warning message:
In readChar(con, 5L, useBytes = TRUE) :
  cannot open compressed file 'example_binary_fullGRM.rda', probable reason 'No such file or directory'
> library(SAIGE)
> load("example_binary_fullGRM.rda")
Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection
In addition: Warning message:
In readChar(con, 5L, useBytes = TRUE) :
  cannot open compressed file 'example_binary_fullGRM.rda', probable reason 'No such file or directory'
> load("output/example_binary_fullGRM.rda")
> obj.glmm.null = modglmm
> SAIGE:::setSparseSigma_new(sparseGRMFile, sparseGRMSampleIDFile, obj.glmm.null$sampleID, obj.glmm.null$theta, obj.glmm.null$mu2,  obj.glmm.null$traitType)
> new_ss_list = SAIGE:::setSparseSigma_new(sparseGRMFile, sparseGRMSampleIDFile, obj.glmm.null$sampleID, obj.glmm.null$theta, obj.glmm.null$mu2,  obj.glmm.null$traitType)
