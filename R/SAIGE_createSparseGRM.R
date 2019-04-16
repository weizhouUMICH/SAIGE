#' Construct a sparse GRM for a given data set 
#'
#' @param plinkFile character. Path to plink file to be used for calculating the sparse GRM
#' @param outputPrefix character. Path to the output files with prefix
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples for sparse GRM. By default, 1000
#' @param relatednessCutoff float. The threshold to treat two samples as unrelated if IsSparseKin is TRUE. By default, 0.125
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param isDiagofKinSetAsOne  logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param nThreads integer. Number of threads to be used. By default, 1 
#' @return a file ended with sampleIDs.txt that contains sample IDs for the sparse GRM and a file ended with .sparseGRM.mtx that contains the sparse GRM 
#' @export
createSparseGRM = function(plinkFile = "", 
		outputPrefix="",
                numRandomMarkerforSparseKin = 1000,
                relatednessCutoff = 0.125,
		memoryChunk = 2,
	        isDiagofKinSetAsOne = FALSE,
		nThreads = 1
                ){

  if(nThreads > 1){
    RcppParallel:::setThreadOptions(numThreads = nThreads)
    cat(nThreads, " threads are set to be used ", "\n")
  }

  cat("sparse GRM will be created\n")
  #  
  famFile = paste0(plinkFile, ".fam")
  fam = data.frame(data.table:::fread(famFile, header=F, stringsAsFactors=FALSE))
  sparseGRMSampleID = fam[,2]
  sparseGRMSampleIDFile = paste0(outputPrefix,"_relatednessCutoff_",relatednessCutoff,"_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")

  cat("write sample IDs for the sparse GRM to ", sparseGRMSampleIDFile ,"\n")
  write.table(sparseGRMSampleID, sparseGRMSampleIDFile, quote=F, col.names=F, row.names=F)


  genoSampleIndex = seq(1, nrow(fam))
  setgeno(plinkFile, genoSampleIndex, memoryChunk, isDiagofKinSetAsOne)

    freqVec = getAlleleFreqVec()
    MAFindex = which(freqVec >= 0.01 & freqVec <= 0.99)
    cat(numRandomMarkerforSparseKin, "genetic markers are randomly selected to decide which samples are related\n")
    if(length(MAFindex) < numRandomMarkerforSparseKin){
      stop("ERROR! not enough genetic markers with MAF >= 1% to detect which samples are related\n","Try include at least ", numRandomMarkerforSparseKin, " genetic markers with MAF >= 1% in the plink file\n")
    }

    markerIndexforSparseM = sample(MAFindex, size = numRandomMarkerforSparseKin, replace=FALSE)

    cat("Start detecting related samples for the sparse GRM\n")
    ta = proc.time()
    setSubMarkerIndex(markerIndexforSparseM -1)
    tb = proc.time()
    cat("tb-ta\n")
    print(tb-ta)


    cat("Start creating sparse GRM\n")
    ta = proc.time()
    sparseMList = createSparseKinParallel(nblocks = nThreads, ncore = nThreads, relatednessCutoff)
    tb = proc.time()
    cat("tb-ta\n")
    print(tb-ta)



    cat("length(sparseMList$iIndex): ", length(sparseMList$iIndex), "\n")
    print(sparseMList$iIndex[1:102])
    cat("length(sparseMList$jIndex): ", length(sparseMList$jIndex), "\n")
    print(sparseMList$jIndex[1:102])
    cat("length(sparseMList$kinValue): ", length(sparseMList$kinValue), "\n")
    print(sparseMList$kinValue[1:102])
    sparseGRM = Matrix:::sparseMatrix(i = as.vector(sparseMList$iIndex), j = as.vector(sparseMList$jIndex), x = as.vector(sparseMList$kinValue), symmetric = TRUE)
    cat("nrow(sparseGRM): ", nrow(sparseGRM), "\n")
    cat("ncol(sparseGRM): ", ncol(sparseGRM), "\n")
    cat("ncol(sparseGRM): ", sum(sparseGRM != 0), "\n")

    tc = proc.time()
    cat("tc-tb\n")
    print(tc-tb)


    sparseGRMFile = paste0(outputPrefix,"_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx")

  cat("write sparse GRM to ", sparseGRMFile ,"\n")
  Matrix:::writeMM(sparseGRM, sparseGRMFile)
  return(1)
}
