getGenoOfGene_bgen = function(markerIndicesVec, minMAF=0, maxMAF=0.5, minInfo=0, isSparseDosages = TRUE){
  setIsSparseDosage_bgen(isSparseDosages)
  setMarkerIndicesToInclude(markerIndicesVec)
  Mtest = length(markerIndicesVec)
  cat("Mtest: ", Mtest, "\n")
  Gvec = NULL
  iIndex = NULL
  jIndex = NULL
  markerIDs = NULL
  markerAFs = NULL
  indexforMissing = NULL
  cnt = 0
  result = list()
  MACs = NULL
  print(gc())
  if(Mtest > 0){    
    for(i in 1:Mtest){
      Gx = getOneMarker(markerIndicesVec[i])	    
      AF = Gx$variants$AF
      AC = Gx$variants$AC
      markerInfo = Gx$info
      if(AF >= 0.5){
        MAF = 1 - AF
        MAC = 1 - AC
      }else{
        MAF = AF
        MAC = AC
      }
    if(MAF >= minMAF & MAF < maxMAF & markerInfo >= minInfo){
        Gvec = c(Gvec, Gx$dosages)
        if(isSparseDosages){
	  iIndex = c(iIndex, Gx$iIndexforMarker)
	  jIndex = c(jIndex, rep((cnt+1), length(Gx$dosages)))
        }
	#print(i)
        markerIDs = c(markerIDs, Gx$variants$rsid)
        markerAFs = c(markerAFs, AF)
        MACs = c(MACs, MAC)
	indexforMissing = c(indexforMissing, Gx$indexforMissing)
        cnt = cnt + 1
      }
    }
    indexforMissing = unique(indexforMissing)
    result$dosages = Gvec
    result$indexforMissing = indexforMissing
    if(isSparseDosages){
      result$iIndex = iIndex
      result$jIndex = jIndex
    }
    result$markerIDs = markerIDs
    result$markerAFs = markerAFs
    result$MACs = MACs
    result$cnt = cnt
  }else{
    result$cnt = 0
  }
  #print(gc(T))
  closetestGenoFile_bgenDosage()
  return(result)
}
