
getGenoOfGene_bgen = function(bgenFile,bgenFileIndex,marker_group_line, minMAF=0, maxMAF=0.5, minInfo=0){
  ids_to_exclude = as.character(vector())
  ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
  ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
  idslist = strsplit(marker_group_line, split="\t")[[1]]
  ids_to_include = idslist[-1]	
  print("ids_to_include")
  print(ids_to_include)
  Mtest = setgenoTest_bgenDosage_v2(bgenFile,bgenFileIndex, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
  print(gc())
  rm(marker_group_line)
  rm(ids_to_include)
  rm(idslist)
  print(gc())
  cat("Mtest: ", Mtest, "\n")
  Gvec = NULL
  markerIDs = NULL
  markerAFs = NULL
  indexforMissing = NULL
  cnt = 0
  result = list()
  MACs = NULL
  print(gc())
  if(Mtest > 0){    
    for(i in 1:Mtest){
      Gx = getDosage_bgen_withquery()
      AF = Gx$variants$AF
      AC = Gx$variants$AC
      markerInfo = getMarkerInfo()
      if(AF >= 0.5){
        MAF = 1 - AF
        MAC = 1 - AC
      }else{
        MAF = AF
        MAC = AC
      }
    if(MAF >= minMAF & MAF < maxMAF & markerInfo >= minInfo){
        Gvec = c(Gvec, Gx$dosages)
	#print(object.size(Gvec))
        markerIDs = c(markerIDs, Gx$variants$rsid)
        markerAFs = c(markerAFs, AF)
        MACs = c(MACs, MAC)
        cnt = cnt + 1
	indexforMissing = c(indexforMissing, Gx$indexforMissing)
      }
    }
    indexforMissing = unique(indexforMissing)	
    result$dosages = Gvec
    result$markerIDs = markerIDs
    result$markerAFs = markerAFs
    result$MACs = MACs
    result$cnt = cnt
    result$indexforMissing = indexforMissing
    print("indexforMissing test")	
    print(indexforMissing)
  }else{
    result$cnt = 0
  }
  print(gc())
  closetestGenoFile_bgenDosage()
  return(result)
}


getGenoOfGene_bgen_Sparse = function(bgenFile,bgenFileIndex,marker_group_line, minMAF=0, maxMAF=0.5, minInfo=0){
  ids_to_exclude = as.character(vector())
  ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
  ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
  idslist = strsplit(marker_group_line, split="\t")[[1]]
  ids_to_include = idslist[-1]	
  print("ids_to_include")
  print(ids_to_include)
  Mtest = setgenoTest_bgenDosage_v2(bgenFile,bgenFileIndex, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
  print(gc())
  rm(marker_group_line)
  rm(ids_to_include)
  rm(idslist)
  print(gc())
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
      Gx = getDosage_bgen_withquery_Sparse()
      AF = Gx$variants$AF
      AC = Gx$variants$AC
      markerInfo = getMarkerInfo()
      if(AF >= 0.5){
        MAF = 1 - AF
        MAC = 1 - AC
      }else{
        MAF = AF
        MAC = AC
      }
    if(MAF >= minMAF & MAF < maxMAF & markerInfo >= minInfo){
        Gvec = c(Gvec, Gx$dosages)
	iIndex = c(iIndex, Gx$iIndexforMarker)
	jIndex = c(jIndex, rep((cnt+1), length(Gx$dosages)))
	print(i)
	print("object.size(Gvec)")
	print(object.size(Gvec))
	print("object.size(iIndex)")
	print(object.size(iIndex))
	print("object.size(jIndex)")
	print(object.size(jIndex))
	

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
    result$iIndex = iIndex
    result$jIndex = jIndex
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
