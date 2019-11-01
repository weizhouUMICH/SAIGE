
getGenoOfGene_bgen = function(bgenFile,bgenFileIndex,marker_group_line, minMAF=0, maxMAF=0.5, minInfo=0){
  ids_to_exclude = as.character(vector())
  ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
  ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
  idslist = strsplit(marker_group_line, split="\t")[[1]]
  ids_to_include = idslist[-1]	
  print("ids_to_include")
  print(ids_to_include)
  Mtest = setgenoTest_bgenDosage_v2(bgenFile,bgenFileIndex, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
  cat("Mtest: ", Mtest, "\n")
  Gvec = NULL
  markerIDs = NULL
  markerAFs = NULL
  cnt = 0
  result = list()
  MACs = NULL
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
        markerIDs = c(markerIDs, Gx$variants$rsid)
        markerAFs = c(markerAFs, MAF)
        MACs = c(MACs, MAC)
        cnt = cnt + 1
      }
    }
    result$dosages = Gvec
    result$markerIDs = markerIDs
    result$markerAFs = markerAFs
    result$MACs = MACs
    result$cnt = cnt
  }else{
    result$cnt = 0
  }

  closetestGenoFile_bgenDosage()
  return(result)
}
