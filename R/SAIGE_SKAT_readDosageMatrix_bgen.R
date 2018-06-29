
getGenoOfGene_bgen = function(bgenFile,bgenFileIndex,marker_group_line, minMAF=0, maxMAF=0.5){
  ids_to_exclude = as.character(vector())
  ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
  ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
  ids_to_include = strsplit(marker_group_line, split=",")[[1]]   
  Mtest = setgenoTest_bgenDosage(bgenFile,bgenFileIndex, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
  Gvec = NULL
  markerIDs = NULL
  markerAFs = NULL
  cnt = 0
  result = list()
  if(Mtest > 0){    
    for(i in 1:Mtest){
      Gx = getDosage_bgen_withquery()
      AF = Gx$variants$AF
      if(AF >= 0.5){
        MAF = 1 - AF
      }else{
        MAF = AF
      }
      if(MAF >= minMAF && MAF < maxMAF){
        Gvec = c(Gvec, Gx$dosages)
        markerIDs = c(markerIDs, Gx$variants$SNPID)
        markerAFs = c(markerIDs, Gx$variants$AF)
        cnt = cnt + 1
      }
    }
    result$dosages = Gvec
    result$markerIDs = markerIDs
    result$markerAFs = markerAFs
    result$cnt = cnt
  }else{
    result$cnt = 0
  }

  closetestGenoFile_bgenDosage()
  return(result)
}
