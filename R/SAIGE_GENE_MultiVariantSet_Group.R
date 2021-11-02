# read group file
# should be 
# GeneID GroupName Markers1,...
ReadGroupFile<-function(groupFile){
  
  group_info_list<-list()
  line=0
  gf = file(groupFile, "r")
  while (TRUE) {
    line = line+1
    marker_group_line = readLines(gf, n = 1)
    if (length(marker_group_line) == 0) {
      break
    }
    marker_group_line_list = strsplit(marker_group_line, split=c(" +", "\t"))[[1]]
    if (length(marker_group_line_list) < 3) {
      stop("Error, group file line:",line ,"Each line should have geneID GroupName  \n")
    }
    
    geneID = marker_group_line_list[1]
    GroupName = marker_group_line_list[2]
    
    geneinfo = list(geneID = geneID, GroupName=GroupName, markerID= marker_group_line_list[-c(1:2)])
    if(is.null(group_info_list[[geneID]])){
      group_info_list[[geneID]]<-list()
    }
    group_info_list[[geneID]][[GroupName]]<-geneinfo
                              
  }
  close(gf)
  
  # Find set with all groups
  ngroup<-length(group_info_list)
  for(i in 1:ngroup){
    n1<-length(group_info_list[[i]])
    all_marker<-NULL
    geneID=group_info_list[[i]][[1]]$geneID
    if(!is.null(group_info_list[[i]][["all"]])){
      stop("Error: all should not be used as a GroupName!")
    }
    
    for(j in 1:n1){
      all_marker<-union(all_marker, group_info_list[[i]][[j]]$markerID)
    }
    group_info_list[[geneID]][["all"]]<-list(geneID = geneID, GroupName="all", markerID= all_marker)
  }
  
  return(group_info_list)
}

  