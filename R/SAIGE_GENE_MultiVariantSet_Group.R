# read group file
# should be 
# GeneID type Markers1,...
#	type: var - makerID
#	type: anno - annotation
#	type: weight - weight for each variant, weight does not need to be in the file
ReadGroupFile<-function(groupFile){
  
  #groupFile = "./test_multiset/group_multiSets.txt";
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
    type = marker_group_line_list[2]
    
    # check type 
    if(!(type %in% c("var", "anno"))){
    	stop("Error, group file, 2nd column should be either var or anno, line: ", line, " \n")	 
    }
    
    if(is.null(group_info_list[[geneID]])){
      group_info_list[[geneID]]<-list(geneID = geneID)
    } 
    group_info_list[[geneID]][[type]] = marker_group_line_list[-c(1:2)]      
                  
  }
  close(gf)
  
  # check 
  ngroup<-length(group_info_list)
  for(i in 1:ngroup){
  	n1<-length(group_info_list[[i]])
  	if(n1 < 3){
  		stop("Error: group file, each gene should have var and anno:", i , " \n")	
  	}
  	
  	var_g = group_info_list[[i]][["var"]]
  	anno_g = group_info_list[[i]][["anno"]]
  	n_var = length(var_g)
    n_anno = length(anno_g)
	
  	if(n_var != n_anno){
  		stop("Error: group file, the length of var and anno are different:", i, n1, n2 , " \n")	
  	}
  	anno_name = unique(anno_g)
  	for(j in 1:length(anno_name)){
  		anno_name_1 = anno_name[j]
  		markerID = var_g[anno_g == anno_name_1]
  		group_info_list[[i]][[anno_name_1]]<-list(markerID= markerID)
  	}
  }  
  
  return(group_info_list)
}

  