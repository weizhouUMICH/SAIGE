# read group file
# should be 
# GeneID type Markers1,...
#	type: var - makerID
#	type: anno: - annotation
#	type: weight: - weight for each variant, weight does not need to be in the file
ReadGroupFile_New<-function(groupFile){
  
  #groupFile = "~/SAIGE_GENE_200WES_Pipeline/SAIGE_GENE_Test/group_multi.txt"
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
    
    # marker_group_line_list<-c("DEFB125", "anno:lof", "0", "0", "0")

    geneID = marker_group_line_list[1]
    type = marker_group_line_list[2]
    type.all = strsplit(type, ":")[[1]]
    type1 = type.all[1]
    
    # check type 
    if(!(type.all[1] %in% c("var", "anno", "weight"))){
    	stop("Error, group file, 2nd column should be either var/anno:/weight:, line: ", line, " \n")	 
    }
    
    if(is.null(group_info_list[[geneID]])){
      group_info_list[[geneID]]<-list(geneID = geneID)
    } 
    if(type1=="var"){
    	group_info_list[[geneID]][[type]] = marker_group_line_list[-c(1:2)]      
    } else {
    	type2 = type.all[2]
    	group_info_list[[geneID]][[type2]]= list(type=type1, val=as.numeric(marker_group_line_list[-c(1:2)])) 
    }          
  }
  close(gf)  
  
  return(group_info_list)
}



#############################
#  
ReadGroupFile<-function(groupFile){
  
  #groupFile = "./SAIGE_GENE_Test/group_multiSets_2.txt"
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

  
  
#
# info is the function_group_marker_list[[function_group]]
Extract_Weight<-function(markerIDs_extract, markerIDs, function_group_marker_list, function_group){
	# markerIDs_extract = info$marker_collapse_list
	if(length(markerIDs_extract)==0){
		return(NULL)	
	}
	idx_markerIDs_weight = match(markerIDs_extract, function_group_marker_list$var)
	markerIDs_weight = function_group_marker_list[[function_group]]$val[idx_markerIDs_weight]
	idx<-match(markerIDs_extract, markerIDs)
	weight = markerIDs_weight[idx]
	return(weight)
}

# Function return list of markers to collapse and not to collapse
# input 
#   markerIDs: IDs of all markers in Gmat
#   function_group_marker_list: marker IDs in each functional group. Markers here should not be in markerIDs
#   MACvec: MAC vector
#   MAF: MAF vector
#   function_group_test: function gruops to test. Should be ordered by functional importance
#   MAF_cutoff: MAF cutoff
#   MACCutoff_to_CollapseUltraRare: MACCutoff for collapsing
#
# output
#   FuncMAF_list: IDs of non-collapsed markers, by func group and MAF cutoff
#   marker_collapse_list: makers collapsed for each group. high priority group variants will be included in low priority groups
Get_MultiSet_Id_New<-function(markerIDs, function_group_marker_list, MACvec, MAF, MAF_cutoff=c(0.0001, 0.001, 0.01), MACCutoff_to_CollapseUltraRare = 10){
	

  #function_group_marker_list = group_info_list[[1]]; markerIDs = group_info_list[[1]]$var; MAF_cutoff=c(0.001, 0.01); MACCutoff_to_CollapseUltraRare = 10;function_group_test=c("lof", "lofmis", "w2")
  #MACvec = rep(20, length(markerIDs)); MACvec[1:10]<-1; MAF = MACvec / 10000
	n_cutoff<-length(MAF_cutoff)
		
	marker_collapse_all<-markerIDs[MACvec <= MACCutoff_to_CollapseUltraRare]
	FuncMAF_list<-list()
	MAF_group_list<-list()
	for(i in 1:n_cutoff){
		MAF_group_list[[i]]<-markerIDs[MAF <=MAF_cutoff[[i]]]
	}
	
	marker_collapse_list<-list()
	idx_all<-NULL
	group_without_collapse<-list()

	for(i in 1:length(function_group_test)){
	  
		function_group<-function_group_test[i]
		info<-function_group_marker_list[[function_group]]
		if(is.null(info)){
			marker = NULL
		} else {
			idx = which(info$val > 0)
			marker = function_group_marker_list[["var"]][idx]
		}
		marker_collapse_list<-intersect(marker, marker_collapse_all)
		marker_without_collapse<-setdiff(marker, marker_collapse_list)
		FuncMAF_list<-list()
		for(j in 1:n_cutoff){
			FuncMAF_list[[j]]<-intersect(marker_without_collapse, MAF_group_list[[j]])
		}
		function_group_marker_list[[function_group]]$marker_collapse_list<-marker_collapse_list
		function_group_marker_list[[function_group]]$marker_without_collapse<-marker_without_collapse
		function_group_marker_list[[function_group]]$FuncMAF_list <- FuncMAF_list
		
	}
	

	return(function_group_marker_list)
}


Get_MultiSet_Id<-function(markerIDs, function_group_marker_list, MACvec, MAF, 
                           function_group_test=c("lof", "missense"), MAF_cutoff=c(0.0001, 0.001, 0.01), MACCutoff_to_CollapseUltraRare = 10){
	

  #function_group_marker_list = group_info_list[[1]]; markerIDs = Gx$markerIDs; MAF_cutoff=c(0.001, 0.01); MACCutoff_to_CollapseUltraRare = 10;function_group_test=c("lof", "missense")
  
	n_cutoff<-length(MAF_cutoff)
		
	marker_collapse_all<-markerIDs[MACvec <= MACCutoff_to_CollapseUltraRare]
	FuncMAF_list<-list()
	MAF_group_list<-list()
	for(i in 1:n_cutoff){
		MAF_group_list[[i]]<-markerIDs[MAF <=MAF_cutoff[[i]]]
	}
	

	marker_collapse_list<-list()
	idx_all<-NULL
	group_without_collapse<-list()
	
	marker_collapse_prev_group<-NULL
	for(i in 1:length(function_group_test)){
	  
		function_group<-function_group_test[[i]]
		info<-function_group_marker_list[[function_group]]
		marker<-intersect(info$markerID, markerIDs)
		marker_collapse_list[[i]]<-union(marker_collapse_prev_group, intersect(marker, marker_collapse_all))
		group_without_collapse[[i]]<-setdiff(marker, marker_collapse_list[[i]])
		marker_collapse_prev_group<-marker_collapse_list[[i]]
	}
	
	# Make union without collapsing
	group<-NULL
	for(i in 1:length(function_group_test)){
		group<-c(group, group_without_collapse[[i]])
		FuncMAF_list[[i]]<-list()
		for(j in 1:n_cutoff){
			FuncMAF_list[[i]][[j]]<-intersect(group, MAF_group_list[[j]])
		}
	}

	re_group_id<-list(FuncMAF_list=FuncMAF_list, marker_collapse_list=marker_collapse_list)
	return(re_group_id)
}

  


#	input
#		m: number of markers
# 		re_group_idx: output from Get_MultiSet_Idx
# 		fun_group_name: function group name
# 	output
#		Gmat: new genotype matrix with collapsed variants
#		markerIDs_new: new marker ID
#		ncollapse: number of collapsed variants
#		collapse_indicator: list of indicators for each functional groups
Get_Collapsed_Genotype_New<-function(Gmat, markerIDs, m, function_group_marker_list, function_group_test, DosageCutoff_for_UltraRarePresence){

	
	marker_collapse_all_idx = NULL
	GCollapsing = NULL
	Collapsing_ID = NULL
	ncollapse=0
	collapse_indicator=list()
	
	collapse_indicator = rep(-1,length(function_group_test))

  	for(i in 1:length(function_group_test)){

		function_group<-function_group_test[i]
		info<-function_group_marker_list[[function_group]]
		macle10Index = match(info$marker_collapse_list, markerIDs)
		
		if(length(macle10Index) > 0 && info$type =="anno"){
			G1rare = Gmat[, macle10Index, drop = F]
			Gnew = qlcMatrix::rowMax(G1rare)
			ID1<-which(as.vector(Gnew < (1 + DosageCutoff_for_UltraRarePresence)))
			ID2<-which(as.vector(Gnew > DosageCutoff_for_UltraRarePresence))
      		ID3<-which(as.vector(Gnew >= (1 + DosageCutoff_for_UltraRarePresence)))
      
      		Gnew[intersect(ID1, ID2)] = 1
      		Gnew[ID3] = 2
      		Gnew = as(Gnew, "sparseMatrix")

    	} else {
    		weights = Extract_Weight(info$marker_collapse_list, markerIDs, function_group_marker_list, function_group)
    		G1rare = Gmat[, macle10Index, drop = F]
    		Gnew = Matrix::colSums(t(G1rare) * (weights))	
 			Gnew = as(Gnew, "sparseMatrix")
    	
    	}
      	
      	GCollapsing<-cbind(GCollapsing, Gnew)
        	
      	Collapsing_ID = c(Collapsing_ID, sprintf("C_%s", function_group_test[i] ))
      	ncollapse = ncollapse+1
      	collapse_indicator[i] = ncollapse
		marker_collapse_all_idx = union(marker_collapse_all_idx, macle10Index )
	}

	
	if(length(marker_collapse_all_idx)==0){
		# no markers for collapsing
		Gmat=Gmat
	  	markerIDs_new=markerIDs
	} else if (length(marker_collapse_all_idx) < m) {
		Gmat = cbind(GCollapsing, Gmat[, -marker_collapse_all_idx, drop = F])
	  	markerIDs_nocol = markerIDs[-marker_collapse_all_idx]
	  	markerIDs_new = c(Collapsing_ID, markerIDs[-marker_collapse_all_idx])
  	} else {
    	Gmat = cbind(Gnew)
    	markerIDs_new = Collapsing_ID
  	}
	
	return(list(Gmat=Gmat, markerIDs_new=markerIDs_new, ncollapse= ncollapse, collapse_indicator=collapse_indicator))

}

    

#	input
#		m: number of markers
# 		re_group_idx: output from Get_MultiSet_Idx
# 		fun_group_name: function group name
# 	output
#		Gmat: new genotype matrix with collapsed variants
#		markerIDs_new: new marker ID
#		ncollapse: number of collapsed variants
#		collapse_indicator: list of indicators for each functional groups
Get_Collapsed_Genotype<-function(Gmat, markerIDs, m, re_group_id, function_group_test, DosageCutoff_for_UltraRarePresence){

	
	marker_collapse_list = re_group_id$marker_collapse_list
	marker_collapse_all_idx = NULL
	GCollapsing = NULL
	Collapsing_ID = NULL
	ncollapse=0
	collapse_indicator=list()
	
	collapse_indicator = rep(-1,length(marker_collapse_list))

  for(i in 1:length(marker_collapse_list)){
	
		macle10Index = match(marker_collapse_list[[i]], markerIDs)
		
		if(length(macle10Index) > 0){
			G1rare = Gmat[, macle10Index, drop = F]
			Gnew = qlcMatrix::rowMax(G1rare)
			ID1<-which(as.vector(Gnew < (1 + DosageCutoff_for_UltraRarePresence)))
			ID2<-which(as.vector(Gnew > DosageCutoff_for_UltraRarePresence))
      ID3<-which(as.vector(Gnew >= (1 + DosageCutoff_for_UltraRarePresence)))
      
      Gnew[intersect(ID1, ID2)] = 1
      Gnew[ID3] = 2
      Gnew = as(Gnew, "sparseMatrix")
      GCollapsing<-cbind(GCollapsing, Gnew)
        	
      Collapsing_ID = c(Collapsing_ID, sprintf("C_%s", function_group_test[i] ))
      ncollapse = ncollapse+1
      collapse_indicator[i] = ncollapse
    }
		marker_collapse_all_idx = union(marker_collapse_all_idx, macle10Index )
	}

	if(length(marker_collapse_all_idx)==0){
	  Gmat=Gmat
	  markerIDs_new=markerIDs
	} else if (length(marker_collapse_all_idx) < m) {
	  Gmat = cbind(GCollapsing, Gmat[, -marker_collapse_all_idx, drop = F])
	  markerIDs_nocol = markerIDs[-marker_collapse_all_idx]
	  markerIDs_new = c(Collapsing_ID, markerIDs[-marker_collapse_all_idx])
  } else {
    Gmat = cbind(Gnew)
    markerIDs_new = Collapsing_ID
  }
	
	return(list(Gmat=Gmat, markerIDs_new=markerIDs_new, ncollapse= ncollapse, collapse_indicator=collapse_indicator))

}

    
    