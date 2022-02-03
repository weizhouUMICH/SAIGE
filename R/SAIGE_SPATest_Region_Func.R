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
  #print("group_info_list")
  #print(group_info_list)
  return(group_info_list)
}


check_close= function(groupFile){
 if(isOpen(groupFile)){
    close(groupFile)
  }	 
}


#check Group file format and count the number of regions
checkGroupFile<-function(groupFile){

  cat("Start extracting marker-level information from 'groupFile' of", groupFile, "....\n")	
  Check_File_Exist(groupFile, "RegionFile")	
  gf = file(groupFile, "r")
  marker_group_line = readLines(gf, n = 1)
  marker_group_line = readLines(gf, n = 1)
  is_weight_included = FALSE
  a = 2
  marker_group_line = readLines(gf, n = 1)
  if(length(marker_group_line) == 1){
    if (length(marker_group_line) == 0) {
      line = line - 1
      stop("Error, group file has emply lines\n")
    }
    marker_group_line_list = strsplit(marker_group_line, split=c(" +", "\t"))[[1]]
    if (length(marker_group_line_list) < 3) {
      check_close(gf)
      stop("Error, group file line:",line ,"Each line should have a region name and a catergory (var, anno or weight) \n")
    }
    geneID = marker_group_line_list[1]
    type = marker_group_line_list[2]
    if(type == "weight"){
        is_weight_included = TRUE
        print("weights are included for markers")
        a = 3
    }
  }

  close(gf)

  gf = file(groupFile, "r")
  line=0
  nregion = 0
  while (TRUE) {

    marker_group_line = readLines(gf, n = a)
    #line = line+a

    if (length(marker_group_line) == 0 ){	    
    	break
    }else{
	if(length(marker_group_line) < a){
		marker_group_line_list = strsplit(marker_group_line[1], split=c(" +", "\t"))[[1]]	
		if (length(marker_group_line_list) < 3) {
      			stop("Error, group file line:",line-2 ," is incomplete.\n")
    		}
		geneID = marker_group_line_list[1]
		stop("Group file is incomplete for ", geneID,".\n")
	}	
    }	    


    marker_group_line_list = strsplit(marker_group_line[1], split=c(" +", "\t"))[[1]]
    line = line + 1
    if (length(marker_group_line_list) < 3) {
      	stop("Error, group file line:",line ," is incomplete.\n")
    }    
    geneID = marker_group_line_list[1]
    type = marker_group_line_list[2]

    if(type != "var"){
	stop("Error! No var is specified for ", geneID, "\n")
    }
    geneID0 = geneID
    numMarkers = length(marker_group_line_list) - 2	

    marker_group_line_list = strsplit(marker_group_line[2], split=c(" +", "\t"))[[1]]
    line = line + 1
   if (length(marker_group_line_list) < 3) {
        stop("Error, group file line:",line ," is incomplete.\n")
    }
    geneID = marker_group_line_list[1]
    type = marker_group_line_list[2]
	
    if(type != "anno"){
	stop("Error! No anno is specified for ", geneID, "\n")
    }
    if(geneID != geneID0){
	stop("anno for ", geneID0, " is missing.\n")
    }	    
    numAnnos = length(marker_group_line_list) - 2
    if(numAnnos != numMarkers){
	stop("The length of annotations for markers in region ", geneID, " is not equal to the length of marker IDs\n")
    }	

    if(is_weight_included){

   	marker_group_line_list = strsplit(marker_group_line[3], split=c(" +", "\t"))[[1]]
        line = line + 1
	if (length(marker_group_line_list) < 3) {
        	stop("Error, group file line:",line ," is incomplete.\n")
    	}
	geneID = marker_group_line_list[1]
    	type = marker_group_line_list[2]
	if(type != "weight"){
        	stop("Error! No weight is specified for ", geneID, "\n")
    	}
    	if(geneID != geneID0){
        	stop("weight for ", geneID0, " is missing.\n")
    	}
	numWeights = length(marker_group_line_list) - 2
        if(numWeights != numMarkers){
                stop("The length of weights for markers in region ", geneID, " is not equal to the length of marker IDs\n")
        }
    }
   nregion = nregion + 1
   }    
    #if(is.null(group_info_list[[geneID]])){
    #  group_info_list[[geneID]]<-list(geneID = geneID)
    #}
    #group_info_list[[geneID]][[type]] = marker_group_line_list[-c(1:2)]

  close(gf)	
  return(list(nRegions = nregion, is_weight_included = is_weight_included))
}	

SPA_ER_kernel_related_Phiadj_fast_new<-function(p.new, Score, Phi, p.value_burden, regionTestType){ 
	p.m = length(Score)
        zscore.all_0 = Score
	zscore.all_1<-rep(0, p.m)
	if(regionTestType != "BURDEN"){
        	VarS_org=diag(as.matrix(Phi))
	}else{
		VarS_org = as.vector(Phi)
	}	
        stat.qtemp =Score^2/VarS_org
	

	idx_0<-which(VarS_org >0)
        idx_p0<-which(p.new >0)
        idx_p1<-which(p.new <0)


	#if(length(idx_0) > 0){
 	#	zscore.all_1[idx_0]= qnorm(p.new[idx_0]/2, mean = 0, sd =sqrt(VarS_org[idx_0]),lower.tail = FALSE, log.p = FALSE)*sign(Score)
        #}
	VarS = zscore.all_0^2/500

	if(length(idx_p0) > 0){
                VarS[idx_p0]= zscore.all_0[idx_p0]^2/qchisq(p.new[idx_p0], 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
        }
	vars_inf=which(VarS==Inf)
	if(regionTestType != "BURDEN"){
        	if (length(vars_inf)>0){
                	VarS[vars_inf] = 0
                	zscore.all_1[vars_inf]=0
                	Phi[vars_inf,]=0
                	Phi[,vars_inf]=0
        	}
	}else{
		if (length(vars_inf)>0){
			VarS[vars_inf] = 0
                	zscore.all_1[vars_inf]=0
                	Phi[vars_inf]=0
		}
	}	


	scaleFactor = sqrt(VarS/VarS_org)
	###################################
        # Burden test
	#print("SPA_ER_kernel_related_Phiadj_fast_new 5")
        zscore.all_1 = zscore.all_0
	#print(as.matrix(Phi))
	#print(diag(as.vector(VarS/VarS_org)))
	VarStoorg = as.vector(VarS/VarS_org)

	if(regionTestType != "BURDEN"){
	#t(sqrt(a)*t(b)*sqrt(a))
		G2_adj_n=t(t(Phi * sqrt(VarStoorg)) * sqrt(VarStoorg))
        	#v1=rep(1,nrow(G2_adj_n))
	}else{
		G2_adj_n = Phi * VarStoorg	
	}
        	VarQ = sum(G2_adj_n)
        	Q_b=sum(zscore.all_1)^2
        	VarQ_2=Q_b/qchisq(p.value_burden, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
        	if (VarQ_2== 0) {
                	r=1
        	}else{
                	r=VarQ/VarQ_2
        	}
        	r=min(r,1)

	if(regionTestType != "BURDEN"){

        	Phi_ccadj=as.matrix(G2_adj_n * 1/r)
	}else{
		Phi_ccadj=as.vector(G2_adj_n * 1/r)
	}	

        outlist=list();
        outlist$val = Phi_ccadj

        scaleFactor = scaleFactor /sqrt(r)
        outlist$scaleFactor = scaleFactor
        outlist$p.new = p.new
        return(outlist)
}


get_newPhi_scaleFactor = function(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType){
	p.value_burden<-SPAtest:::Saddle_Prob(q.sum , mu=mu.a, g=g.sum, Cutoff=2,alpha=2.5*10^-6)$p.value
        re_phi= SPA_ER_kernel_related_Phiadj_fast_new(p.new, Score, Phi, p.value_burden, regionTestType)
	#print("here2")
	return(re_phi)
}

get_SKAT_pvalue = function(Score, Phi, r.corr, regionTestType){
	if(regionTestType == "BURDEN"){
		Q = try((SKAT:::SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r), silent = TRUE)
		Q.res = NULL
        	a <- as.matrix(sum(Phi))
        	out_SKAT_List <- try(SKAT:::Get_Liu_PVal(Q, a, Q.res), silent = TRUE)
	}else{

                out_SKAT_List = try(SKAT:::Met_SKAT_Get_Pvalue(Score = Score,
                                                         Phi = Phi,
                                                         r.corr = r.corr,
                                                         method = "optimal.adj",
                                                         Score.Resampling = NULL),
                              silent = TRUE)
	}


	if(regionTestType != "BURDEN"){
                BETA_Burden = sum(Score)/(sum(diag(Phi)))
                if(class(out_SKAT_List) == "try-error"){
                        Pvalue = c(NA, NA, NA)
                        error.code = 2
                        BETA_Burden = NA
                        SE_Burden = NA
                }else if(!any(c(0,1) %in% out_SKAT_List$param$rho & !is.null(out_SKAT_List$p.value))){
                        #Pvalue = c(NA, NA, NA)
                        Pvalue = c(out_SKAT_List$p.value, out_SKAT_List$p.value, out_SKAT_List$p.value)
                        #error.code = 3
                        error.code = 0
                        #BETA_Burden = NA
                        #SE_Burden = NA
			SE_Burden = abs(BETA_Burden/qnorm((out_SKAT_List$p.value)/2))
                }else{
                        pos00 = which(out_SKAT_List$param$rho == 0)
                        pos01 = which(out_SKAT_List$param$rho == 1)
                        Pvalue = c(out_SKAT_List$p.value,                # SKAT-O
                                out_SKAT_List$param$p.val.each[pos00],   # SKAT
                                out_SKAT_List$param$p.val.each[pos01])   # Burden Test
                        error.code = 0
                        SE_Burden = abs(BETA_Burden/qnorm((out_SKAT_List$param$p.val.each[pos01])/2))

                }
		return(list(Pvalue_SKATO = Pvalue[1], Pvalue_Burden = Pvalue[3], Pvalue_SKAT = Pvalue[2], BETA_Burden = BETA_Burden, SE_Burden = SE_Burden))
	}else{
		BETA_Burden = sum(Score)/(sum(Phi))
		Pvalue_Burden = out_SKAT_List$p.value
		SE_Burden = abs(BETA_Burden/qnorm(Pvalue_Burden/2))
		Pvalue_SKATO = NA
		Pvalue_SKAT = NA
		return(list(Pvalue_SKATO = NA, Pvalue_Burden = Pvalue_Burden, Pvalue_SKAT = NA, BETA_Burden = BETA_Burden, SE_Burden = SE_Burden))
	}	

}

get_CCT_pvalue = function(pvalue){
   pvals = pvalue
   notna = which(!is.na(pvals))
   if(length(notna) > 0){
     pvals = pvals[!is.na(pvals)]
     cctpval = CCT(pvals)
   }else{
     cctpval = NA
   }
   return(cctpval)
}	
