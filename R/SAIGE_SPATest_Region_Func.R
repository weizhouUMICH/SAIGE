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
  print("group_info_list")
  print(group_info_list)
  return(group_info_list)
}




SPA_ER_kernel_related_Phiadj_fast_new<-function(p.new, Score, Phi, p.value_burden){ 
	#print("SPA_ER_kernel_related_Phiadj_fast_new 0")
	#cat("p.value_burden ", p.value_burden, "\n")
	p.m = length(Score)
	#cat("Score ", Score, "\n")
	#cat("p.m ", p.m, "\n")
        zscore.all_0 = Score
	zscore.all_1<-rep(0, p.m)
	#print("zscore.all_1")
	#print(zscore.all_1)
        VarS_org=diag(as.matrix(Phi))
        stat.qtemp =Score^2/VarS_org
	#print("SPA_ER_kernel_related_Phiadj_fast_new 1")
	

	idx_0<-which(VarS_org >0)
        idx_p0<-which(p.new >0)
        idx_p1<-which(p.new <0)
	#print("SPA_ER_kernel_related_Phiadj_fast_new 2")


	if(length(idx_0) > 0){
 		zscore.all_1[idx_0]= qnorm(p.new[idx_0]/2, mean = 0, sd =sqrt(VarS_org[idx_0]),lower.tail = FALSE, log.p = FALSE)*sign(Score)
        }
	VarS = zscore.all_0^2/500
        if(length(idx_p0) > 0){
                VarS[idx_p0]= zscore.all_0[idx_p0]^2/qchisq(p.new[idx_p0], 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
        }
	#print("VarS")
	#print(VarS)
	#print("VarS_org")
	#print(VarS_org)
	#print("SPA_ER_kernel_related_Phiadj_fast_new 3")

	vars_inf=which(VarS==Inf)
        if (length(vars_inf)>0){
                VarS[vars_inf] = 0
                zscore.all_1[vars_inf]=0
                Phi[vars_inf,]=0
                Phi[,vars_inf]=0
        }

	#print("SPA_ER_kernel_related_Phiadj_fast_new 4")
	scaleFactor = sqrt(VarS/VarS_org)

	###################################
        # Burden test
	#print("SPA_ER_kernel_related_Phiadj_fast_new 5")

	#print(as.matrix(Phi))
	#print(diag(as.vector(VarS/VarS_org)))
	VarStoorg = as.vector(VarS/VarS_org)
	#t(sqrt(a)*t(b)*sqrt(a))
	G2_adj_n=t(t(Phi * sqrt(VarStoorg)) * sqrt(VarStoorg))
        v1=rep(1,nrow(G2_adj_n))
        VarQ = sum(G2_adj_n)
	#print("zscore.all_1")
	#print(zscore.all_1)
        Q_b=sum(zscore.all_1)^2
	#print("SPA_ER_kernel_related_Phiadj_fast_new 6")
	#cat("p.value_burden ", p.value_burden, "\n")
	#cat("Q_b ", Q_b, "\n")
        VarQ_2=Q_b/qchisq(p.value_burden, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
        if (VarQ_2== 0) {
                r=1
        } else {
                r=VarQ/VarQ_2
        }
        r=min(r,1)

        Phi_ccadj=as.matrix(G2_adj_n * 1/r)

        outlist=list();
        outlist$val = Phi_ccadj
        scaleFactor = scaleFactor /sqrt(r)

        outlist$scaleFactor = scaleFactor
        outlist$p.new = p.new
        return(outlist)
}


get_newPhi_scaleFactor = function(q.sum, mu.a, g.sum, p.new, Score, Phi){
	#print("here0")
	#q=q.sum
	#m1 <- sum(mu.a * g.sum)
        #var1 <- sum(mu.a * (1 - mu.a) * g.sum^2)
    	#cat("q ", q, "\n")
    	#cat("m1 ", m1, "\n")
    	#cat("var1 ", var1, "\n")
	p.value_burden<-SPAtest:::Saddle_Prob(q.sum , mu=mu.a, g=g.sum, Cutoff=2,alpha=2.5*10^-6)$p.value
	#print("here1")
	#print(p.new)
	#print(Score)
	#print(Phi)
	#print(p.value_burden)
        re_phi= SPA_ER_kernel_related_Phiadj_fast_new(p.new, Score, Phi, p.value_burden)
	#print("here2")
	return(re_phi)
}

get_SKAT_pvalue = function(Score, Phi, r.corr){
                out_SKAT_List = try(SKAT:::Met_SKAT_Get_Pvalue(Score = Score,
                                                         Phi = Phi,
                                                         r.corr = r.corr,
                                                         method = "optimal.adj",
                                                         Score.Resampling = NULL),
                              silent = TRUE)
                BETA_Burden = sum(Score)/(sum(diag(Phi)))
                if(class(out_SKAT_List) == "try-error"){
                        Pvalue = c(NA, NA, NA)
                        error.code = 2
                        BETA_Burden = NA
                        SE_Burden = NA
                }else if(!any(c(0,1) %in% out_SKAT_List$param$rho)){
                        Pvalue = c(NA, NA, NA)
                        error.code = 3
                        BETA_Burden = NA
                        SE_Burden = NA
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

}	
