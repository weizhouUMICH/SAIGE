

# Currently, conditional analysis is not supported
# Only absence_or_presence are supported
# idx_function_group has the most functional (ex. loss of function) to least functional


MultiSets_GroupTest = function (Gmat, obj.model, obj_cc, y, X, tauVec, traitType, 
                                function_group_marker_list, MAF_cutoff=c(0.0001, 0.001, 0.01), function_group_test=c("lof", "missense", "synonymous"),
                                cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, 
                                kernel, method, weights.beta.rare,  weightMAFcutoff, weights.beta.common, 
                                r.corr, sparseSigma, IsSingleVarinGroupTest, markerIDs, IsSparse, weights_specified,
                                method_to_CollapseUltraRare = "absence_or_presence",
                                MACCutoff_to_CollapseUltraRare = 10, DosageCutoff_for_UltraRarePresence = 0.5)
{
  # function_group_marker_list = group_info_list[[1]]; markerIDs = Gx$markerIDs; MAF_cutoff=c(0.001, 0.01); MACCutoff_to_CollapseUltraRare = 10;function_group_test=c("lof", "missense")
  
  obj.model$theta = tauVec
  obj.model$residuals = as.vector(y - obj.model$mu)
  adjustCCratioinGroupTest=FALSE
  if(traitType=="binary"){
    adjustCCratioinGroupTest=TRUE
  }
  
  ##########################
  # 	Process G, flipping 
  MACvec = Matrix::colSums(Gmat)
  m = ncol(Gmat)
  n = nrow(Gmat)
  AF = MACvec/(2 * n)
  flipindex = which(AF > 0.5)
  if (length(flipindex) > 0) {
    Gmat[, flipindex] = 2 - Gmat[, flipindex]
    MACvec[flipindex] = 2 * n - MACvec[flipindex]
    cat("Note the ", flipindex, "th variants were flipped to use dosages for the minor alleles in gene-based tests\n")
  }
  MAF = Matrix::colMeans(Gmat)/2
  flipInd = (AF > 0.5)
  
  if(traitType=="binary"){
    caseIndex = which(y==1)
    ctrlIndex = which(y==0)
    nCase = length(caseIndex)
    nCtrl = length(ctrlIndex)
  }
  #############################
  # Check method
  out.method = SKAT:::SKAT_Check_Method(method=method, r.corr=r.corr) 
  method = out.method$method
  r.corr = out.method$r.corr
  
  
  ###########################
  # Collapsing, collapsing by group
  re_group_id = Get_MultiSet_Id(markerIDs, function_group_marker_list, MACvec, MAF, 
                                function_group_test= function_group_test, MAF_cutoff=MAF_cutoff, 
                                MACCutoff_to_CollapseUltraRare = MACCutoff_to_CollapseUltraRare)
  
  
  re_collapsed = Get_Collapsed_Genotype(Gmat=Gmat, markerIDs=markerIDs, m=m, re_group_id=re_group_id, 
                                        function_group_test=function_group_test, DosageCutoff_for_UltraRarePresence=DosageCutoff_for_UltraRarePresence)
  
  ##############################
  #
  G1 = re_collapsed$Gmat
  
  #obj_cc = NULL;weights_specified=NULL;adjustCCratioinGroupTest=FALSE;mu=NULL;mu2=NULL;weightMAFcutoff=0.5
  
  re_phi_score = Get_Phi_Score(G1=G1, obj=obj.model, obj_cc=obj_cc, y = y, X = X, 
                               cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec,  
                               kernel= kernel, weights.beta.rare=weights.beta.rare, weights.beta.common=weights.beta.common, weightMAFcutoff = weightMAFcutoff,
                               sparseSigma = sparseSigma, weights_specified = weights_specified, adjustCCratioinGroupTest=adjustCCratioinGroupTest)
  
  
  Phi=NULL
  p.new=NULL
  if(adjustCCratioinGroupTest){
    Phi = re_phi_score$Phi_ccadj$val
    pval.new = re_phi_score$Phi_ccadj$p.new
  } else {
    Phi=re_phi_score$Phi 
  }
  
  
  #################################
  # Test for different groups
  
  re_test_gene_base<-list()
  re_test_single<-NULL
  k = 1
  for(i in 1:length(function_group_test)){
    collapse_index = re_collapsed$collapse_indicator[i]
    if(collapse_index < 0){
      collapse_index=NULL
    }
    
    for(j in  1:length(MAF_cutoff)){
      
      FuncMAF_marker_id = re_group_id$FuncMAF_list[[i]][[j]]
      if(length(FuncMAF_marker_id) > 0){
        SNP_Index = match(FuncMAF_marker_id, re_collapsed$markerIDs_new)
      } else {
        SNP_Index=NULL
      }
      # Add collapse variants
      SNP_Index<-c(collapse_index, SNP_Index) 
      index_test = setdiff(SNP_Index, re_phi_score$indexNeg)
      
      # Run the test
      re = Run_Genebase_Test(re_phi_score$Score, Phi, 
                                                 index_test=index_test, method=method, r.corr=r.corr, IsOutputBETASEinBurdenTest=IsOutputBETASEinBurdenTest)
      re_test_gene_base[[k]]<-list()
      re_test_gene_base[[k]]$re=re
      re_test_gene_base[[k]]$group = function_group_test[i]
      re_test_gene_base[[k]]$cutoff = MAF_cutoff[j]
      
      
      # Calculate MAC by case/control 
      MACg=0
      MAC_caseg=0
      MAC_ctrlg=0
      if(length(index_test)>0){
        MACg = sum(G1[,index_test])
        if(traitType=="binary"){
          MAC_caseg = sum(G1[caseIndex,index_test])
          MAC_ctrlg = sum(G1[ctrlIndex,index_test])
        }
      }
      re_test_gene_base[[k]]$outvec = c(re$p.value,  re$m, re$pval_Burden, 
                                        MACg, MAC_caseg, MAC_ctrlg)
      k = k+1
    }
  }
  if (IsSingleVarinGroupTest) {
    
    re_test_single = Run_Single_Test(re_phi_score$Score, Phi, re_phi_score$weights, pval=p.new)
    MACs =  Matrix::colSums(G1)
    MAFs = MACs/n/2
    MAC_cases=NA
    MAC_ctrls=NA
    MAF_cases=NA
    MAF_ctrls=NA
    
    if(traitType=="binary"){
      MAC_cases = Matrix::colSums(G1[caseIndex, ])
      MAC_ctrls = Matrix::colSums(G1[ctrlIndex, ])

      MAF_cases = MAC_case/nCase/2
      MAF_ctrls = MAC_ctrl/nCtrl/2
      
    } 
    re_test_single$out_vecs_df = data.frame(MarkerIDs=re_collapsed$markerIDs_new, Beta=re_test_single$Beta_single,
                          SE=re_test_single$SE_single, Pval= re_test_single$Pval_single, 
                          MACs=MACs, MAFs=MAFs,
                          MAC_cases=MAC_cases, MAC_ctrls=MAC_ctrls, MAF_cases=MAF_cases,MAF_ctrls=MAF_ctrls)
  }
  
  return(list(re_test_gene_base=re_test_gene_base, re_test_single=re_test_single))
  
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
Get_MultiSet_Id<-function(markerIDs, function_group_marker_list, MACvec, MAF, 
                           function_group_test=c("lof", "missense"), MAF_cutoff=c(0.0001, 0.001, 0.01), MACCutoff_to_CollapseUltraRare = 10){
	

  #function_group_marker_list = group_info_list[[1]]; markerIDs = Gx$markerIDs; MAF_cutoff=c(0.001, 0.01); MACCutoff_to_CollapseUltraRare = 10;function_group_test=c("lof", "missense")
  
	n_cutoff<-length(MAF_cutoff)
	n_group<-length(function_group_marker_list)
		
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
      GCollapsing<-cbind(Gnew, GCollapsing)
        	
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

    

Get_Weight<-function(MAF, weightMAFcutoff, weights.beta.rare, weights.beta.common, weights_specified=NULL){
  
  #print(MAF)
  if(is.null(weights_specified)){
    if(length(MAF) > 1){
      weights=rep(0,length(MAF))
      index1 = which(MAF<=weightMAFcutoff)
      if(length(index1) > 0) {weights[index1] = SKAT:::Beta.Weights(MAF[index1],weights.beta.rare)}
      index2 = which(MAF>weightMAFcutoff)
      if(length(index2) > 0) {weights[index2] = SKAT:::Beta.Weights(MAF[index2],weights.beta.common)}	
    }else{
      if(MAF<=weightMAFcutoff){
        weights = SKAT:::Beta.Weights(MAF,weights.beta.rare)
      }else{
        weights = SKAT:::Beta.Weights(MAF,weights.beta.common)
      }
    }
  }else{
    weights = weights_specified
    cat("weights is specified in the group file.\n")
  }
  
  return(weights)
  
}



#obj is the rda. file output from SAIGE step 1
#G1 is genotypes for testing gene, which contains m markers
Get_Phi_Score  = function(G1, obj, obj_cc, y, X, 
	cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec,  
	kernel= "linear.weighted", weights.beta.rare=c(1,25), weights.beta.common=c(0.5,0.5), weightMAFcutoff = 0.01,
	sparseSigma = NULL,  weights_specified = NULL, adjustCCratioinGroupTest=TRUE){

  #obj=obj.model; weights_specified=NULL;kernel= "linear.weighted"; weights.beta.rare=c(1,25); weights.beta.common=c(0.5,0.5); weightMAFcutoff = 0.01;
  obj.noK = obj$obj.noK
	obj.noK$y = y
  m = ncol(G1)
  n = nrow(G1)
	MACvec = Matrix::colSums(G1)
	MAF = MACvec/(2*n)
  AF = MAF
  mu=obj$mu
  mu2=obj$mu2
	

	weights = Get_Weight(MAF, weightMAFcutoff, weights.beta.rare, weights.beta.common, weights_specified)


	indexNeg = NULL
	MACvec_indVec_Zall = getCateVarRatio_indVec(MACvector=MACvec, cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude)
	GratioMatrixall = getGratioMatrix(MACvec_indVec_Zall, ratioVec)

    if (kernel == "linear.weighted") {
    	G1 = t(t(G1) * (weights[1:m]))
	}

	Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])

	Phi_ccadj=NULL
	if(is.null(obj$P)){
		G1_tilde_Ps_G1_tilde = getCovM_nopcg(G1=G1, G2=G1, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
		Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])
		
		if(adjustCCratioinGroupTest){
			Phi_ccadj = SPA_ER_kernel_related_Phiadj(G1_org, obj_cc, obj.noK, Cutoff=2, Phi, weights[1:m], VarRatio_Vec = as.vector(GratioMatrixall[1:m,1]), mu, sparseSigma)
		}
	}
	
	# Calculate indexNeg which identify Phi with very small Phi
	indexNeg = which(diag(as.matrix(Phi)) <= (.Machine$double.xmin)^(1/4)) 
	
	re = list(Score=Score, Phi=Phi, Phi_ccadj=Phi_ccadj, indexNeg=indexNeg, 
	MACvec_indVec_Zall=MACvec_indVec_Zall, GratioMatrixall=GratioMatrixall, MAF=MAF, weights=weights)
	return(re)

}

Run_Genebase_Test<-function(Score, Phi, index_test, method, r.corr, IsOutputBETASEinBurdenTest){

	#Score=re_phi_score$Score; Phi=re_phi_score$Phi; index_test=index_test
	m_test = length(index_test)
	if(m_test ==0){
		re = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, m_test=0)
		re$param = list(p.val.each = c(1, 1), rho=c(1,0))
		re$pval_Burden=1
		return(re)
	}


  Phi = Phi[index_test, index_test]
  Score = Score[index_test]
	
	re = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL))
	if(class(re) == "try-error"){
		re_btemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=1, method=method, Score.Resampling=NULL))
		re_stemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=0, method=method, Score.Resampling=NULL))
		
		if(class(re_btemp) == "try-error" | class(re_stemp) == "try-error"){	
			re = list(p.value = NA)
		}else{
			re = list(p.value = 2*min(re_btemp$p.value, re_stemp$p.value, 0.5), param = list())
			re$param = list(p.val.each = c(re_btemp$p.value, re_stemp$p.value), rho=c(1,0))

		}
	}
	

	#re$Phi_sum = sum(Phi)
	#re$Score_sum = sum(Score)
	#re$BETA_Burden  = re$Score_sum/re$Phi_sum

	re$IsMeta=TRUE
	re$m = m_test
	
	re$pval_Burden = NA
	re$SE_Burden = NA
	#P-value of burden
	if(!is.null(re$param$rho) ){
	  idx = which(re$param$rho > 0.999)
	  if(length(idx)>0){
	    re$pval_Burden = re$param$p.val.each[idx[1]]
	   # re$SE_Burden = abs(re$BETA_Burden/qnorm(re$pval_Burden/2))
	  }  
	} 

	return(re)	
	
}

Run_Single_Test<-function(Score, Phi, weights, pval=NULL){

	#Score =re_phi_score$Score; Phi=re_phi_score$Phi; weights; pval=NULL
  re<-list()
  m<-length(Score)
	re$Phi_single=diag(Phi)/(weights[1:m]^2)
  re$Score_single=Score/weights[1:m]
	re$Beta_single=(re$Score_single)/(re$Phi_single)
	
	if(is.null(pval)){
		re$Pval_single= pchisq((re$Score_single)^2/(re$Phi_single), lower.tail = FALSE, df=1, log.p=IsOutputlogPforSingle)
	} else {
		re$Pval_single = pval
	}
	if(IsOutputlogPforSingle){
		re$SE_single = abs((re$Beta_single)/(qnorm(re$Pval_single - log(2), log.p=T, lower.tail = F)))	
	}else{
		re$SE_single = abs((re$Beta_single)/qnorm(re$Pval_single/2))
	}		
	
	return(re)

}

				