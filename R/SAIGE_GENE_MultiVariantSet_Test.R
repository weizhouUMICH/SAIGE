

# Currently, conditional analysis is not supported
# Only absence_or_presence are supported
# idx_function_group has the most functional (ex. loss of function) to least functional


MultiSets_GroupTest = function (Gmat, obj.model, obj_cc, y, X, tauVec, traitType, 
                                function_group_marker_list, MAF_cutoff=c(0.0001, 0.001, 0.01), function_group_test=c("lof", "missense", "synonymous"),
                                cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, 
                                kernel, method, weights.beta.rare,  weightMAFcutoff, weights.beta.common, 
                                r.corr, sparseSigma, IsSingleVarinGroupTest, markerIDs, IsSparse, weights_specified,
                                method_to_CollapseUltraRare = "absence_or_presence",
                                MACCutoff_to_CollapseUltraRare = 10, DosageCutoff_for_UltraRarePresence = 0.5,
                                IsFastApprox=TRUE)
{
  # function_group_marker_list = group_info_list[[1]]; markerIDs = Gx$markerIDs; MAF_cutoff=c(0.001, 0.01); MACCutoff_to_CollapseUltraRare = 10;function_group_test=c("lof", "missense")
  

  obj.model$theta = tauVec
  obj.model$residuals = as.vector(y - obj.model$mu)
  adjustCCratioinGroupTest=FALSE
  if(traitType=="binary"){
    adjustCCratioinGroupTest=TRUE
  }
  
  
  ###########################
  # 	Process G, remove any variants outside of the analysis. 
  MaxMAF_Cutoff = max(MAF_cutoff)
  MACvec = Matrix::colSums(Gmat)
  m = ncol(Gmat)
  n = nrow(Gmat)
  
  MAF = MACvec/(2 * n)
  flipindex = which(MAF > 0.5)
  if (length(flipindex) > 0) {
  	MAF[flipindex] = 1-MAF[flipindex]
  }	
  
  idx_MAC0 = union(which(MACvec==0), which(MACvec==2*n))
  idx_OutRange = which(MAF > MaxMAF_Cutoff)
  idx_exclude_marker = union(idx_MAC0, idx_OutRange)
  
  if(length(idx_exclude_marker) > 0){
  	Gmat = Gmat[,-idx_exclude_marker, drop = F]
  	markerIDs=markerIDs[-idx_exclude_marker]
  	MACvec= MACvec[-idx_exclude_marker]
  	MAF = MAF[-idx_exclude_marker]
  }

  ##########################
  # 	Process G, flipping 

  flipindex = which(MACvec > n)
  if (length(flipindex) > 0) {
    Gmat[, flipindex] = 2 - Gmat[, flipindex]
    MACvec[flipindex] = 2 * n - MACvec[flipindex]
    cat("Note the ", flipindex, "th variants were flipped to use dosages for the minor alleles in gene-based tests\n")
  }

      
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
  
  #re_group_id1<<-re_group_id; markerIDs1<<-markerIDs; function_group_marker_list1<<-function_group_marker_list; MACvec1<-MACvec; MAF1<<-MAF
  
  re_collapsed = Get_Collapsed_Genotype(Gmat=Gmat, markerIDs=markerIDs, m=m, re_group_id=re_group_id, 
                                        function_group_test=function_group_test, DosageCutoff_for_UltraRarePresence=DosageCutoff_for_UltraRarePresence)
  
  #re_collapsed1<<-re_collapsed
  ##############################
  #
  G1 = re_collapsed$Gmat
  
  #obj_cc = NULL;weights_specified=NULL;adjustCCratioinGroupTest=FALSE;mu=NULL;mu2=NULL;weightMAFcutoff=0.5
  
  #re_phi_score = Get_Phi_Score(G1=G1, obj=obj.model, obj_cc=obj_cc, y = y, X = X, 
  #                             cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec,  
  #                             kernel= kernel, weights.beta.rare=weights.beta.rare, weights.beta.common=weights.beta.common, weightMAFcutoff = weightMAFcutoff,
  #                             sparseSigma = sparseSigma, weights_specified = weights_specified, adjustCCratioinGroupTest=adjustCCratioinGroupTest)
  
  re_phi_score = Get_Phi_Score(G1=G1, obj=obj.model, obj_cc=obj_cc, y = y, X = X, 
			cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec,    
			sparseSigma = sparseSigma, adjustCCratioinGroupTest=adjustCCratioinGroupTest, IsFastApprox=IsFastApprox)

  # beta weight 
  weights = Get_Weight(MAF=MAF, weightMAFcutoff=weightMAFcutoff, weights.beta.rare=weights.beta.rare, weights.beta.common=weights.beta.common, weights_specified=NULL)
 
  
  #re_phi_score1<<-re_phi_score 
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
                            index_test=index_test, method=method, r.corr=r.corr, IsOutputBETASEinBurdenTest=IsOutputBETASEinBurdenTest, weights=weights,
                            IsFastApprox=IsFastApprox)
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
    
    re_test_single = Run_Single_Test(re_phi_score$Score, Phi, pval=p.new)
    MACs =  Matrix::colSums(G1)
    MAFs = MACs/n/2
    MAC_cases=NA
    MAC_ctrls=NA
    MAF_cases=NA
    MAF_ctrls=NA
    
    if(traitType=="binary"){
      MAC_cases = Matrix::colSums(G1[caseIndex, ])
      MAC_ctrls = Matrix::colSums(G1[ctrlIndex, ])

      MAF_cases = MAC_cases/nCase/2
      MAF_ctrls = MAC_ctrls/nCtrl/2
      
    } 
    re_test_single$out_vecs_df = data.frame(MarkerIDs=re_collapsed$markerIDs_new, Beta=re_test_single$Beta_single,
                          SE=re_test_single$SE_single, Pval= re_test_single$Pval_single, 
                          MACs=MACs, MAFs=MAFs,
                          MAC_cases=MAC_cases, MAC_ctrls=MAC_ctrls, MAF_cases=MAF_cases,MAF_ctrls=MAF_ctrls)
  }
  
  return(list(re_test_gene_base=re_test_gene_base, re_test_single=re_test_single))
  
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

	#Gmat<-Gmat1; markerIDs<-markerIDs1; re_group_id<-re_group_id1; function_group_test=c("lof", "missense", "synonymous");DosageCutoff_for_UltraRarePresence=0.5
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
    Gmat = cbind(GCollapsing)
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
# Weight is not applied to....
Get_Phi_Score  = function(G1, obj, obj_cc, y, X, 
	cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec,  
	sparseSigma = NULL,  adjustCCratioinGroupTest=TRUE, IsFastApprox=FALSE){

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
	

	indexNeg = NULL
	MACvec_indVec_Zall = getCateVarRatio_indVec(MACvector=MACvec, cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude)
	GratioVecall = getGratioVector(MACvec_indVec_Zall, ratioVec)
	
	Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])

	Phi_ccadj=NULL
	if(is.null(obj$P)){
		G1_tilde_Ps_G1_tilde = getCovM_nopcg_fast(G1=G1, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2, IsFastApprox=IsFastApprox)
		GratioVecall_sqrt = sqrt(GratioVecall)
		#G1_tilde_Ps_G1_tilde1<<-G1_tilde_Ps_G1_tilde; GratioVecall_sqrt1<<-GratioVecall_sqrt
		Phi = t(t(G1_tilde_Ps_G1_tilde * GratioVecall_sqrt) * GratioVecall_sqrt)
		
		if(adjustCCratioinGroupTest){
			#G2<<-G1; obj_cc<<-obj_cc; obj.noK<<-obj.noK; Phi11<<-Phi; VarRatio_Vec<<-GratioVecall; mu<<-mu; sparseSigma<<-sparseSigma
			
			Phi_ccadj= SPA_ER_kernel_related_Phiadj_fast(G1, obj_cc, obj.noK, Cutoff=2, Phi, Score, VarRatio_Vec= GratioVecall, mu, sparseSigma)
			
			#Phi_ccadj = SPA_ER_kernel_related_Phiadj_fast(G1, obj_cc, obj.noK, Cutoff=2, Phi, rep(1,m), VarRatio_Vec = as.vector(GratioMatrixall[1:m,1]), mu, sparseSigma)
			#Phi_ccadj = SPA_ER_kernel_related_Phiadj(G1, obj_cc, obj.noK, Cutoff=2, Phi, rep(1,m), VarRatio_Vec = as.vector(GratioMatrixall[1:m,1]), mu, sparseSigma)
		}
	}
	
	
	
	# Calculate indexNeg which identify Phi with very small Phi
	indexNeg = which(diag(as.matrix(Phi)) <= (.Machine$double.xmin)^(1/4)) 
	
	re = list(Score=Score, Phi=Phi, Phi_ccadj=Phi_ccadj, indexNeg=indexNeg, 
	MACvec_indVec_Zall=MACvec_indVec_Zall, GratioVecall=GratioVecall, MAF=MAF)
	return(re)

}

Run_SKAT_Burden_Only<-function(Score, Phi, method, IsFastApprox=TRUE){
	
	re = list(p.value = NA)
	re_btemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=1, method=method, Score.Resampling=NULL, IsFast=IsFastApprox))
	re_stemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=0, method=method, Score.Resampling=NULL, IsFast=IsFastApprox))
		
	if(class(re_btemp) == "try-error" | class(re_stemp) == "try-error"){	
		re = list(p.value = NA)
	}else{
		re = list(p.value = 2*min(re_btemp$p.value, re_stemp$p.value, 0.5), param = list())
		re$param = list(p.val.each = c(re_btemp$p.value, re_stemp$p.value), rho=c(1,0))
	}

	return(re)

}

Run_Genebase_Test<-function(Score, Phi, index_test, method, r.corr, IsOutputBETASEinBurdenTest, weights=NULL, IsFastApprox=TRUE){

	#Score=re_phi_score$Score; Phi=re_phi_score$Phi; index_test=index_test
	
	m_test = length(index_test)
	if(m_test ==0){
		re = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, m_test=0)
		re$param = list(p.val.each = c(1, 1), rho=c(1,0))
		re$pval_Burden=1
		return(re)
	}


  	Phi_test = Phi[index_test, index_test]
  	Score_test = Score[index_test]
	
	if(!is.null(weights)){
		weights1 = weights[index_test]
		Phi_test = t(t(Phi_test * weights1)*weights1)
		Score_test = Score_test * weights1
	}
	
	p1 = length(Score_test)
	
	# When p1 > 1000, only run Burden and SKAT
	if(IsFastApprox && p1 > 1000){
		re = Run_SKAT_Burden_Only(Score=Score_test, Phi=Phi_test, method=method, IsFastApprox=IsFastApprox)
	} else {
		re = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_test, Phi=Phi_test, r.corr=r.corr, method=method, Score.Resampling=NULL))
		if(class(re) == "try-error"){
			re = Run_SKAT_Burden_Only(Score=Score_test, Phi=Phi_test, method=method, IsFastApprox=IsFastApprox)
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
	
	if(m_test==1 && is.na(re$pval_Burden)){
		re$pval_Burden = re$p.value
	}

	return(re)	
	
}

Run_Single_Test<-function(Score, Phi, pval=NULL, IsOutputlogPforSingle=FALSE){

	#Score =re_phi_score$Score; Phi=re_phi_score$Phi; weights; pval=NULL
	re<-list()
	m<-length(Score)
	re$Phi_single=diag(Phi)
	re$Score_single=Score
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


#obj is the rda. file output from SAIGE step 1
#G1 is genotypes for testing gene, which contains m markers
Get_Phi_Score_OLD  = function(G1, obj, obj_cc, y, X, 
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
	
	#GratioMatrixall1<<-GratioMatrixall; m1<<-m; MACvec_indVec_Zall1<<-MACvec_indVec_Zall; ratioVec1<<-ratioVec
    if (kernel == "linear.weighted") {
    	G1 = t(t(G1) * (weights[1:m]))
	}

	Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])

	Phi_ccadj=NULL
	if(is.null(obj$P)){
		G1_tilde_Ps_G1_tilde = getCovM_nopcg(G1=G1, G2=G1, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
		Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])
		
		if(adjustCCratioinGroupTest){
			Phi_ccadj = SPA_ER_kernel_related_Phiadj(G1, obj_cc, obj.noK, Cutoff=2, Phi, rep(1,m), VarRatio_Vec = as.vector(GratioMatrixall[1:m,1]), mu, sparseSigma)
		}
	}
	
	# Calculate indexNeg which identify Phi with very small Phi
	indexNeg = which(diag(as.matrix(Phi)) <= (.Machine$double.xmin)^(1/4)) 
	
	re = list(Score=Score, Phi=Phi, Phi_ccadj=Phi_ccadj, indexNeg=indexNeg, 
	MACvec_indVec_Zall=MACvec_indVec_Zall, GratioMatrixall=GratioMatrixall, MAF=MAF, weights=weights)
	return(re)

}

				