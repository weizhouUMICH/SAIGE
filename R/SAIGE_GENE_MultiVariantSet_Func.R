


splitfun_weight = function(x) {
    return(strsplit(x, split = ";")[[1]][2])
}
splitfun_markerID = function(x) {
    return(strsplit(x, split = ";")[[1]][1])
}
    

Check_File_Exist<-function(file, filetype){
    if (!file.exists(file)) {
        stop("ERROR! ", filetype , file, " does not exsit\n")
    }
}

Check_Genotypes_and_Samples<-function(bgenFile, vcfFile, vcfField, vcfFileIndex, savFile, sampleFile, sampleID, chrom){

    # genotype data
    if (bgenFile != "") {
        Check_File_Exist(bgenFile, "bgenFile")
        dosageFileType = "bgen"
    }
    else if (vcfFile != "") {
        Check_File_Exist(vcfFile, "vcfFile")
        
        if (!grepl("\\.sav$", vcfFile) && !file.exists(paste(vcfFile, ".csi", sep = ""))) {
            stop("ERROR! vcfFileIndex ", paste(vcfFile, ".csi", sep = ""), " does not exist\n")
        }
        dosageFileType = "vcf"
        if (chrom == "") {
            stop("ERROR! chrom needs to be specified for the vcf file\n")
        }
    }
    else if (savFile != "") {
        Check_File_Exist(savFile, "savFile")
        vcfFile = savFile
        dosageFileType = "vcf"
    }
	
	  # Check sample file
    sampleListinDosage = NULL
    if (!file.exists(sampleFile)) {
        if (dosageFileType == "bgen") {
            stop("ERROR! The dosage file type is bgen but sampleFile ", sampleFile, " does not exsit\n")
        }
    }
    else {
        sampleListinDosage = data.frame(data.table:::fread(sampleFile, header = F, stringsAsFactors = FALSE, colClasses = c("character")))
        sampleListinDosage$IndexDose = seq(1, nrow(sampleListinDosage), by = 1)
        cat(nrow(sampleListinDosage), " sample IDs are found in sample file\n")
        colnames(sampleListinDosage)[1] = "IIDDose"
    }

	# SampleInModel 
    sampleInModel = NULL
    sampleInModel$IID = sampleID
    sampleInModel = data.frame(sampleInModel)
    sampleInModel$IndexInModel = seq(1, length(sampleInModel$IID), by = 1)
    cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")

	# When VCF
   	if (dosageFileType == "vcf") {
        vcffileopen = FALSE
        isVariant = SAIGE:::setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
        sampleListinDosage_vec = SAIGE:::getSampleIDlist_vcfMatrix()
        
        if (is.null(sampleListinDosage)) {
            sampleListinDosage = data.frame(IIDDose = sampleListinDosage_vec)
            sampleListinDosage$IndexDose = seq(1, nrow(sampleListinDosage), by = 1)
            cat(nrow(sampleListinDosage), " sample IDs are found in the vcf file\n")
        }
    }
    dataMerge = merge(sampleInModel, sampleListinDosage, by.x = "IID", by.y = "IIDDose")
    dataMerge_sort = dataMerge[with(dataMerge, order(IndexInModel)), ]

  	if (nrow(dataMerge_sort) < nrow(sampleInModel)) {
        stop("ERROR!", nrow(sampleInModel) - nrow(dataMerge_sort), " samples used in glmm model fit do not have dosages\n")
    } else {
        dataMerge_v2 = merge(dataMerge_sort, sampleListinDosage, by.x = "IID", by.y = "IIDDose", all.y = TRUE)
        dataMerge_v2_sort = dataMerge_v2[with(dataMerge_v2, order(IndexDose.y)), ]
        sampleIndex = dataMerge_v2_sort$IndexInModel
        N = sum(!is.na(sampleIndex))
        cat(N, " samples were used in fitting the NULL glmm model and are found in sample file\n")
        sampleIndex[is.na(sampleIndex)] = -10
        sampleIndex = sampleIndex - 1
        rm(dataMerge)
        rm(dataMerge_v2)
        rm(dataMerge_sort)
        rm(dataMerge_v2_sort)
    }
    
    rm(sampleInModel)


	re<-list(dosageFileType=dosageFileType, vcfFile=vcfFile, sampleIndex=sampleIndex)
	return(re)
	
}


Get_Variance_Ratio<-function(varianceRatioFile, sparseSigmaFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude){

	ratioVec = c(1)
	
    # check variance ratio
    if (!file.exists(varianceRatioFile)) {
        if (sparseSigmaFile == "") {
            stop("ERROR! varianceRatioFile ", varianceRatioFile, " does not exsit but sparseSigmaFile also does not exist \n")
        }
        else {
            cat("varianceRatioFile is not specified so variance ratio won't be used\n")
        }
        ratioVec = c(1)
    }
    else {
        varRatioData = data.frame(data.table:::fread(varianceRatioFile, header = F, stringsAsFactors = FALSE))
        ln = length(cateVarRatioMinMACVecExclude)
        hn = length(cateVarRatioMaxMACVecInclude)
        if (nrow(varRatioData) == 1) {
            stop("ERROR! To perform gene-based tests, categorical variance ratios are required\n")
        }
        else {
            ratioVec = varRatioData[, 1]
            nrv = length(ratioVec)
            if (nrv != ln) {
                stop("ERROR! The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude\n")
            }
            if (ln != (hn + 1)) {
                stop("ERROR! The length of cateVarRatioMaxMACVecInclude does not match with the lenght of cateVarRatioMinMACVecExclude (-1)\n")
            }
        }
        cat("variance Ratio is ", ratioVec, "\n")
    }
    
    return(ratioVec)
    
}


getGratioVector<-function(MACvec_indVec, ratioVec) {

    numCate = length(ratioVec)
    if (numCate > 1) {
        indMatrix = contr.sum(numCate, contrasts = FALSE)
        GindMatrix = NULL
        for (i in MACvec_indVec) {
            GindMatrix = rbind(GindMatrix, indMatrix[i, ])
        }
        GratioVec = GindMatrix %*% matrix(ratioVec, ncol = 1)
      
    }
    else {
        GratioVec = rep(ratioVec[1], length(MACvec_indVec))
    }
    
    GratioVec = as.vector(  GratioVec )
    return(GratioVec)
}


getCovM_nopcg_fast<-function(G1, XV, XXVX_inv, sparseSigma=NULL, mu2, IsFastApprox=TRUE){

	#G21=G2
	#G1=G2; XV=obj.noK$XV; XXVX_inv=obj.noK$XXVX_inv; mu2 = mu21
	# XV<-obj.noK$XV; XXVX_inv<-obj.noK$XXVX_inv
 	if(!IsFastApprox){
 		re = SAIGE:::getCovM_nopcg(G1=G1, G2=G1, XV=XV, XXVX_inv=XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
 		return(re)
 	}
        nSNP1<-ncol(G1)
	n<-nrow(G1)
	
	XV_G1 = XV %*% G1
	A_G1<-XXVX_inv %*% XV_G1
	
    if(!is.null(sparseSigma)){ 
    	SI_XXVX_inv = solve(sparseSigma, XXVX_inv ,sparse=TRUE)
    	SI_G1<-solve(sparseSigma, G1, sparse = TRUE)
    	SI_A_G1<-SI_XXVX_inv %*% XV_G1

		var1<-colSums(G1 * SI_G1) - colSums(G1 *SI_A_G1 ) *2 + colSums(A_G1 *SI_A_G1)
		G1_sum  = colSums(G1)
    	G1_cov<-crossprod(G1, G1) -  G1_sum%*% t(G1_sum)/n
    	
    	diag_cov<-diag(G1_cov)
    	G1_cor<-t(t(G1_cov/sqrt(diag_cov))/sqrt(diag_cov))
    	Mat<- t(t(G1_cor * sqrt(var1))* sqrt(var1))

    } else {
    	G2=G1
    	G2 = G2 * mu2
        XV_G2 = XV %*% G2
        SI_A_G2<-XXVX_inv %*% XV_G2
        A1<- t(G1) %*% (G2 - SI_A_G2)
        A2<- t(XXVX_inv) %*% (G2 - SI_A_G2)
        Mat<-(A1 - t(XV_G1) %*% A2)
        Mat<-as.matrix(Mat)
    	                 
    }
    
    return(Mat)
	

}




#
#	Some changes in the function
#		Use Score instead of calculating q and mu1
#
#	
SPA_ER_kernel_related_Phiadj_fast<-function(G, obj, obj.noK, Cutoff=2, Phi, Score, VarRatio_Vec, mu.a, sparseSigma){
	
	#G<-G2; obj<-obj_cc; obj.noK<-obj.noK; Phi<-Phi1; weight<-rep(1,ncol(G2)); mu.a<-mu; Cutoff=2

	p.m<-ncol(G)
	n<-nrow(G)
	zscore.all_0<-rep(0, p.m)
	zscore.all_1<-rep(0, p.m)

	MAFsum=Matrix::colSums(G)
	mu2.a = mu.a *(1-mu.a)
	
	VarS_org=diag(Phi)	
	stat.qtemp =Score^2/VarS_org
	p.new = pchisq(stat.qtemp, lower.tail = FALSE, df = 1)  
	zscore.all_0 = Score
		
	#XXVX_XV_G = obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G)
	#q<-colSums(G *obj.noK$y) - colSums(XXVX_XV_G *obj.noK$y )
	#mu.qtemp=mu.a;    
	#mu1 <- colSums(mu.qtemp * G) - colSums(mu.qtemp * XXVX_XV_G)

 	#stat.qtemp<-(q - mu1)^2/VarS_org
    #p.new<-pchisq(stat.qtemp, lower.tail = FALSE, df = 1)  
	#zscore.all_0=(q-mu1)  ##sum(G[,jj]*(obj.noK$y-mu.a))  		
	id1<-which(stat.qtemp > Cutoff^2)
	if(length(id1)> 0){
	
		for(jj in id1){
			if (MAFsum[jj]<=10){
				p_temp1=SKAT::SKATBinary((G[,jj,drop=F]),obj, method.bin="Hybrid")$p.value
			} else {
			
				n.g<-MAFsum[jj]
				NAset<-which(G[,jj]==0)
				
				AC = MAFsum[jj]
				AF = AC/n
				MAF = AF
		
				p_temp1=scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G[,jj], AC, AF, MAF, IsSparse=TRUE, 
					obj.noK, mu.a = mu.a, mu2.a = mu2.a, obj.noK$y, varRatio=VarRatio_Vec[jj], Cutoff = Cutoff, 
					rowHeader=NULL, sparseSigma=sparseSigma)$p.value
				p_temp1 = unlist(p_temp1)[1]
			
			}
			p.new[jj]=p_temp1
		}
	}
	
	
	idx_0<-which(VarS_org >0)
	idx_p0<-which(p.new >0)
	idx_p1<-which(p.new <0)
	
	if(length(idx_0) > 0){
		zscore.all_1[idx_0]= qnorm(p.new[idx_0]/2, mean = 0, sd =sqrt( VarS_org[idx_0]),lower.tail = FALSE, log.p = FALSE)*sign(Score)
	}
	VarS = zscore.all_0^2/500
	if(length(idx_p0) > 0){
		VarS[idx_p0]= zscore.all_0[idx_p0]^2/qchisq(p.new[idx_p0], 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
	}	 
	

	###################################
	# Two different types of burden comparison
		
	
	vars_inf=which(VarS==Inf)
	if (length(vars_inf)>0){
		VarS[vars_inf] = 0
		zscore.all_1[vars_inf]=0
		Phi[vars_inf,]=0
		Phi[,vars_inf]=0
	}
	
	scaleFactor = sqrt(VarS/VarS_org)	

	###################################
	# Burden test
		
	G_Burden = rowSums(G)
	g.sum<-G_Burden  - obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G_Burden)
	q.sum<-sum(Score)
	
	p.value_burden<-SPAtest:::Saddle_Prob(q.sum , mu=mu.a, g=g.sum, Cutoff=2,alpha=2.5*10^-6)$p.value

	
	###################################
	# Compare two burden test approach and calculate adjust ratio r
	
	G2_adj_n=t(t(Phi * sqrt(VarS/VarS_org)) * sqrt(VarS/VarS_org))
	v1=rep(1,nrow(G2_adj_n))
	VarQ = sum(G2_adj_n)
	Q_b=sum(zscore.all_1)^2
	
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



Get_Results_DF<-function(groupTestResult, geneID){
  
  #groupTestResult1<<-groupTestResult
  if(is.null(groupTestResult)){
  	return(list(gene_base_test_df=NULL, single_test_df=NULL))
  }
  
  re_test_gene_base = groupTestResult$re_test_gene_base
  re_test_single = groupTestResult$re_test_single
  
  nSets<-length(re_test_gene_base)
  group.a<-rep("", nSets)
  cutoff.a<-rep("", nSets)
  outvecs<-NULL
  for(i in 1:nSets){
  	outvec = re_test_gene_base[[i]]$outvec
    group.a[i]<-re_test_gene_base[[i]]$group
    cutoff.a[i]<-re_test_gene_base[[i]]$cutoff
    outvecs = rbind(outvecs, outvec)
  }
  
  pval = outvecs[,1]
  pval_burden = outvecs[,3]
  idx1<-which(!is.na(pval))
  idx2<-which(!is.na(pval_burden))  
  
  
  cc_pval=NA
  if(length(idx1)> 0){
  	cc_pval = SAIGE:::CCT(pval[idx1])
  }
  
  cc_pval_burden=NA
  if(length(idx2)> 0){
  	cc_pval_burden = SAIGE:::CCT(pval_burden[idx2])
  }
  
    
  cc_vec = c(cc_pval, NA, cc_pval_burden, NA, NA, NA)
  outvecs = rbind(outvecs, cc_vec)
  gene_base_test_df = data.frame(GeneID = geneID, group=c(group.a, "CC-all"), cutoff=c(cutoff.a, "CC-all"), outvecs=outvecs )
  
  ## Add CC results
  
    
  #re$p.value,  re$m, re$BETA_Burden, re$SE_Burden, re$pval_Burden,  MACg, MAC_caseg, MAC_ctrlg)
  colnames(gene_base_test_df)<-c("GeneID", "FuncGroup", "Cutoff", "pval", "m",
                                 "pval_Burden", "MAC", "MAC_case", "MAC_ctrl")
  
  single_test_df = data.frame(GeneID=geneID, re_test_single$out_vecs_df )
  return(list(gene_base_test_df=gene_base_test_df, single_test_df=single_test_df))
}



