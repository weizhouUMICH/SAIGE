setSparseSigma = function(sparseSigmaFile){

  Check_File_Exist(sparseSigmaFile, "sparseSigmaFile")
  sparseSigma = Matrix:::readMM(sparseSigmaFile)
  locations = rbind(sparseSigma@i, sparseSigma@j)
  values = sparseSigma@x
  nSubj = dim(sparseSigma)[1]
  sigmaMatListR = list(locations = locations,
                     values = values,
                     nSubj = nSubj)
  return(sigmaMatListR)	
}	




##working
SAIGE.Region = function(objNull,
			objGeno,
			sparseSigma,
			OutputFile,
			method_to_CollapseUltraRare,
                        MACCutoff_to_CollapseUltraRare,
                        DosageCutoff_for_UltraRarePresence, 
			groupFile, 
			annolist, 
			maxMAFlist,
		        max_markers_region,	
			genoType, 
			markerInfo,
			traitType,
			isImputation,
			isCondition, 
			r.corr){
  OutputFileIndex = NULL	
  if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index")

  outList = checkOutputFile(OutputFile, OutputFileIndex, "Region", 1) # Check 'Util.R'

  indexChunk = outList$indexChunk
  Start = outList$Start
  End = outList$End	

  if(End)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results are saved in '", OutputFile, "'. ",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    return(message)
  }

  n = length(objNull$y) #sample size 

  ## annotation in region
  ##need to revise
  if(is.null(r.corr)){
    out.method = SKAT:::SKAT_Check_Method(method="optimal.adj", r.corr=0)
    method=out.method$method
    r.corr=out.method$r.corr
    cat("SKAT-O test will be performed. P-values for BUTDEN and SKAT tests will also be output\n.")
    regionTestType = "SKAT-O"
  }else if(r.corr == 0){
    method = NULL
    cat("SKAT test will be performed\n")
    regionTestType = "SKAT"
  }else if(r.corr == 1){	  
    method = NULL
    cat("BURDEN test will be performed\n")
    regionTestType = "BURDEN"
  }


  RegionList = SAIGE.getRegionList(groupFile, annolist, markerInfo)
  nRegions = length(RegionList)


  P1Mat = matrix(0, max_markers_region, n);
  P2Mat = matrix(0, n, max_markers_region);



  chrom1 = "FakeCHR";

  for(i in (indexChunk+1):nRegions){

    pval.Region = NULL
    region = RegionList[[i]]
    regionName = names(RegionList)[i]
    annolist = region$annoVec 
    if(length(annolist) == 0){
      next
    }	    
    SNP = region$SNP
    if(genoType != "vcf"){
    	genoIndex = region$genoIndex
    }else{
        SNPlist = paste(c(regionName, SNP), collapse = "\t") 
        set_iterator_inVcf(SNPlist, chrom1, 1, 200000000)
        isVcfEnd =  check_Vcf_end()
    	if(!isVcfEnd){
		genoIndex = rep(-1, length(SNP))
    	}else{
        	warning("No markers in region ", regionName, " are found in the VCF file")
		next
    	}	
    }	    
    annoIndicatorMat = region$annoIndicatorMat
    chrom = region$chrom

    print(paste0("Analyzing Region of ", regionName, " (",i,"/",nRegions,")."))
    #print(paste(SNP, collapse = ", "))

    outList0 = mainRegion(genoType, genoIndex, annoIndicatorMat, maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, isImputation, annolist, isCondition, regionTestType)
    if(length(outList0) == 0){
      next
    } 
    outList = outList0$outList
    #print("outList$VarMat")
    #print(outList$VarMat)

    if(isCondition){
      weightMat_G2_G2 = outList$G2_Weight_cond %*% t(outList$G2_Weight_cond)
    }	  
    info.Region = outList0$info.Region
    ### Get annotation maf indicators
    annoMAFIndicatorMat = outList$annoMAFIndicatorMat
    #print("dim(annoMAFIndicatorMat)")
    #print(dim(annoMAFIndicatorMat))
    #annoMAFIndicatorMat = annoMAFIndicatorMat[complete.cases(annoMAFIndicatorMat),,drop=F]
    #annoMAFIndicatorMat = annoMAFIndicatorMat[which(rowSums(annoMAFIndicatorMat) != 0),,drop=F]
    #print(dim(annoMAFIndicatorMat))

    ### 3. Adjust for saddlepoint approximation
    notNAindice = which(!is.na(outList$TstatVec))
    StatVec = outList$TstatVec[notNAindice]
    #outList$VarMat = outList$VarMat[notNAindice, notNAindice, drop=F]
    if(regionTestType == "BURDEN"){
	VarSVec = outList$VarVec
    }else{	    
        VarSVec = diag(outList$VarMat)
    }
    VarSVec = VarSVec[!is.na(VarSVec)]
    adjPVec = outList$pvalVec
    adjPVec = adjPVec[!is.na(adjPVec)]
    
    varTestedIndices = which(rowSums(annoMAFIndicatorMat) > 0)
    annoMAFIndicatorMat = annoMAFIndicatorMat[varTestedIndices, , drop=F]
    MAFVec = outList$MAFVec[varTestedIndices]
    #print("MAFVec")
    #print(MAFVec)
    AnnoWeights = dbeta(MAFVec,1,25)
    weightMat = AnnoWeights %*% t(AnnoWeights)
    if(isCondition){    
    	weightMat_G1_G2 = AnnoWeights %*% t(outList$G2_Weight_cond)
    }	
    wStatVec = StatVec * AnnoWeights

    if(regionTestType == "BURDEN"){
      wadjVarSVec = outList$VarVec * AnnoWeights
      if(isCondition){
	  wStatVec_cond = wStatVec - outList$TstatAdjCond
          wadjVarSVec_cond = wadjVarSVec - diag(outList$VarVecAdjCond)
      }
    }else{	    
      wadjVarSMat = outList$VarMat * weightMat
      if(isCondition){
	  wStatVec_cond = wStatVec - outList$TstatAdjCond
          wadjVarSMat_cond = wadjVarSMat - outList$VarMatAdjCond	
      }
    }



    annoMAFIndVec = c()
    for(j in 1:length(annolist)){
	AnnoName = annolist[j]    
	for(m in 1:length(maxMAFlist)){
		jm = (j-1)*(length(maxMAFlist)) + m
		maxMAFName = maxMAFlist[m]
		tempPos = which(annoMAFIndicatorMat[,jm] != 0)
#	        print("tempPos")
#	        print(tempPos)
#		print("length(tempPos)")
#	        print(length(tempPos))
	       if(length(tempPos) > 0){
		annoMAFIndVec = c(annoMAFIndVec, jm)
	        if(regionTestType == "BURDEN"){
			Phi = wadjVarSVec[tempPos]
	        }else{		
			Phi = wadjVarSMat[tempPos, tempPos, drop=F]
		}
		Score = wStatVec[tempPos]
		if(traitType == "binary"){
			p.new = adjPVec[tempPos]
			g.sum = outList$genoSumMat[,jm]
			q.sum<-sum(outList$gyVec[tempPos] * AnnoWeights[tempPos])
			mu.a = objNull$mu
			re_phi = get_newPhi_scaleFactor(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType)
		        Phi = re_phi$val
                }
#	        print("r.corr")
#		print(r.corr)


		groupOutList = get_SKAT_pvalue(Score, Phi, r.corr, regionTestType)
#		print("Score")
#		print(Score)
#		 print("Phi")
 #               print(Phi)
#		print("groupOutList")
#		print(groupOutList)
		resultDF = data.frame(Region = regionName,
                                                    Group = AnnoName,
                                                    max_MAF = maxMAFName,
                                                    Pvalue = groupOutList$Pvalue_SKATO,
                                                    Pvalue_Burden = groupOutList$Pvalue_Burden,
                                                    Pvalue_SKAT = groupOutList$Pvalue_SKAT,
                                                    BETA_Burden = groupOutList$BETA_Burden,
                                                    SE_Burden = groupOutList$SE_Burden)
#		print("Phi")
#		print(Phi)
	      if(isCondition){
		if(traitType == "binary"){
			#weightMat_G1_G2_sub = weightMat_G1_G2[tempPos,]
			#weightMat_sub = weightMat[tempPos,tempPos]
			G1tilde_P_G2tilde_Mat_scaled = t(t((outList$G1tilde_P_G2tilde_Weighted_Mat[tempPos,,drop=F]) * sqrt(as.vector(re_phi$scaleFactor))) * sqrt(as.vector(outList$scalefactor_G2_cond)))
#t(t(b * sqrt(a1)) * sqrt(a2))
			#G1tilde_P_G2tilde_Mat_scaled = diag(as.vector(re_phi$scaleFactor)) %*% (outList$G1tilde_P_G2tilde_Weighted_Mat) %*% diag(as.vector(outList$scalefactor_G2_cond))
			#VarInvMat_cond_scaled = diag(as.vector(1/(outList$scalefactor_G2_cond))) %*% outList$VarInvMat_G2_cond %*% diag(as.vector(1/(outList$scalefactor_G2_cond))) * (1/weightMat_G2_G2)
			#print("dim(G1tilde_P_G2tilde_Mat_scaled)")
			#print(dim(G1tilde_P_G2tilde_Mat_scaled))
			#print("dim(outList$VarInvMat_G2_cond_scaled)")
			#print(dim(outList$VarInvMat_G2_cond_scaled))
		        adjCondTemp = G1tilde_P_G2tilde_Mat_scaled %*% outList$VarInvMat_G2_cond_scaled	
			#print("dim(adjCondTemp)")
			#print(dim(adjCondTemp))
		#outList$VarMatAdjCond = adjCondTemp %*% t(G1tilde_P_G2tilde_Mat_scaled)
			VarMatAdjCond = adjCondTemp %*% t(G1tilde_P_G2tilde_Mat_scaled)
		#outList$TstatAdjCond = adjCondTemp %*% (outList$Tstat_G2_cond * outList$G2_Weight_cond)
			TstatAdjCond = adjCondTemp %*% (outList$Tstat_G2_cond * outList$G2_Weight_cond)
			Phi_cond = re_phi$val - diag(VarMatAdjCond)
		#wadjVarSMat = re_phi$val
			Score_cond = Score - TstatAdjCond
		#wStatVec = wStatVec - outList$TstatAdjCond
	        #wadjVarSMat = wadjVarSMat - outList$VarMatAdjCond

			
		}else{
			Score_cond = wStatVec_cond[tempPos]
			Phi_cond = wadjVarSMat_cond[tempPos, tempPos]
		}


			

		groupOutList_cond = get_SKAT_pvalue(Score_cond, Phi_cond, r.corr, regionTestType)

		resultDF$Pvalue_cond = groupOutList_cond$Pvalue_SKATO
		resultDF$Pvalue_Burden_cond = groupOutList_cond$Pvalue_Burden
		resultDF$Pvalue_SKAT_cond = groupOutList_cond$Pvalue_SKAT
		resultDF$BETA_Burden_cond = groupOutList_cond$BETA_Burden
		resultDF$SE_Burden_cond = groupOutList_cond$SE_Burden
	    }
		pval.Region = rbind.data.frame(pval.Region, resultDF)

	   }  
	}
    }
    pval.Region$MAC = outList$MAC_GroupVec[annoMAFIndVec]
    if(traitType == "binary"){
    	pval.Region$MAC_case = outList$MACCase_GroupVec[annoMAFIndVec]
    	pval.Region$MAC_control = outList$MACCtrl_GroupVec[annoMAFIndVec]
    }
    pval.Region$Number_rare = outList$NumRare_GroupVec[annoMAFIndVec]
    pval.Region$Number_ultra_rare = outList$NumUltraRare_GroupVec[annoMAFIndVec]
	
   if(sum(!is.na(pval.Region$Pvalue)) > 0){
   ##Combined using the Cauchy combination
   #pvals = pval.Region$Pvalue
   #pvals = pvals[!is.na(pvals)]
   #cctpval = CCT(pvals)
   cctpval = get_CCT_pvalue(pval.Region$Pvalue
   #pvals = pval.Region$Pvalue_Burden
   #pvals = pvals[!is.na(pvals)]
   #cctpval_Burden = CCT(pvals) 
   cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden)
   #pvals = pval.Region$Pvalue_SKAT
   #pvals = pvals[!is.na(pvals)]
   #cctpval_SKAT = CCT(pvals)
   cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT)

   cctVec = c(regionName, "Cauchy", NA, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA)
   if(isCondition){
        #pvals = pval.Region$Pvalue_cond
   	#pvals = pvals[!is.na(pvals)]
   	#cctpval = CCT(pvals)
   	#pvals = pval.Region$Pvalue_Burden_cond
   	#pvals = pvals[!is.na(pvals)]
   	#cctpval_Burden = CCT(pvals)
   	#pvals = pval.Region$Pvalue_SKAT_cond
   	#pvals = pvals[!is.na(pvals)]
   	#cctpval_SKAT = CCT(pvals)
	cctpval = get_CCT_pvalue(pval.Region$Pvalue_cond
   	cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)
	cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)

	 cctVec = c(cctVec, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA)
   }
   cctVec = c(cctVec, NA)
   #cctVec = c(regionName, "Cauchy", NA, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA, NA)
   if(traitType == "binary"){
     cctVec = c(cctVec, NA, NA)
   }
   cctVec = c(cctVec, NA, NA)

   pval.Region = rbind(pval.Region, cctVec)


   #print(pval.Region)
   #print(info.Region)
   ##remove columns with all NA
   pval.Region <- pval.Region[,which(unlist(lapply(pval.Region, function(x)!all(is.na(x))))),with=F]


   writeOutputFile(Output = list(pval.Region,  info.Region),
                    OutputFile = list(OutputFile, paste0(OutputFile, ".markerInfo")),
                    OutputFileIndex = OutputFileIndex,
                    AnalysisType = "Region",
                    nEachChunk = 1,
                    indexChunk = i,
                    Start = (i==1),
                    End = (i==nRegions))
   }


  }

  message = paste0("Analysis done! The results have been saved to '", OutputFile,"' and '",
                   paste0(OutputFile, ".markerInfo"),"'.")
  return(message)

}


# extract region-marker mapping from regionFile
SAIGE.getRegionList = function(groupFile,
			 annoVec, #c("lof","lof;missense" 
                         markerInfo)
{
  cat("Start extracting marker-level information from 'groupFile' of", groupFile, "....\n")

  Check_File_Exist(groupFile, "RegionFile")

  # read group file
  group_info_list =  ReadGroupFile(groupFile) #SAIGE_GENE_MultiVariantSet_Group.R
  # Main Test
  RegionData = NULL
  ngroup<-length(group_info_list)
  for(i in 1:ngroup){
    gene=group_info_list[[i]]$geneID	  
    var=group_info_list[[i]]$var
    anno=group_info_list[[i]]$anno
    RegionData = rbind(RegionData, cbind(rep(gene, length(var)), var, anno))
  }
  #print("RegionData")
  #print(RegionData)
  colnames(RegionData) = c("REGION", "SNP", "ANNO")
  RegionData = as.data.frame(RegionData)
  #print(RegionData)

##for non-VCF files  
if(!is.null(markerInfo)){
  # updated on 2021-08-05
  colnames(markerInfo)[3] = "MARKER"
  colnames(markerInfo)[6] = "genoIndex"
  colnames(markerInfo)[1] = "CHROM"
  RegionData = merge(RegionData, markerInfo, by.y = "MARKER", by.x = "SNP", all.x = T, sort = F)
  posNA = which(is.na(RegionData$genoIndex))
 # print("RegionData")
 # print(RegionData)
 # print("posNA")  
 # print(posNA)

  if(length(posNA) != 0){
    #print(head(RegionData[posNA,1:2]))
    stop("Total ",length(posNA)," markers in 'RegionFile' are not in 'GenoFile'.
         Please remove these markers before region-level analysis.")
  }
}

  #HeaderInRegionData = colnames(RegionData)
  HeaderInRegionData = unique(RegionData$ANNO)
  RegionAnnoHeaderList = list()
  if(length(annoVec) == 0){	
	stop("At least one annotation is required\n")
  }
  for(q in 1:length(annoVec)){
  	RegionAnnoHeaderList[[q]] = strsplit(annoVec[q],";")[[1]]
  }	  

  RegionList = list()
  uRegion = unique(RegionData$REGION)
  RegionData = as.data.frame(RegionData)

  for(r in uRegion){
    #print(paste0("Analyzing region ",r,"...."))
    #print(RegionData$REGION)
    #print(r)
    #which(as.vector(RegionData$REGION) == r)
    posSNP = which(RegionData$REGION == r)
    SNP = RegionData$SNP[posSNP]


    if(any(duplicated(SNP)))
      stop("Please check RegionFile: in region ", r,": duplicated SNPs exist.")

    if(!is.null(markerInfo)){
      genoIndex = as.numeric(RegionData$genoIndex[posSNP])
      chrom = RegionData$CHROM[posSNP]
      uchrom = unique(chrom)
      if(length(uchrom) != 1)
        stop("In region ",r,", markers are from multiple chromosomes.")
    }

    annoIndicatorMat = matrix(0, nrow=length(posSNP), ncol=length(annoVec)) 	 
    annoVecNew = c() 
  for(q in 1:length(annoVec)){
	indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]])
 	if(length(indiceVec) > 0){
	        annoVecNew = c(annoVecNew, annoVec[q])
		annoIndicatorMat[indiceVec, q] = 1			
		#annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
	}	
    }

  RegionAnnoHeaderListNew = list()
  if(length(annoVecNew) == 0){
	warning("No markers are found for at least one annotation, so this group will be skipped\n")
        #stop("At least one annotation is required\n")
  }else{
   if(length(annoVecNew) < length(annoVec)){
        annoIndicatorMat = matrix(0, nrow=length(posSNP), ncol=length(annoVecNew))	   
    for(q in 1:length(annoVecNew)){
        RegionAnnoHeaderListNew[[q]] = strsplit(annoVecNew[q],";")[[1]]
        indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderListNew[[q]])
        #if(length(indiceVec) > 0){
        annoIndicatorMat[indiceVec, q] = 1
                #annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
    }
  }else{
	annoVecNew = annoVec
  }	    


  }

   if(!is.null(markerInfo)){ 
    RegionList[[r]] = list(SNP = SNP,
			   annoIndicatorMat = annoIndicatorMat,
                           genoIndex = genoIndex,
#                           chrom = uchrom, 
			   annoVec = annoVecNew)
   }else{
    RegionList[[r]] = list(SNP = SNP,
			   annoIndicatorMat = annoIndicatorMat,
 #                          chrom = uchrom, 
			   annoVec = annoVecNew)

   }	   

  }


  return(RegionList)
}



##Working
mainRegion = function(genoType, genoIndex, annoIndicatorMat, maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, isImputation, annolist, isCondition, regionTestType)
{


  OutList = mainRegionInCPP(genoType, genoIndex, annoIndicatorMat, maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType)
  
 
 if(length(OutList) > 1){
  noNAIndices = which(!is.na(OutList$pvalVec))
  #print("dim(OutList$VarMat)")
  #print(dim(OutList$VarMat))

  URindices = which(OutList$infoVec == "UR")
  #print(length(which(OutList$infoVec == "")))
  OutList$infoVec[which(OutList$infoVec == "")] = "NA:NA:NA:NA"
  if(length(URindices) > 1){
  	OutList$infoVec[URindices] = "UR:UR:UR:UR"
  	#URindices_adj = URindices - length(genoIndex)
  }
  for(j in 1:length(annolist)){
        AnnoName = annolist[j]
        for(m in 1:length(maxMAFlist)){
                jm = (j-1)*(length(maxMAFlist)) + m
                maxMAFName = maxMAFlist[m]
		OutList$markerVec[URindices][jm] = paste0(AnnoName,":",maxMAFName) 
   	}
  }


  OutList$CHR = sapply(strsplit(as.character(OutList$infoVec), ":"), "[[", 1)
  OutList$POS = sapply(strsplit(as.character(OutList$infoVec), ":"), "[[", 2)
  OutList$Allele1 = sapply(strsplit(as.character(OutList$infoVec), ":"), "[[", 3)
  OutList$Allele2 = sapply(strsplit(as.character(OutList$infoVec), ":"), "[[", 4)

  obj.mainMarker = data.frame(CHR = OutList$CHR,
                               POS = OutList$POS,
                               MarkerID = OutList$markerVec,
                               Allele1 = OutList$Allele1,
                               Allele2 = OutList$Allele2,
                               AC_Allele2 = OutList$altCountsVec,
                               AF_Allele2 = OutList$altFreqVec)

  if(isImputation){
    obj.mainMarker$imputationInfo = OutList$imputationInfo
  }else{
    obj.mainMarker$MissingRate = OutList$missingRateVec
  }

  obj.mainMarker$BETA = OutList$BetaVec
  obj.mainMarker$SE = OutList$seBetaVec
  obj.mainMarker$Tstat = OutList$TstatVec
  obj.mainMarker$var = OutList$varTVec
  obj.mainMarker$p.value = OutList$pvalVec

  if(traitType == "binary"){
    obj.mainMarker$p.value.NA = OutList$pvalNAVec
    obj.mainMarker$Is.SPA.converge = OutList$isSPAConvergeVec
    obj.mainMarker$AF_caseVec = OutList$AF_caseVec
    obj.mainMarker$AF_ctrlVec = OutList$AF_ctrlVec
    obj.mainMarker$N_caseVec = OutList$N_caseVec
    obj.mainMarker$N_ctrlVec = OutList$N_ctrlVec

  }else if(traitType == "quantitative"){
    obj.mainMarker$N = obj.mainMarker$N_Vec
  }

    if(isCondition){
       obj.mainMarker$BETA_c = OutList$Beta_cVec
       obj.mainMarker$SE_c = OutList$seBeta_cVec
       obj.mainMarker$Tstat_c = OutList$Tstat_cVec
       obj.mainMarker$var_c = OutList$varT_cVec
       obj.mainMarker$p.value_c = OutList$pval_cVec
       if(traitType == "binary"){
                obj.mainMarker$p.value.NA_c = OutList$pvalNA_cVec
       }
    }


  obj.mainMarker = obj.mainMarker[which(obj.mainMarker$MarkerID != ""), ]
  if(!is.null(OutList$VarMatAdjCond)){
	OutList$VarMatAdjCond = OutList$VarMatAdjCond[noNAIndices,noNAIndices]
  	OutList$TstatAdjCond = OutList$TstatAdjCond[noNAIndices]
	OutList$G1tilde_P_G2tilde_Weighted_Mat = OutList$G1tilde_P_G2tilde_Weighted_Mat[noNAIndices,]
  }	  

  if(regionTestType == "BURDEN"){
    OutList$VarVec = OutList$varTVec
  }

  obj.mainMarker = obj.mainMarker[noNAIndices,]
  OutList$gyVec = OutList$gyVec[noNAIndices]
  return(list(outList = OutList,
              info.Region = obj.mainMarker))
 }else{
  return(list())
 }	 
}



mainRegionURV = function(NullModelClass = "SAIGE_NULL_Model",
                         genoType,
                         genoIndex,
                         n)
{
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainRegionURV = mainRegionURVInCPP("SAIGE", genoType, genoIndex, n)

  return(obj.mainRegionURV)
}

