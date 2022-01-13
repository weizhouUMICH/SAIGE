SPAGMMATtest = function(bgenFile = "",
                 bgenFileIndex = "",
                 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
                 savFileIndex = "",
                 sampleFile = "",
                 bedFile="",
                 bimFile="",
                 famFile="",
                 AlleleOrder = "alt-first", #new
                 idstoExcludeFile = NULL,
                 idstoIncludeFile = NULL,
                 rangestoExcludeFile = NULL,
                 rangestoIncludeFile = NULL,
                 chrom = "",
                 start = 1,
                 end = 250000000,
                 max_missing = 0.15,  #new
		 impute_method = "mean",  #"drop", "mean", "minor"     #new
                 min_MAC = 0.5,
                 min_MAF = 0,
                 min_Info = 0,
		 is_imputed_data = FALSE, #new
                 GMMATmodelFile = "",
                 LOCO=TRUE,
                 varianceRatioFile = "",
                 cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5),
                 cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5),

		 SPAcutoff=2,
                 SAIGEOutputFile = "",
                 numLinesOutput = 10,
                 condition="",
                 sparseSigmaFile="",
                 groupFile="",
                 kernel="linear.weighted",
                 method="optimal.adj",
                 weights.beta.rare = c(1,25),
                 weights.beta.common = c(1,25),
                 weightMAFcutoff = 0.01,
                 weightsIncludeinGroupFile=FALSE,
                 weights_for_G2_cond = NULL,
                 r.corr=0,

		 dosage_zerod_cutoff = 0.2,
		 dosage_zerod_MAC_cutoff = 10,
		 is_output_moreDetails = FALSE, #new
                 X_PARregion="60001-2699520,154931044-155270560",
                 is_rewrite_XnonPAR_forMales=FALSE,
                 sampleFile_male="",
                 method_to_CollapseUltraRare="absence_or_presence",
		 MACCutoff_to_CollapseUltraRare = 10,
                 DosageCutoff_for_UltraRarePresence = 0.5,
                 function_group_test =c("lof", "missense", "synonymous"),  #new
                 maxMAFforGroupTest = c(0.1),
                 max_markers_region = 100   #new
		 ){

   if(!(impute_method %in% c("mean","minor","drop"))){
     stop("impute_method should be 'mean', 'minor', or 'drop'.")
   }


   checkArgsListBool(is_imputed_data = is_imputed_data,
                     LOCO = LOCO,
		     is_output_moreDetails = is_output_moreDetails,
		     is_rewrite_XnonPAR_forMales = is_rewrite_XnonPAR_forMales)



   checkArgsListNumeric(start = start,
                     end = end,
		     max_missing = max_missing,
                     min_MAC = min_MAC,
                     min_MAF = min_MAF,
                     min_Info = min_Info,
                     SPAcutoff = SPAcutoff,
                     numLinesOutput = numLinesOutput,
                     dosage_zerod_cutoff = dosage_zerod_cutoff,
		     dosage_zerod_MAC_cutoff = dosage_zerod_MAC_cutoff)
    if(file.exists(SAIGEOutputFile)) {print("ok -2 file exist")} 



    ##check and create the output file
    #Check_OutputFile_Create(SAIGEOutputFile)
    OutputFile = SAIGEOutputFile
    OutputFileIndex=NULL
    if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index") 

    ##check the variance ratio file and extract the variance ratio vector

    if(groupFile == ""){
      isGroupTest = FALSE
      cat("single-variant association test will be performed\n")
	  setMarker_GlobalVarsInCPP(impute_method,
				    max_missing,
                            min_MAF,
                            min_MAC,
                            min_Info,
                            1,
			    is_output_moreDetails,
			    numLinesOutput,
			    dosage_zerod_cutoff,
			    dosage_zerod_MAC_cutoff
                            )

    }else{
      isGroupTest = TRUE
      Check_File_Exist(groupFile, "groupFile")
      cat("group-based test will be performed\n")
      checkArgsList_for_Region(method_to_CollapseUltraRare,
                                    MACCutoff_to_CollapseUltraRare,
                                    DosageCutoff_for_UltraRarePresence,
                                    maxMAFforGroupTest = maxMAFforGroupTest,
				    max_markers_region = max_markers_region)




    if(file.exists(SAIGEOutputFile)) {print("ok -1 file exist")} 

        IsOutputlogPforSingle = FALSE   #to check
        OUT_Filename_Single<-sprintf("%s.single",SAIGEOutputFile)
        Check_OutputFile_Create(OUT_Filename_Single)
      if (sum(weights.beta.rare != weights.beta.common) > 0) {
        cat("WARNING:The option for weights.beta.common is not fully developed\n")
        cat("weights.beta.common is set to be equal to weights.beta.rare\n")
        weights.beta.common = weights.beta.rare
      }

      setRegion_GlobalVarsInCPP(impute_method,
                                max_missing,
				maxMAFforGroupTest,
				max_markers_region,
				1,
				method_to_CollapseUltraRare,
				MACCutoff_to_CollapseUltraRare,
				DosageCutoff_for_UltraRarePresence,
				dosage_zerod_cutoff,
                            	dosage_zerod_MAC_cutoff
                            )

    }
    
    ratioVec = Get_Variance_Ratio(varianceRatioFile, sparseSigmaFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest) #readInGLMM.R
  
    sparseSigmaRList = list()
    isSparseGRM = TRUE
    if(sparseSigmaFile != ""){ 
      sparseSigmaRList = setSparseSigma(sparseSigmaFile)
      isSparseGRM = TRUE
    }else{
      sparseSigmaRList = list(nSubj = 0, locations = matrix(0,nrow=2,ncol=2), values = rep(0,2))  
      isSparseGRM = FALSE 
    }	    

    obj.model = ReadModel(GMMATmodelFile, chrom, LOCO) #readInGLMM.R


    nsample = length(obj.model$y)
    cateVarRatioMaxMACVecInclude = c(cateVarRatioMaxMACVecInclude, nsample)	
    #print(names(obj.model$obj.noK))


     #in Geno.R
    objGeno = setGenoInput(bgenFile = bgenFile,
                 bgenFileIndex = bgenFileIndex,
                 vcfFile = vcfFile,   #not activate yet
                 vcfFileIndex = vcfFileIndex,
                 vcfField = vcfField,
                 savFile = savFile,
                 savFileIndex = savFileIndex,
                 sampleFile = sampleFile,
                 bedFile=bedFile,
                 bimFile=bimFile,
                 famFile=famFile,
                 idstoExcludeFile = idstoExcludeFile,
                 idstoIncludeFile = idstoIncludeFile,
                 rangestoExcludeFile = rangestoExcludeFile,
                 rangestoIncludeFile = rangestoIncludeFile,
                 chrom = chrom,
                 start = start,
                 end = end,
                 AlleleOrder = AlleleOrder,
                 sampleInModel = obj.model$sampleID)

    markerInfo = objGeno$markerInfo
    genoIndex = markerInfo$genoIndex
    #print("genoIndex")
    #print(genoIndex[1:10])
    genoType = objGeno$dosageFileType


    cat("isSparseGRM ", isSparseGRM, "\n")

   if (condition != "") {
        isCondition = TRUE
        #n = length(obj.model$y) #sample size
        #print(n)
        #assign_conditionMarkers_factors(genoType, condition_genoIndex,  n)
    }
    else {
        isCondition = FALSE
    }
    
    condition_genoIndex = c(-1)
    cat("isCondition ", isCondition, "\n") 
    #set up the SAIGE object based on the null model results
    setSAIGEobjInCPP(t_XVX=obj.model$obj.noK$XVX,
		     t_XXVX_inv=obj.model$obj.noK$XXVX_inv,
		     t_XV=obj.model$obj.noK$XV,
		     t_XVX_inv_XV=obj.model$obj.noK$XVX_inv_XV,
		     t_X=obj.model$X,
		     t_S_a=obj.model$obj.noK$S_a,
		     t_res=obj.model$residuals,
		     t_mu2=obj.model$mu2,
		     t_mu=obj.model$mu,
		     t_varRatio = as.vector(ratioVec),
		     t_cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
		     t_cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
		     t_SPA_Cutoff = SPAcutoff,
		     t_tauvec = obj.model$theta,
		     t_traitType = obj.model$traitType,
		     t_y = obj.model$y,
		     t_impute_method = impute_method, 
		     t_flagSparseGRM = isSparseGRM,
        	     t_locationMat = as.matrix(sparseSigmaRList$locations),
        	     t_valueVec = sparseSigmaRList$values,
        	     t_dimNum = sparseSigmaRList$nSubj, 
		     t_isCondition = isCondition,
		     t_condition_genoIndex = condition_genoIndex)


   #process condition
    if (isCondition) {
	 #print("OK1")
        n = length(obj.model$y) #sample size
        condition_genoIndex = extract_genoIndex_condition(condition, markerInfo)
	assign_conditionMarkers_factors(genoType, condition_genoIndex,  n)
	 #print("OK2")
	if(obj.model$traitType == "binary" & isGroupTest){
		outG2cond = RegionSetUpConditional_binary_InCPP()
	 #print("OK3")


	G2condList = get_newPhi_scaleFactor(q.sum = outG2cond$qsum_G2_cond, mu.a = obj.model$mu, g.sum = outG2cond$gsum_G2_cond, p.new = outG2cond$pval_G2_cond, Score = outG2cond$Score_G2_cond, Phi = outG2cond$VarMat_G2_cond)
	#print(G2condList)
	scaleFactorVec = as.vector(G2condList$scaleFactor)
	#print(scaleFactorVec)
		assign_conditionMarkers_factors_binary_region(scaleFactorVec)
	 #print("OK5")
	}	
    }


    if(file.exists(SAIGEOutputFile)) {print("ok 0 file exist")} 


    #cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
    #cat("Number of markers in each chunk:\t", numLinesOutput, "\n")
    #cat("Number of chunks for all markers:\t", nChunks, "\n")
    #}

    if(!isGroupTest){
    OutputFile = SAIGEOutputFile

    nMarkersEachChunk = numLinesOutput
    if(file.exists(SAIGEOutputFile)) {print("ok 2 file exist")}
        SAIGE.Marker(obj.model,
                   objGeno,
                   OutputFile,
                   OutputFileIndex,
                   nMarkersEachChunk,
                   is_output_moreDetails,
		   is_imputed_data,
                   LOCO,
                   chrom,
		   isCondition)
    }else{
	SAIGE.Region(obj.model,
		     objGeno,
		     sparseSigma,
		     OutputFile,
		     method_to_CollapseUltraRare,
		     MACCutoff_to_CollapseUltraRare,
                     DosageCutoff_for_UltraRarePresence,
                     groupFile,
                     function_group_test,
                     maxMAFforGroupTest,
                     max_markers_region,
                     genoType,
                     markerInfo,
		     obj.model$traitType,
		     is_imputed_data,
		     isCondition)
    }	    
}
