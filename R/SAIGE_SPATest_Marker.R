mainMarker = function(genoType, genoIndex, traitType, isMoreOutput, isImputation, isCondition)
{
    # Check 'Main.cpp'
   #time_mainMarkerInCPP = system.time({OutList = mainMarkerInCPP(genoType, traitType, genoIndex, isMoreOutput)})
   OutList = mainMarkerInCPP(genoType, traitType, genoIndex, isMoreOutput)

   #print("time_mainMarkerInCPP")
   #print(time_mainMarkerInCPP)

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
    if(isCondition){
       obj.mainMarker$BETA_c = OutList$Beta_cVec
       obj.mainMarker$SE_c = OutList$seBeta_cVec
       obj.mainMarker$Tstat_c = OutList$Tstat_cVec
       obj.mainMarker$var_c = OutList$varT_cVec
       obj.mainMarker$p.value_c = OutList$pval_cVec
      obj.mainMarker$p.value.NA_c = OutList$pvalNA_cVec
    }	    
    obj.mainMarker$AF_caseVec = OutList$AF_caseVec
    obj.mainMarker$AF_ctrlVec = OutList$AF_ctrlVec
    obj.mainMarker$N_caseVec = OutList$N_caseVec
    obj.mainMarker$N_ctrlVec = OutList$N_ctrlVec

    if(isMoreOutput){
      obj.mainMarker$N_case_homVec = OutList$N_case_homVec
      obj.mainMarker$N_case_hetVec = OutList$N_case_hetVec
      obj.mainMarker$N_ctrl_homVec = OutList$N_ctrl_homVec
      obj.mainMarker$N_ctrl_hetVec = OutList$N_ctrl_hetVec
    }

  }else if(traitType == "quantitative"){
    if(isCondition){
       obj.mainMarker$BETA_c = OutList$Beta_cVec
       obj.mainMarker$SE_c = OutList$seBeta_cVec
       obj.mainMarker$Tstat_c = OutList$Tstat_cVec
       obj.mainMarker$var_c = OutList$varT_cVec
       obj.mainMarker$p.value_c = OutList$pval_cVec
    }	    
    obj.mainMarker$N = obj.mainMarker$N_Vec
  } 	  

  noNAIndices = which(!is.na(OutList$pvalVec))
  obj.mainMarker = obj.mainMarker[noNAIndices,]
  rm(OutList)
  return(obj.mainMarker)
}



setMarker = function(objNull, chrom, impute_method, missing_cutoff, min_maf_marker, min_mac_marker, min_info_marker, omp_num_threads, isOutputMoreDetails)
{
  # Check Main.cpp
  setMarker_GlobalVarsInCPP(impute_method,
                            missing_cutoff,
                            min_maf_marker,
                            min_mac_marker,
			    min_info_marker,
                            omp_num_threads,
			    isOutputMoreDetails
                            )

  # Check SAIGE.R
  obj.setMarker = setMarker.SAIGE(objNull, control)

  return(obj.setMarker)
}




SAIGE.Marker = function(objNull,
			objGeno,
                        OutputFile,
                        OutputFileIndex = NULL,
                        nMarkersEachChunk, 
			isMoreOutput,
			isImputation,
			LOCO,
			chrom,
			isCondition)
{

  if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index")
  
  genoType = objGeno$genoType

  outIndex = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F))    # this function is in 'Util.R'
  outIndex = outIndex$indexChunk
  if(outIndex != 1)
    cat("Restart the analysis from chunk:\t", outIndex, "\n")


  ## set up an object for genotype
  if(genoType != "vcf"){
    markerInfo = objGeno$markerInfo
    if(LOCO){
      markerInfo = markerInfo[which(markerInfo$CHROM == chrom),]  
    }
    CHROM = markerInfo$CHROM
    genoIndex = markerInfo$genoIndex
    ##only for one chrom
    # all markers were split into multiple chunks,
    genoIndexList = splitMarker(genoIndex, nMarkersEachChunk, CHROM);
    nChunks = length(genoIndexList)
    cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
    cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
    cat("Number of chunks for all markers:\t", nChunks, "\n")
    if(outIndex > nChunks){
      stop("The analysis has been finished! Please delete ", OutputFileIndex, " if the analysis needs to be run again")
      is_marker_test = FALSE 
    }else{
      is_marker_test = TRUE
      i = outIndex    
    }
    
  }else{
    if(outIndex > 1){
	move_forward_iterator_Vcf(outIndex*nMarkersEachChunk)    
    }
    isVcfEnd =  check_Vcf_end()
    if(!isVcfEnd){
    	outIndex = 1
    	genoIndex = rep(-1, nMarkersEachChunk) 
	nChunks = outIndex + 1
	is_marker_test = TRUE
    }else{
	is_marker_test = FALSE    
	stop("No markers are left in VCF")
    }
  }

  chrom = "InitialChunk"

  while(is_marker_test){
  #for(i in outIndex:nChunks)
  #{
#time_left = system.time({
    if(genoType != "vcf"){	  
      tempList = genoIndexList[[i]]
      genoIndex = as.integer(tempList$genoIndex)
      tempChrom = tempList$chrom
    }

    #print("tempList")
    #print(tempList)
    #print(tempList$genoIndex)
#})
#print("time_left")
#print(time_left)
    #print("genoIndex here")
    #print(genoIndex)
    # set up objects that do not change for different variants
    #if(tempChrom != chrom){
    #  setMarker("SAIGE", objNull, control, chrom, Group, ifOutGroup)
    #  chrom = tempChrom
    #}
    #ptm <- proc.time()
    #print(ptm)
    #print("gc()")
    #print(gc())
    cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \n"))

    # main function to calculate summary statistics for markers in one chunk
    #time_mainMarker = system.time({resMarker = mainMarker(genoType, genoIndex, objNull$traitType, isMoreOutput, isImputation, isCondition)})
    resMarker = as.data.frame(mainMarkerInCPP(genoType, objNull$traitType, genoIndex, isMoreOutput, isImputation)) 
    resMarker = resMarker[which(!is.na(resMarker$BETA)), ]

#    print("time_mainMarker")
#   print(time_mainMarker)    


    #timeoutput=system.time({writeOutputFile(Output = list(resMarker),
  if(nrow(resMarker) > 0){
  writeOutputFile(Output = list(resMarker),
                    OutputFile = list(OutputFile),
                    OutputFileIndex = OutputFileIndex,
                    AnalysisType = "Marker",
                    nEachChunk = format(nMarkersEachChunk, scientific=F),
                    indexChunk = i,
                    Start = (i==1),
                    End = (i==nChunks))

  }
                    #End = (i==nChunks))})
    #print("timeoutput")
    #print(timeoutput)
    ptm <- proc.time()
    print(ptm)
    print("gc()")
    print(gc())
    #rm(resMarker)


  if(genoType == "vcf"){
    isVcfEnd =  check_Vcf_end()
    if(isVcfEnd){
	is_marker_test = FALSE	     
    }
  }else{
    i = i + 1
    if(i > nChunks){
      is_marker_test = FALSE
    }	    
  }
	  
  } #while(is_marker_test){

  # information to users
  output = paste0("Analysis done! The results have been saved to '", OutputFile,"'.")

  return(output)
}


splitMarker = function(genoIndex, nMarkersEachChunk, CHROM)
{
  genoIndexList = list()
  iTot = 1;

  uCHROM = unique(CHROM)
  for(chrom in uCHROM){
    pos = which(CHROM == chrom)
    gIdx = genoIndex[pos]
    M = length(gIdx)

    idxStart = seq(1, M, nMarkersEachChunk)
    idxEnd = idxStart + nMarkersEachChunk - 1

    nChunks = length(idxStart)
    idxEnd[nChunks] = M

    for(i in 1:nChunks){
      idxMarker = idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] = list(chrom = chrom,
                                   genoIndex = gIdx[idxMarker])
      iTot = iTot + 1;
    }
  }

  return(genoIndexList)
}
