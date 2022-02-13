getChromNumber = function(chrom = ""){
  if(chrom == ""){
    stop("chrom is not specified\n")
  }else{	  
    chrom_v2 = as.character(chrom)
    chrom_v2 = gsub("CHR", "", chrom_v2, ignore.case = T)
    chrom_v3 = as.numeric(gsub("[^0-9.]", "",chrom_v2))
    if(chrom_v3 > 22 | chrom_v3 < 1) {
      stop("chromosome ", chrom, " is out of the range of null model LOCO results\n")
    }else {
      cat("Leave chromosome ", chrom_v3, " out will be applied\n")
    }
  }
  return(chrom_v3)
}	

removeLOCOResult = function(chromList, obj.glmm.null){
  for (chr in 1:22) {
    if(chr %in% chromList){    
      obj.glmm.null$LOCOResult[chr] = list(NULL)
      cat("chromosome ", chr, " model results are removed to save memory\n")
      gc()
    }
  }
  return(obj.glmm.null)
}	



ReadModel = function(GMMATmodelFile = "", chrom="", LOCO=TRUE, is_Firth_beta=FALSE){	
  # Check file existence
  Check_File_Exist(GMMATmodelFile, "GMMATmodelFile")
  if(!LOCO %in% c(TRUE, FALSE))
    stop("LOCO should be TRUE or FALSE.")
  # load GMMATmodelFile
  load(GMMATmodelFile)
  if(!modglmm$LOCO){
	if(LOCO){
		LOCO=FALSE
	}	
  }
	  
  obj.glmm.null = modglmm
  obj.glmm.null$Y = NULL
  #obj.glmm.null$offset = obj.glmm.null$linear.predictors - obj.glmm.null$coefficients[1]
  obj.glmm.null$linear.predictors = NULL
  obj.glmm.null$coefficients = NULL
  obj.glmm.null$cov = NULL
  rm(modglmm)
  gc()
  #traitType = obj.glmm.null$traitType
  #y = obj.glmm.null$y
  #X = obj.glmm.null$X
  #N = length(y)
  #tauVec = obj.glmm.null$theta
  #indChromCheck = FALSE
  chrom_v3=NULL

  if(!LOCO) {
    print("Leave-one-chromosome-out is not applied")
    if(obj.glmm.null$LOCO) {
      for (chr in 1:22) {
        obj.glmm.null$LOCOResult[chr] = list(NULL)
        cat("chromosome ", chr, " model results are removed to save memory\n")
        gc()
      }
    }
  }else{
    if (!obj.glmm.null$LOCO){
      stop("LOCO is TRUE but the null model file .rda does not contain LOCO results. In order to apply Leave-one-chromosome-out, please run Step 1 using LOCO. Otherwise, please set LOCO=FALSE in this step (Step 2).\n")
    }else{
        if(chrom == ""){
          stop("chrom needs to be specified in order to apply Leave-one-chromosome-out on gene- or region-based tests")
        }else{
          chrom_v3 = getChromNumber(chrom)
        }
   }

 if(!is.null(chrom_v3)){
   chromList = c(1:22)
   chromList = chromList[which(chromList != chrom_v3)]   
   obj.glmm.null = removeLOCOResult(chromList, obj.glmm.null)
   obj.glmm.null$fitted.values = obj.glmm.null$LOCOResult[[chrom_v3]]$fitted.values
   obj.glmm.null$residuals = obj.glmm.null$LOCOResult[[chrom_v3]]$residuals
   #obj.glmm.null$offset = obj.glmm.null$LOCOResult[[chrom_v3]]$linear.predictors -  obj.glmm.null$LOCOResult[[chrom_v3]]$coefficients[1]
   obj.glmm.null$obj.noK = obj.glmm.null$LOCOResult[[chrom_v3]]$obj.noK
   if(is_Firth_beta){
     if(!is.null(obj.glmm.null$LOCOResult[[chrom_v3]]$offset)){	   
       obj.glmm.null$offset = obj.glmm.null$LOCOResult[[chrom_v3]]$offset
     }
   }
   obj.glmm.null$LOCOResult[chrom_v3] = list(NULL)
   gc() 
 }
  }

 obj.glmm.null$mu = as.vector(obj.glmm.null$fitted.values)
 tau = obj.glmm.null$theta
 N = length(obj.glmm.null$mu)
     if(obj.glmm.null$traitType == "binary"){
             obj.glmm.null$mu2 = (obj.glmm.null$mu) *(1-obj.glmm.null$mu)
           }else if(obj.glmm.null$traitType == "quantitative"){
             obj.glmm.null$mu2 = (1/tau[1])*rep(1,N)
           }
 #if(FALSE){
 if(is_Firth_beta){
	if(obj.glmm.null$traitType == "binary"){
		if(is.null(obj.glmm.null$offset)){
			#if(FALSE){
			print("WARNING. is_Firth_beta = TRUE, but offset was not computed in Step 1. Please re-run Step 1 using the more rencet version of SAIGE/SAIGE-GENE.")
			cat("Applying is_Firth_beta = TRUE and note the results would be more accurate withe Step 1 results using the more rencet version of SAIGE/SAIGE-GENE. \n")

			if(ncol(obj.glmm.null$X) == 1){
				covoffset = rep(0,nrow(obj.glmm.null$X))
			}else{	
				formulastr = paste0("y ~ ", paste(colnames(obj.glmm.null$X)[-1], collapse="+"))
				obj.glmm.null$X = cbind(obj.glmm.null$X, as.vector(obj.glmm.null$y))
				colnames(obj.glmm.null$X)[ncol(obj.glmm.null$X)] = "y"
				obj.glmm.null$X = as.data.frame(obj.glmm.null$X)
				formula.new = as.formula(formulastr)
				modwitcov = glm(formula.new, data = obj.glmm.null$X, family = binomial)
				obj.glmm.null$X = obj.glmm.null$X[,-ncol(obj.glmm.null$X)]
				covoffset = as.matrix(obj.glmm.null$X[,-1]) %*%  modwitcov$coefficients[-1]
				obj.glmm.null$X = as.matrix(obj.glmm.null$X)
			}
			obj.glmm.null$offset = covoffset
			#}
		}
	}	
  }
 #}	 
 return(obj.glmm.null)    
}


Get_Variance_Ratio<-function(varianceRatioFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest){

    ratioVec = c(1)
    # check variance ratio
    if (!file.exists(varianceRatioFile)) {
        #if (sparseSigmaFile == "") {
        #    stop("ERROR! varianceRatioFile ", varianceRatioFile, " does not exsit but sparseSigmaFile also does not exist \n")
        #}else {
        cat("varianceRatioFile is not specified so variance ratio won't be used\n")
        #}
        ratioVec = c(1)
    }else {
        varRatioData = data.frame(data.table:::fread(varianceRatioFile, header = F, stringsAsFactors = FALSE))
        ln = length(cateVarRatioMinMACVecExclude)
        hn = length(cateVarRatioMaxMACVecInclude)
	if(isGroupTest){
        if (nrow(varRatioData) == 1) {
            stop("ERROR! To perform gene-based tests, categorical variance ratios are required\n")
        }else {
            ratioVec = varRatioData[, 1]
            nrv = length(ratioVec)
            if (nrv != ln) {
                stop("ERROR! The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude\n")
            }
            if (ln != (hn + 1)) {
                stop("ERROR! The length of cateVarRatioMaxMACVecInclude does not match with the lenght of cateVarRatioMinMACVecExclude (-1)\n")
            }
        }
	
	}else{
		ratioVec = varRatioData[, 1]

	 }	
        cat("variance Ratio is ", ratioVec, "\n")
    }
    return(ratioVec)
}


