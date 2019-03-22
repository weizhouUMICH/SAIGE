outputDiagGRM = function(
		GRMdiagfile = "",
		GRMdiagSamplefile = "",
		plinkFile = "", 
                phenoFile = "",
                phenoCol = "",
                traitType = "binary",
                invNormalize = FALSE,
                covarColList = NULL,
                qCovarCol = NULL,
                sampleIDColinphenoFile = "",
                tol=0.02,
                maxiter=20,
                tolPCG=1e-5,
                maxiterPCG=500,
                nThreads = 1, 
                SPAcutoff = 2, 
                numMarkers = 30, 
                skipModelFitting = FALSE,
		memoryChunk = 2,
		tauInit = c(0,0),
		LOCO = FALSE,
		traceCVcutoff = 0.0025,
		ratioCVcutoff = 0.001, 
                outputPrefix = "",
		outputPrefix_varRatio = NULL,
		IsSparseKin = FALSE,
		sparseGRMFile=NULL,
                sparseGRMSampleIDFile=NULL,
		numRandomMarkerforSparseKin = 1000,
		relatednessCutoff = 0.125, 
		isCateVarianceRatio = FALSE,
		cateVarRatioIndexVec = NULL,
		cateVarRatioMinMACVecExclude = c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5),
		cateVarRatioMaxMACVecInclude = c(1.5,2.5,3.5,4.5,5.5,10.5,20.5),
		isCovariateTransform = TRUE,
		isDiagofKinSetAsOne = FALSE,
		useSparseSigmaConditionerforPCG = FALSE,
		useSparseSigmaforInitTau = FALSE){


  if(useSparseSigmaConditionerforPCG){
    cat("sparse sigma will be used as the conditioner for PCG\n")
    if(!file.exists(sparseGRMFile)){
      stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
    }
  }


  if(useSparseSigmaforInitTau){
    cat("sparse sigma will be used to estimate the inital tau\n")
    if(!file.exists(sparseGRMFile)){
      stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
    }
  }

  
  if(nThreads > 1){
    RcppParallel:::setThreadOptions(numThreads = nThreads)
    cat(nThreads, " threads are set to be used ", "\n")
  }

  #check and read files
  #output file
  modelOut=paste0(outputPrefix, ".rda")
  SPAGMMATOut=paste0(outputPrefix, "_", numMarkers,"markers.SAIGE.results.txt")

  if (is.null(outputPrefix_varRatio)){
    outputPrefix_varRatio = outputPrefix
  }

  varRatioFile=paste0(outputPrefix_varRatio,".varianceRatio.txt")

  if(!file.exists(varRatioFile)){
    file.create(varRatioFile, showWarnings = TRUE)
  }else{
    stop("WARNING: The variance ratio file ", varRatioFile, " already exists. The new variance ratios will be output to ", varRatioFile,". In order to avoid over-writting, please remove the ", varRatioFile, " or use the argument outputPrefix_varRatio to specify a different prefix to output the variance ratio(s)\n")
  }


  if(!file.exists(modelOut)){
    file.create(modelOut, showWarnings = TRUE)
  }

  if(!file.exists(paste0(plinkFile, ".bed"))){
    stop("ERROR! ", plinkFile, ".bed does not exsit\n")
  }

  if(!file.exists(paste0(plinkFile, ".bim"))){
    stop("ERROR! ", plinkFile, ".bim does not exsit\n")
  }else{
      chromosomeStartIndexVec = NULL
      chromosomeEndIndexVec = NULL
    ###if LOCO, record the indices of markers on each chromosome
    if(LOCO){
      cat("leave-one-chromosome-out is activated! Note this option will only be applied to autosomal variants\n")
   
      bimData = data.table:::fread(paste0(plinkFile,".bim"),  header=F)
      for(i in 1:22){
	if(length(which(bimData[,1] == i)) > 0){
          chromosomeStartIndexVec = c(chromosomeStartIndexVec, min(which(bimData[,1] == i))-1)
	  chromosomeEndIndexVec = c(chromosomeEndIndexVec, max(which(bimData[,1] == i))-1)
	  if(chromosomeStartIndexVec[i] <= chromosomeStartIndexVec[i-1] | chromosomeEndIndexVec[i] <= chromosomeEndIndexVec[i-1]){
		stop(paste0("ERROR! chromosomes need to be ordered from 1 to 22 in ", plinkFile, ".bim\n"))
	  }

	}else{
	  chromosomeStartIndexVec = c(chromosomeStartIndexVec, NA)
	  chromosomeEndIndexVec = c(chromosomeEndIndexVec, NA)

        }   	
      }
      cat("chromosomeStartIndexVec: ", chromosomeStartIndexVec, "\n")
      cat("chromosomeEndIndexVec: ", chromosomeEndIndexVec, "\n")

     # setChromosomeIndicesforLOCO(chromosomeStartIndexVec, chromosomeEndIndexVec, chromosomeVecVec) 
    }else{
      chromosomeStartIndexVec = rep(NA, 22)
      chromosomeEndIndexVec = rep(NA, 22)
    }	
  }

  if(!file.exists(paste0(plinkFile, ".fam"))){
    stop("ERROR! ", plinkFile, ".fam does not exsit\n")
  }else{
    sampleListwithGenov0 = data.table:::fread(paste0(plinkFile,".fam"),  header=F)
    sampleListwithGenov0 = data.frame(sampleListwithGenov0)
    colnames(sampleListwithGenov0) = c("FIDgeno", "IIDgeno", "father", "mother", "sex", "phe")
    sampleListwithGeno = NULL
    sampleListwithGeno$IIDgeno = sampleListwithGenov0$IIDgeno
    sampleListwithGeno = data.frame(sampleListwithGeno)
    sampleListwithGeno$IndexGeno = seq(1,nrow(sampleListwithGeno), by=1)
    cat(nrow(sampleListwithGeno), " samples have genotypes\n")
  }


  #phentoype file
  if(!file.exists(phenoFile)){
    stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
  }else{
    #ydat = data.table:::fread(phenoFile, header=T, stringsAsFactors=FALSE)
    if( grepl(".gz$",phenoFile) | grepl(".bgz$",phenoFile) ) {
      ydat = data.table:::fread(cmd=paste0("gunzip -c ", phenoFile), header=T, stringsAsFactors=FALSE)
    } else {
      ydat = data.table:::fread(phenoFile, header=T, stringsAsFactors=FALSE)
    }

    data = data.frame(ydat)

    for(i in c(phenoCol, covarColList, qCovarCol, sampleIDColinphenoFile)){
      if(!(i %in% colnames(data))){
        stop("ERROR! column for ", i, " does not exsit in the phenoFile \n")
      }
    }

    #update the categorical variables

    if(length(covarColList) > 0){
      #qCovarColUpdated = NULL
      #for(i in qCovarCol){
      #  j = paste0("factor(", i, ")")
      #  qCovarColUpdated = c(qCovarColUpdated, j)
      #}
      formula = paste0(phenoCol,"~", paste0(covarColList,collapse="+"))
      #formula = paste0(phenoCol,"~",paste0(c(covarColList,qCovarColUpdated),collapse="+"))
      hasCovariate = TRUE
    }else{
      formula = paste0(phenoCol,"~ 1")
      hasCovariate = FALSE
    }    

    cat("formula is ", formula,"\n")
    formula.null = as.formula(formula)
    mmat = model.frame(formula.null, data, na.action=NULL)
    mmat$IID = data[,which(sampleIDColinphenoFile == colnames(data))]
    mmat_nomissing = mmat[complete.cases(mmat),]
    mmat_nomissing$IndexPheno = seq(1,nrow(mmat_nomissing), by=1)
    cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")

    dataMerge = merge(mmat_nomissing, sampleListwithGeno, by.x = "IID", by.y = "IIDgeno")
    dataMerge_sort = dataMerge[with(dataMerge, order(IndexGeno)), ]

    if(nrow(dataMerge_sort) < nrow(sampleListwithGeno)){
      cat(nrow(sampleListwithGeno) - nrow(dataMerge_sort), " samples in geno file do not have phenotypes\n")
    }
    cat(nrow(dataMerge_sort), " samples will be used for analysis\n")
#    cat("dataMerge_sort$IID ", dataMerge_sort$IID, "\n")
  }

  if(invNormalize){
      cat("Perform the inverse nomalization for ", phenoCol, "\n")
      invPheno = qnorm((rank(dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)], na.last="keep")-0.5)/sum(!is.na(dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)])))
      dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)] = invPheno
  }


  #check for perfect separation
  if(traitType == "binary"){
    out_checksep = checkPerfectSep(formula.null, data=dataMerge_sort)
    covarColList <- covarColList[!(covarColList %in% out_checksep)]
    formula = paste0(phenoCol,"~", paste0(covarColList,collapse="+"))
    formula.null = as.formula(formula)
    dataMerge_sort <- dataMerge_sort[, !(names(dataMerge_sort) %in% out_checksep)]
  }


  if(isCovariateTransform){
    cat("qr transformation has been performed on covariates\n")
    out.transform<-Covariate_Transform(formula.null, data=dataMerge_sort)
    formulaNewList = c("Y ~ ", out.transform$Param.transform$X_name[1])
    if(length(out.transform$Param.transform$X_name) > 1){
      for(i in c(2:length(out.transform$Param.transform$X_name))){
        formulaNewList = c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
      }
    }

    formulaNewList = paste0(formulaNewList, collapse="")
    formulaNewList = paste0(formulaNewList, "-1")
    formula.new = as.formula(paste0(formulaNewList, collapse=""))
    data.new = data.frame(cbind(out.transform$Y, out.transform$X1))
    colnames(data.new) = c("Y",out.transform$Param.transform$X_name)
    cat("colnames(data.new) is ", colnames(data.new), "\n")
    cat("out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr), "\n")

  }else{ #if(isCovariateTransform) else
    formula.new = formula.null
    data.new = dataMerge_sort
  }

  setgeno(plinkFile, dataMerge_sort$IndexGeno, memoryChunk, isDiagofKinSetAsOne)
  diagGRM = get_GRMdiagVec()
  if(GRMdiagfile != "" & GRMdiagSamplefile != ""  ){
  write.table(diagGRM, GRMdiagfile, col.names=F, row.names=F, quote=F)
  print(colnames(dataMerge_sort))
  write.table(dataMerge_sort$IID, GRMdiagSamplefile, col.names=F, row.names=F, quote=F) 
  data2 = cbind(data.new,dataMerge_sort$IID)
  colnames(data2)[length(colnames(data2))] = "IID"  
  write.table(data2, paste0(GRMdiagSamplefile, ".data.new.txt"), col.names=T, row.names=F, quote=F) 
  print(formula.new)
  write.table(as.character(formula.new), paste0(GRMdiagSamplefile, ".formula.new.txt"), col.names=F, row.names=F, quote=F) 

  }else{
    stop("GRMdiagfile or GRMdiagSamplefile is not specified\n")
  }

}

