#' Fit the null logistic mixed model and estimate the variance ratio by a set of randomly selected variants 
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM)
#' @param phenoFile character. Path to the phenotype file
#' @param phenoCol character. Column name for the trait e.g. "CAD"
#' @param traitType character. e.g. "binary" or "quantitative". By default, "binary"
#' @param invNormalize logical. Whether to perform the inverse normalization of the trait or not. E.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the glm model e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the glm model (NOT work yet)
#' @param sampleIDColinphenoFile character.  Column name for the sample IDs in the phenotype file e.g. "IID".  
#' @param nThreads integer. Number of threads to be used. By default, 1 
#' @param numMarkers integer (>0). Number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed 
#' @param tauInit vector of numbers. e.g. c(1,1), Unitial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is not specified, the second element will be 0.5 for binary traits.  
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 4
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. 
#' @param traceCVcutoff float. The threshold for coefficient of variantion (CV) for the trace estimator to increase nrun
#' @param outputPrefix character. Path to the output files with prefix.
#' @param isCateVarianceRatio logical. Whether to estimate variance ratio based on different MAC categories. If yes, six categories will be used MAC = 1, 2, 3, 4, 5, >5. Currently, if isCateVarianceRatio=TRUE, then LOCO=FALSE 
#' @param IsSparseKin logical. Whether to exploit the sparsity of GRM to estimate the variance ratio. By default, TRUE
#' @param numRandomMarkerforSparseKin integer (>0). Number of markers to be used for first estimating the relatedness between each sample pair  if IsSparseKin is TRUE
#' @param relatednessCutoff float. The threshold to treat two samples as unrelated if IsSparseKin is TRUE
#' @param methodforRelatedSample character. The method to fit model for related samples. GMMAT or EMMAX
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio value, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fit_SKAT_NULL = function(kins = NULL, 
                phenoFile = "",
                phenoCol = "",
                traitType = "quantitative",
                invNormalize = FALSE,
                covarColList = NULL,
                qCovarCol = NULL,
                sampleIDColinphenoFile = "",
		outputPrefix = "",
		isCovariateTransform = FALSE,
		sampleFileForDosages="",
		methodforRelatedSample="EMMAX",
		isDiagofKinSetAsOne = FALSE){

  #check and read files

  #output file
  modelOut=paste0(outputPrefix, ".rda")

  if(!file.exists(modelOut)){
    file.create(modelOut, showWarnings = TRUE)
  }

  #phentoype file
  if(!file.exists(phenoFile)){
    stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
  }else{
    ydat = data.table:::fread(phenoFile, header=T, stringsAsFactors=FALSE, colClasses=list(character = sampleIDColinphenoFile))
    data = data.frame(ydat)

    for(i in c(phenoCol, covarColList, qCovarCol, sampleIDColinphenoFile)){
      if(!(i %in% colnames(data))){
        stop("ERROR! column for ", i, " does not exsit in the phenoFile \n")
      }
    }

    if(length(covarColList) > 0){
      formula = paste0(phenoCol,"~", paste0(covarColList,collapse="+"))
      hasCovariate = TRUE
    }else{
      formula = paste0(phenoCol,"~ 1")
      hasCovariate = FALSE
    }    
    cat("formula is ", formula,"\n")
    formula.null = as.formula(formula)
    mmat = model.frame(formula.null, data, na.action=NULL)
    mmat$IID = data[,which(sampleIDColinphenoFile == colnames(data))]

    #mmat_nomissing0 = mmat[complete.cases(mmat),]
    mmat_nomissing = mmat[complete.cases(mmat),]
    #cat(nrow(mmat_nomissing0), " samples have non-missing phenotypes\n")
    cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")
  }

#  cat("mmat_nomissing0\n")
#  print(mmat_nomissing0)

#  if(!file.exists(sampleFileForDosages)){
#    stop("ERROR! sampleFile ", sampleFileForDosages, " does not exsit\n")
#  }else{
#    sampleListinDosage = data.frame(data.table:::fread(sampleFileForDosages, header=F, stringsAsFactors=FALSE))
#    sampleListinDosage$IndexDose = seq(1,nrow(sampleListinDosage), by=1)
#    cat(nrow(sampleListinDosage), " sample IDs are found in sample file\n")
#    colnames(sampleListinDosage)[1] = "IIDDose"
#    mmat_nomissing = merge(mmat_nomissing0, sampleListinDosage, by.x="IID", by.y="IIDDose")	
#  }
#  mmat_nomissing = mmat_nomissing0
#  cat("mmat_nomissing\n")
#  print(mmat_nomissing)


  if(traitType == "quantitative"){
    if(invNormalize){
      cat("Perform the inverse nomalization for ", phenoCol, "\n")
      invPheno = qnorm((rank(mmat_nomissing[,which(colnames(mmat_nomissing) == phenoCol)], na.last="keep")-0.5)/sum(!is.na(mmat_nomissing[,which(colnames(mmat_nomissing) == phenoCol)])))
      mmat_nomissing[,which(colnames(mmat_nomissing) == phenoCol)] = invPheno
    }
  }


  if(!is.null(kins)){
    if(isDiagofKinSetAsOne){
      diag(kins) = 1
    }
    if(methodforRelatedSample == "EMMAX"){
      if(traitType == "quantitative"){
        out.obj = SKAT:::SKAT_NULL_emmaX(formula.null, data = mmat_nomissing, K=kins)
      }else{
        stop("SKAT_NULL_emmaX does not work for binary tratis \n")
      }
    }else if(methodforRelatedSample == "GMMAT"){
      if(traitType == "quantitative"){
        out.obj = GMMAT:::glmmkin(formula.null, data = mmat_nomissing, family=gaussian(link = "identity"), kins=as.matrix(kins), verbose=T)
      }else{
        out.obj = GMMAT:::glmmkin(formula.null, data = mmat_nomissing, family=binomial(link = "logit"), kins=as.matrix(kins), verbose=T) 
      }
    }
  }else{
    if(traitType == "quantitative"){
      out.obj = SKAT:::SKAT_Null_Model(formula.null, data = mmat_nomissing, out_type="C")
    }else{
      out.obj = SKAT:::SKAT_Null_Model(formula.null, data = mmat_nomissing, out_type="D")
    }
  }
  out.obj$sampleID = mmat_nomissing$IID
  out.obj$traitType = traitType

  save(out.obj, file = modelOut)  
}

