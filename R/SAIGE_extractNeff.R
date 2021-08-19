#' Estimate the effectize sample size 
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in he plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios
#' @param phenoFile character. Path to the phenotype file. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null GLMM. Please note that covariates to be used in the NULL GLMM need to specified using the argument covarColList.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary" or "quantitative". By default, "binary"
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. e.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null GLM model e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the glm model (NOT work yet)
#' @param sampleIDColinphenoFile character. Column name for the sample IDs in the phenotype file e.g. "IID".  
#' @param tol numeric. The tolerance for fitting the null GLMMM to converge. By default, 0.02.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param tolPCG numeric. The tolerance for PCG to converge. By default, 1e-5.
#' @param maxiterPCG integer. The maximum number of iterations for PCG. By default, 500. 
#' @param nThreads integer. Number of threads to be used. By default, 1 
#' @param SPAcutoff numeric. The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2.
#' @param numMarkers integer (>0). Minimum number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed 
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param tauInit vector of numbers. e.g. c(1,1), Unitial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is 0,0, the second element will be 0.5 for binary traits and the initial tau vector for quantitative traits is 1,0 
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. By default, TRUE
#' @param traceCVcutoff numeric. The threshold for coefficient of variantion (CV) for the trace estimator to increase nrun. By default, 0.0025
#' @param ratioCVcutoff numeric. The threshold for coefficient of variantion (CV) for the variance ratio estimate. If ratioCV > ratioCVcutoff. numMarkers will be increased by 10. By default, 0.001 
#' @param outputPrefix character. Path to the output files with prefix.
#' @param outputPrefix_varRatio character. Path to the output variance ratio file with prefix. variace ratios will be output to outputPrefix_varRatio.varianceRatio.txt. If outputPrefix_varRatio is not specified, outputPrefix_varRatio will be the same as the outputPrefix
#' @param IsOverwriteVarianceRatioFile logical. Whether to overwrite the variance ratio file if the file exists. By default, FALSE
#' @param IsSparseKin logical. Whether to exploit the sparsity of GRM to estimate the variance ratio. By default, TRUE
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. If not specified and IsSparseKin=TRUE, sparse GRM will be computed
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. 
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples that are included in the sparse GRM. By default, 2000
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated if IsSparseKin is TRUE. By default, 0.125
#' @param cateVarRatioIndexVec vector of integer 0 or 1. The length of cateVarRatioIndexVec is the number of MAC categories for variance ratio estimation. 1 indicates variance ratio in the MAC category is to be estimated, otherwise 0. By default, NULL. If NULL, variance ratios corresponding to all specified MAC categories will be estimated. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(1.5,2.5,3.5,4.5,5.5,10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param isDiagofKinSetAsOne logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param useSparseSigmaforInitTau logical. Whether to use sparse GRM to estimate the initial values for fitting the null GLMM. By default, FALSE
#' @param useSparseSigmaConditionerforPCG logical. Whether to use sparse GRM to construct a conditoner for PCG. By default, FALSE. Current this option is deactivated.   
#' @param useSparseGRMtoFitNULL logical. Whether to use sparse GRM to fit the null GLMM. By default, FALSE
#' @param minCovariateCount integer. If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues. By default, -1 (no covariates will be excluded)
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @param includeNonautoMarkersforVarRatio logical. Whether to allow for non-autosomal markers for variance ratio. By default, FALSE
#' @param FemaleOnly logical. Whether to run Step 1 for females only. If TRUE, sexCol and FemaleCode need to be specified. By default, FALSE
#' @param MaleOnly logical. Whether to run Step 1 for males only. If TRUE, sexCol and MaleCode need to be specified. By default, FALSE
#' @param FemaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for females. By default, '1' 
#' @param MaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for males. By default, '0'
#' @param sexCol character. Coloumn name for sex in the phenotype file, e.g Sex. By default, '' 
#' @param noEstFixedEff logical. Whether to estimate fixed effect coeffciets. By default, FALSE.  
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
getNeff = function(plinkFile = "", 
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
		IsOverwriteVarianceRatioFile=FALSE,
		IsSparseKin=FALSE,
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
		useSparseSigmaforInitTau = FALSE,
		minCovariateCount = -1, 
		minMAFforGRM = 0.01,
		useSparseGRMtoFitNULL=FALSE,
		includeNonautoMarkersforVarRatio = FALSE,
		sexCol = "",
		FemaleCode = 1,
		FemaleOnly = FALSE,
		MaleCode = 0,	
		MaleOnly = FALSE,
		noEstFixedEff = FALSE,
		skipVarianceRatioEstimation = TRUE)
{
    setminMAFforGRM(minMAFforGRM)
    if (minMAFforGRM > 0) {
        cat("Markers in the Plink file with MAF >= ", minMAFforGRM, 
            " will be used to construct GRM\n")
    }
    else {
        cat("Markers in the Plink file with MAF > ", minMAFforGRM, 
            " will be used to construct GRM\n")
    }
    if (useSparseGRMtoFitNULL) {
        cat("sparse GRM will be used to fit the NULL model\n")
        if (!file.exists(sparseGRMFile)) {
            stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
        }
        if (!file.exists(sparseGRMSampleIDFile)) {
            stop("sparseGRMSampleIDFile ", sparseGRMSampleIDFile, 
                " does not exist!")
        }
        useSparseSigmaforInitTau = FALSE
        useSparseSigmaConditionerforPCG = FALSE
        LOCO = FALSE
        #cat("Leave-one-chromosome-out is not applied\n")
    }
    useSparseSigmaConditionerforPCG = FALSE
    if (useSparseSigmaConditionerforPCG) {
        cat("sparse sigma will be used as the conditioner for PCG\n")
        if (!file.exists(sparseGRMFile)) {
            stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
        }
        if (!file.exists(sparseGRMSampleIDFile)) {
            stop("sparseGRMSampleIDFile ", sparseGRMSampleIDFile, 
                " does not exist!")
        }
    }
    if (useSparseSigmaforInitTau) {
        cat("sparse sigma will be used to estimate the inital tau\n")
        if (!file.exists(sparseGRMFile)) {
            stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
        }
        if (!file.exists(sparseGRMSampleIDFile)) {
            stop("sparseGRMSampleIDFile ", sparseGRMSampleIDFile, 
                " does not exist!")
        }
    }
    if (useSparseSigmaConditionerforPCG | useSparseSigmaforInitTau | 
        useSparseGRMtoFitNULL) {
        IsSparseKin = TRUE
    }
    if (nThreads > 1) {
        RcppParallel:::setThreadOptions(numThreads = nThreads)
        cat(nThreads, " threads are set to be used ", "\n")
    }
    if (FemaleOnly & MaleOnly) {
        stop("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job\n")
    }
    if (FemaleOnly) {
        outputPrefix = paste0(outputPrefix, "_FemaleOnly")
        cat("Female-specific model will be fitted. Samples coded as ", 
            FemaleCode, " in the column ", sexCol, " in the phenotype file will be included\n")
    }
    else if (MaleOnly) {
        outputPrefix = paste0(outputPrefix, "_MaleOnly")
        cat("Male-specific model will be fitted. Samples coded as ", 
            MaleCode, " in the column ", sexCol, " in the phenotype file will be included\n")
    }
    modelOut = paste0(outputPrefix, ".rda")
    SPAGMMATOut = paste0(outputPrefix, "_", numMarkers, "markers.SAIGE.results.txt")
    if (is.null(outputPrefix_varRatio)) {
        outputPrefix_varRatio = outputPrefix
    }
    varRatioFile = paste0(outputPrefix_varRatio, ".varianceRatio.txt")
    #if (!file.exists(varRatioFile)) {
    #    file.create(varRatioFile, showWarnings = TRUE)
    #}
    #else {
    #    if (!IsOverwriteVarianceRatioFile) {
    #        stop("WARNING: The variance ratio file ", varRatioFile, 
    #            " already exists. The new variance ratios will be output to ", 
    #            varRatioFile, ". In order to avoid overwriting the file, please remove the ", 
    #            varRatioFile, " or use the argument outputPrefix_varRatio to specify a different prefix to output the variance ratio(s). Otherwise, specify IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten with new variance ratio(s)\n")
    #    }
    #    else {
    #        cat("The variance ratio file ", varRatioFile, " already exists. IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten\n")
    #    }
    #}
    #if (!file.exists(modelOut)) {
    #    file.create(modelOut, showWarnings = TRUE)
    #}
    if (skipVarianceRatioEstimation & useSparseGRMtoFitNULL) {
        print("a sparse GRM will be used to fit the null model and the variance ratio estimation will be skipped, so plink files are not required")
        sampleListwithGenov0 = data.table:::fread(sparseGRMSampleIDFile, 
            header = F, , colClasses = c("character"), data.table = F)
        colnames(sampleListwithGenov0) = c("IIDgeno")
        sampleListwithGeno = NULL
        sampleListwithGeno$IIDgeno = sampleListwithGenov0$IIDgeno
        sampleListwithGeno = data.frame(sampleListwithGeno)
        sampleListwithGeno$IndexGeno = seq(1, nrow(sampleListwithGeno), 
            by = 1)
        cat(nrow(sampleListwithGeno), " samples are in the sparse GRM\n")
    }
    else {
        if (!file.exists(paste0(plinkFile, ".bed"))) {
            stop("ERROR! ", plinkFile, ".bed does not exsit\n")
        }
        if (!file.exists(paste0(plinkFile, ".bim"))) {
            stop("ERROR! ", plinkFile, ".bim does not exsit\n")
        }
        else {
            chromosomeStartIndexVec = NULL
            chromosomeEndIndexVec = NULL
            if (LOCO) {
                cat("WARNING: leave-one-chromosome-out is activated! Note this option will only be applied to autosomal variants\n")
                cat("WARNING: Genetic variants needs to be ordered by chromosome and position in the Plink file\n")
                bimData = data.table:::fread(paste0(plinkFile, 
                  ".bim"), header = F)
                for (i in 1:22) {
                  if (length(which(bimData[, 1] == i)) > 0) {
                    chromosomeStartIndexVec = c(chromosomeStartIndexVec, 
                      min(which(bimData[, 1] == i)) - 1)
                    chromosomeEndIndexVec = c(chromosomeEndIndexVec, 
                      max(which(bimData[, 1] == i)) - 1)
                    if (i > 1) {
                      if (!is.na(chromosomeStartIndexVec[i - 
                        1])) {
                        if (chromosomeStartIndexVec[i] <= chromosomeStartIndexVec[i - 
                          1] | chromosomeEndIndexVec[i] <= chromosomeEndIndexVec[i - 
                          1]) {
                          stop(paste0("ERROR! chromosomes need to be ordered from 1 to 22 in ", 
                            plinkFile, ".bim\n"))
                        }
                      }
                    }
                  }
                  else {
                    chromosomeStartIndexVec = c(chromosomeStartIndexVec, 
                      NA)
                    chromosomeEndIndexVec = c(chromosomeEndIndexVec, 
                      NA)
                  }
                }
                cat("chromosomeStartIndexVec: ", chromosomeStartIndexVec, 
                  "\n")
                cat("chromosomeEndIndexVec: ", chromosomeEndIndexVec, 
                  "\n")
                if (sum(!is.na(chromosomeStartIndexVec)) <= 1 | 
                  sum(!is.na(chromosomeEndIndexVec)) <= 1) {
                  cat("WARNING: The number of autosomal chromosomes is less than 2 and leave-one-chromosome-out can't be conducted! \n")
                  LOCO = FALSE
                }
                chromosomeStartIndexVec_forcpp = chromosomeStartIndexVec
                chromosomeStartIndexVec_forcpp[is.na(chromosomeStartIndexVec_forcpp)] = -1
                chromosomeEndIndexVec_forcpp = chromosomeEndIndexVec
                chromosomeEndIndexVec_forcpp[is.na(chromosomeEndIndexVec_forcpp)] = -1
                setStartEndIndexVec(chromosomeStartIndexVec_forcpp, 
                  chromosomeEndIndexVec_forcpp)
            }
            else {
                chromosomeStartIndexVec = rep(NA, 22)
                chromosomeEndIndexVec = rep(NA, 22)
            }
        }
        if (!file.exists(paste0(plinkFile, ".fam"))) {
            stop("ERROR! ", plinkFile, ".fam does not exsit\n")
        }
        else {
            sampleListwithGenov0 = data.table:::fread(paste0(plinkFile, 
                ".fam"), header = F, , colClasses = list(character = 1:4))
            sampleListwithGenov0 = data.frame(sampleListwithGenov0)
            colnames(sampleListwithGenov0) = c("FIDgeno", "IIDgeno", 
                "father", "mother", "sex", "phe")
            sampleListwithGeno = NULL
            sampleListwithGeno$IIDgeno = sampleListwithGenov0$IIDgeno
            sampleListwithGeno = data.frame(sampleListwithGeno)
            sampleListwithGeno$IndexGeno = seq(1, nrow(sampleListwithGeno), 
                by = 1)
            cat(nrow(sampleListwithGeno), " samples have genotypes\n")
        }
    }
    if (!file.exists(phenoFile)) {
        stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
    }
    else {
        if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
            ydat = data.table:::fread(cmd = paste0("gunzip -c ", 
                phenoFile), header = T, stringsAsFactors = FALSE, 
                colClasses = list(character = sampleIDColinphenoFile))
        }
        else {
            ydat = data.table:::fread(phenoFile, header = T, 
                stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile))
        }
        data = data.frame(ydat)
        for (i in c(phenoCol, covarColList, qCovarCol, sampleIDColinphenoFile)) {
            if (!(i %in% colnames(data))) {
                stop("ERROR! column for ", i, " does not exist in the phenoFile \n")
            }
        }
        if (FemaleOnly | MaleOnly) {
            if (!sexCol %in% colnames(data)) {
                stop("ERROR! column for sex ", sexCol, " does not exist in the phenoFile \n")
            }
            else {
                if (FemaleOnly) {
                  data = data[which(data[, which(colnames(data) == 
                    sexCol)] == FemaleCode), ]
                  if (nrow(data) == 0) {
                    stop("ERROR! no samples in the phenotype are coded as ", 
                      FemaleCode, " in the column ", sexCol, 
                      "\n")
                  }
                }
                else if (MaleOnly) {
                  data = data[which(data[, which(colnames(data) == 
                    sexCol)] == MaleCode), ]
                  if (nrow(data) == 0) {
                    stop("ERROR! no samples in the phenotype are coded as ", 
                      MaleCode, " in the column ", sexCol, "\n")
                  }
                }
            }
        }
        if (length(covarColList) > 0) {
            formula = paste0(phenoCol, "~", paste0(covarColList, 
                collapse = "+"))
            hasCovariate = TRUE
        }
        else {
            formula = paste0(phenoCol, "~ 1")
            hasCovariate = FALSE
        }
        cat("formula is ", formula, "\n")
        formula.null = as.formula(formula)
        mmat = model.frame(formula.null, data, na.action = NULL)
        mmat$IID = data[, which(sampleIDColinphenoFile == colnames(data))]
        mmat_nomissing = mmat[complete.cases(mmat), ]
        mmat_nomissing$IndexPheno = seq(1, nrow(mmat_nomissing), 
            by = 1)
        cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")
        dataMerge = merge(mmat_nomissing, sampleListwithGeno, 
            by.x = "IID", by.y = "IIDgeno")
        dataMerge_sort = dataMerge[with(dataMerge, order(IndexGeno)), 
            ]
        if (nrow(dataMerge_sort) < nrow(sampleListwithGeno)) {
            cat(nrow(sampleListwithGeno) - nrow(dataMerge_sort), 
                " samples in geno file do not have phenotypes\n")
        }
        cat(nrow(dataMerge_sort), " samples will be used for analysis\n")
    }
    if (traitType == "quantitative" & invNormalize) {
        cat("Perform the inverse nomalization for ", phenoCol, 
            "\n")
        invPheno = qnorm((rank(dataMerge_sort[, which(colnames(dataMerge_sort) == 
            phenoCol)], na.last = "keep") - 0.5)/sum(!is.na(dataMerge_sort[, 
            which(colnames(dataMerge_sort) == phenoCol)])))
        dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] = invPheno
    }
    if (traitType == "binary" & (length(covarColList) > 0)) {
        out_checksep = checkPerfectSep(formula.null, data = dataMerge_sort, 
            minCovariateCount)
        covarColList <- covarColList[!(covarColList %in% out_checksep)]
        formula = paste0(phenoCol, "~", paste0(covarColList, 
            collapse = "+"))
        formula.null = as.formula(formula)
        if (length(covarColList) == 1) {
            hasCovariate = FALSE
        }
        else {
            hasCovariate = TRUE
        }
        dataMerge_sort <- dataMerge_sort[, !(names(dataMerge_sort) %in% 
            out_checksep)]
    }
    if (!hasCovariate) {
        noEstFixedEff = FALSE
    }
    if (noEstFixedEff) {
        print("noEstFixedEff=TRUE, so fixed effects coefficnets won't be estimated.")
        if (traitType == "binary") {
            modwitcov = glm(formula.null, data = dataMerge_sort, 
                family = binomial)
            covoffset = modwitcov$linear.predictors
            dataMerge_sort$covoffset = covoffset
        }
        else {
            modwitcov = lm(formula.null, data = dataMerge_sort)
            dataMerge_sort[, which(colnames(dataMerge_sort) == 
                phenoCol)] = modwitcov$residuals
        }
        formula_nocov = paste0(phenoCol, "~ 1")
        formula.null = as.formula(formula_nocov)
        hasCovariate = FALSE
    }
    if (isCovariateTransform & hasCovariate) {
        cat("qr transformation has been performed on covariates\n")
        out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
        formulaNewList = c("Y ~ ", out.transform$Param.transform$X_name[1])
        if (length(out.transform$Param.transform$X_name) > 1) {
            for (i in c(2:length(out.transform$Param.transform$X_name))) {
                formulaNewList = c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
            }
        }
        formulaNewList = paste0(formulaNewList, collapse = "")
        formulaNewList = paste0(formulaNewList, "-1")
        formula.new = as.formula(paste0(formulaNewList, collapse = ""))
        data.new = data.frame(cbind(out.transform$Y, out.transform$X1))
        colnames(data.new) = c("Y", out.transform$Param.transform$X_name)
        cat("colnames(data.new) is ", colnames(data.new), "\n")
        cat("out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr), 
            "\n")
    }else {
        formula.new = formula.null
        data.new = dataMerge_sort
        out.transform = NULL
    }
    if (IsSparseKin) {
        sparseGRMtest = getsubGRM(sparseGRMFile, sparseGRMSampleIDFile, 
        relatednessCutoff, dataMerge_sort$IID)
        m4 = gen_sp_v2(sparseGRMtest)
        print("print m4")
        print(dim(m4))
        A = summary(m4)
        locationMatinR = rbind(A$i - 1, A$j - 1)
        valueVecinR = A$x
        setupSparseGRM(dim(m4)[1], locationMatinR, valueVecinR)
        rm(sparseGRMtest)
    }
    if (traitType == "binary") {
        cat(phenoCol, " is a binary trait\n")
        uniqPheno = sort(unique(dataMerge_sort[, which(colnames(dataMerge_sort) == 
            phenoCol)]))
        if (uniqPheno[1] != 0 | uniqPheno[2] != 1) {
            stop("ERROR! phenotype value needs to be 0 or 1 \n")
        }
        if (!noEstFixedEff) {
            fit0 = glm(formula.new, data = data.new, family = binomial)
        }
        else {
            fit0 = glm(formula.new, data = data.new, offset = covoffset, 
                family = binomial)
        }
	setgeno(plinkFile, dataMerge_sort$IndexGeno, memoryChunk,
                FALSE)
	setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
        cat("glm:\n")
        print(summary(fit0))
        obj.noK = NULL
	tauVec_ss = c(0, 1)
	print(length(fit0$y))
        wVec_ss = rep(1, length(fit0$y))
        bVec_ss = rep(1, length(fit0$y))
        #bVec_ss = c(1, rep(0, length(fit0$y)-1))
        Rinv_1 = getPCG1ofSigmaAndVector(wVec_ss, tauVec_ss,bVec_ss, maxiterPCG, tolPCG)
	#print(Rinv_1)
        t1_Rinv_1 = sum(Rinv_1)
        cat("t1_Rinv_1 is ", t1_Rinv_1, "\n")
        Pn = sum(fit0$y == 1)/(length(fit0$y))
        Nglmm = 4 * Pn * (1 - Pn) * t1_Rinv_1
        cat("Nglmm ", Nglmm, "\n")

        closeGenoFile_plink()
    }
}
