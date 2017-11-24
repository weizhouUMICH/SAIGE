# Functon to get working vector and fixed & random coefficients
# Run iterations to get converged alpha and eta
Get_Coef = function(y, X, tau, family, alpha0, eta0,  offset, maxiterPCG, tolPCG,maxiter, verbose=FALSE){
  tol.coef = 0.1
  mu = family$linkinv(eta0)
  mu.eta = family$mu.eta(eta0)
  Y = eta0 - offset + (y - mu)/mu.eta
  sqrtW = mu.eta/sqrt(family$variance(mu))
  W = sqrtW^2

  for(i in 1:maxiter){
    cat("iGet_Coef: ", i, "\n")
    re.coef = getCoefficients(Y, X, W, tau, maxiter=maxiterPCG, tol=tolPCG)
    alpha = re.coef$alpha
    eta = re.coef$eta + offset

    if(verbose) {
      cat("Tau:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    mu = family$linkinv(eta)
    mu.eta = family$mu.eta(eta)

    Y = eta - offset + (y - mu)/mu.eta
    sqrtW = mu.eta/sqrt(family$variance(mu))
    W = sqrtW^2

    if( max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol.coef))< tol.coef) break
      alpha0 = alpha
    }

    re = list(Y=Y, alpha=alpha, eta=eta, W=W, cov=re.coef$cov, sqrtW=sqrtW, Sigma_iY = re.coef$Sigma_iY, Sigma_iX = re.coef$Sigma_iX, mu=mu)
}


test_stdGeno = function(subSampleInGeno){
  re1 = system.time({setgeno(genofile, subSampleInGeno)})
  for(itest in 1:1){
    cat(itest, " Get_OneSNP_Geno ", Get_OneSNP_Geno(itest), "\nlength ", length(Get_OneSNP_Geno(itest)), "\n")
  }
}

#Fit the null glmm for binary traits
glmmkin.ai_PCG_Rcpp_Binary = function(genofile, fit0, tau = c(0,0), fixtau = c(0,0), maxiter =20, tol = 0.02, verbose = TRUE, Is.Trace.New=TRUE, nrun=30, tolPCG = 1e-5, maxiterPCG = 500, subPheno) {
  subSampleInGeno = subPheno$IndexGeno
  #print(subSampleInGeno[1:100])


  print("Start reading genotype plink file here")
  re1 = system.time({setgeno(genofile, subSampleInGeno)})
  print("Genotype reading is done")

  y = fit0$y
  n = length(y)
  X = model.matrix(fit0)
  offset = fit0$offset
  if(is.null(offset)) offset = rep(0, n)
  family = fit0$family
  eta = fit0$linear.predictors
  mu = fit0$fitted.values
  mu.eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu)/mu.eta

  alpha0 = fit0$coef
  eta0 = eta

  if(family$family %in% c("poisson", "binomial")) {
    tau[1] = 1
    fixtau[1] = 1
  }

  #change, use 0.5 as a default value, and use Get_Coef before getAIScore
  q = 1
  #tau[fixtau == 0] <- var(Y)/(q+1)
  tau[fixtau == 0] = 0.5
  tau0=tau

  re.coef = Get_Coef(y, X, tau, family, alpha0, eta0,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
  re = getAIScore(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG)
  tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)

  if(verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  for (i in seq_len(maxiter)) {

    if(verbose) cat("\nIteration ", i, tau, ":\n")
    alpha0 = re.coef$alpha
    tau0 = tau

    # use Get_Coef before getAIScore        
    re.coef = Get_Coef(y, X, tau, family, alpha0, eta0,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
    fit = fitglmmaiRPCG(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol)

    tau = as.numeric(fit$tau)
    cov = re.coef$cov
    alpha = re.coef$alpha
    eta = re.coef$eta
    Y = re.coef$Y
    mu = re.coef$mu

    if(tau[2] == 0) break
      # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
      if(max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
      if(max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      	i = maxiter
      	break
      }
  }

  if(verbose) cat("\nFinal " ,tau, ":\n")
  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu

  return(list(theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = subPheno$IID))
}




glmmkin.ai_PCG_Rcpp_Quantitative = function(genofile,fit0, tau = c(0,0), fixtau = c(0,0), maxiter = 20, tol = 0.02, verbose = TRUE, Is.Trace.New=TRUE, nrun=30, tolPCG = 1e-5, maxiterPCG = 500, subPheno) {

  subSampleInGeno = subPheno$IndexGeno
  print("Start reading genotype plink file here")
  re1 = system.time({setgeno(genofile, subSampleInGeno)})
  print("Genotype reading is done")

  y = fit0$y
  n = length(y)
  offset = fit0$offset
  if(is.null(offset)) offset = rep(0, n)
  family = fit0$family
  eta = fit0$linear.predictors
  mu = fit0$fitted.values
  mu.eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu)/mu.eta
  sqrtW = mu.eta/sqrt(fit0$family$variance(mu))
  W = sqrtW^2
  X = model.matrix(fit0)

  X1 = SPAtest:::ScoreTest_wSaddleApprox_Get_X1(X)


  #print("X")
  #print(X)
  alpha = fit0$coef
  if(verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }

  if(family$family %in% c("poisson", "binomial")) {
    tau[1] = 1
    fixtau[1] = 1
  }


  q = 1
  cat("tau ",tau,"\n")
  tau[fixtau == 0] = var(Y)/(q+1)
  tau0 = tau

  cat("tauv2 ",tau,"\n")

  print("ok1")
  re = getAIScore_q(Y, X, W, tau, nrun, maxiterPCG, tolPCG)
  tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace[2])/n)
  tau[1] = max(0, tau0[1] + tau0[1]^2 * (re$YPA0PY - re$Trace[1])/n)
  cat("tauv3 ",tau,"\n")
  for (i in seq_len(maxiter)) {
    W = sqrtW^2

    if(verbose) cat("\nIteration ", i, ":\n")
    alpha0 = alpha

    tau0 = tau
    fit = fitglmmaiRPCG_q(Y, X, W, tau, nrun, maxiterPCG, tolPCG, tol)
    tau = as.numeric(fit$tau)
    cov = as.matrix(fit$cov)
    alpha = as.numeric(fit$alpha)
    eta = as.numeric(fit$eta) + offset
    if(verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    mu = family$linkinv(eta)
    mu.eta = family$mu.eta(eta)
    Y = eta - offset + (y - mu)/mu.eta
    sqrtW = mu.eta/sqrt(family$variance(mu))

    if(2*max(max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol)), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
    if(max(tau) > tol^(-2)) {
      warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      i = maxiter
      break
    }
  }

  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu
  Sigma_iy = getSigma_G(W, tau, res, maxiterPCG, tolPCG)
  Sigma_iX = getSigma_X(W, tau, X1, maxiterPCG, tolPCG)

  return(list(theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged, sampleID = subPheno$IID, Sigma_iy = Sigma_iy, Sigma_iX = Sigma_iX))
}


Saddle_Prob_q <-function(q, mu, g, tauVecNew){
  m1 = sum(mu * g)
  var2 = sum(g^2)
  pval.noadj = pchisq(((q - m1)/tauVecNew[1])^2/var2, lower.tail = FALSE, df=1)
  return(list(p.value = pval.noadj, p.value.NA = NA, Is.converge = NA, p1 = NA, p2 = NA))
}



ScoreTest_wSaddleApprox_NULL_Model_q=function (formula, data = NULL){
  X1 = model.matrix(formula, data = data)
  X1 = SPAtest:::ScoreTest_wSaddleApprox_Get_X1(X1)
  glmfit = glm(formula, data = data, family=gaussian(link = "identity"))
  mu = glmfit$fitted.values
  V = 1
  res = glmfit$y - mu
  n1 = length(res)
  XV = t(X1 * V)
  XVX_inv = solve(t(X1) %*% (X1 * V))
  XXVX_inv = X1 %*% XVX_inv
  re = list(y = glmfit$y, mu = mu, res = res, V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv)
  class(re) = "SA_NULL"
  return(re)
}


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
#' @param centerVariables vector of characters.  Covariates to be centered (around the mean) e.g. c("birthYear")
#' @param nThreads integer. Number of threads to be used. By default, 1 
#' @param numMarkers integer (>0). Number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether tp skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed 
#' @param outputPrefix character. Path to the output files with prefix. 
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio value, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM = function(plinkFile = "", 
                phenoFile = "",
                phenoCol = "",
                traitType = "binary",
                invNormalize = FALSE,
                covarColList = NULL,
                qCovarCol = NULL,
                sampleIDColinphenoFile = "",
                centerVariables=NULL,
                tol=0.02,
                maxiter=20,
                tolPCG=1e-5,
                maxiterPCG=500,
                nThreads = 1, 
                Cutoff = 2, 
                numMarkers = 30, 
                skipModelFitting = FALSE,
                outputPrefix = ""){
                #formula, phenoType = "binary",prefix, centerVariables = "", tol=0.02, maxiter=20, tolPCG=1e-5, maxiterPCG=500, nThreads = 1, Cutoff = 2, numMarkers = 1000, skipModelFitting = FALSE){
  if(nThreads > 1){
    RcppParallel:::setThreadOptions(numThreads = nThreads)
    cat(nThreads, " threads are set to be used ", "\n")
  }

  #check and read files
  #output file
  modelOut=paste0(outputPrefix, ".rda")
  SPAGMMATOut=paste0(outputPrefix, "_", numMarkers,"markers.SAIGE.results.txt")
  varRatioFile=paste0(outputPrefix,".varianceRatio.txt")

  if(!file.exists(modelOut)){
    file.create(modelOut, showWarnings = TRUE)
  }

  if(!file.exists(paste0(plinkFile, ".bed"))){
    stop("ERROR! ", plinkFile, ".bed does not exsit\n")
  }

  if(!file.exists(paste0(plinkFile, ".bim"))){
    stop("ERROR! ", plinkFile, ".bim does not exsit\n")
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
    ydat = data.table:::fread(phenoFile, header=T, stringsAsFactors=FALSE)
    data = data.frame(ydat)

    for(i in c(phenoCol, covarColList, qCovarCol, sampleIDColinphenoFile)){
      if(!(i %in% colnames(data))){
        stop("ERROR! column for ", i, " does not exsit in the phenoFile \n")
      }
    }

    #update the categorical variables
    qCovarColUpdated = NULL
    for(i in qCovarCol){
      j = paste0("factor(", i, ")")
      qCovarColUpdated = c(qCovarColUpdated, j)
    }

    formula = paste0(phenoCol,"~",paste0(c(covarColList,qCovarColUpdated),collapse="+"))
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
  }

  #center some covariates
  if(length(centerVariables)!=0){
    for(i in centerVariables){
      if (!(i %in% colnames(dataMerge_sort))){
        stop("ERROR! column for ", i, " does not exsit in the phenoFile \n")
      }else{
        dataMerge_sort[,which(colnames(dataMerge_sort) == i)] = dataMerge_sort[,which(colnames(dataMerge_sort) == i)] - mean(dataMerge_sort[,which(colnames(dataMerge_sort) == i)])
      }
    }
  }


  if(traitType == "binary"){
    cat(phenoCol, " is a binary trait\n")
    uniqPheno = sort(unique(dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)]))
    if (uniqPheno[1] != 0 | uniqPheno[2] != 1){
      stop("ERROR! phenotype value needs to be 0 or 1 \n")
    }
    fit0 = glm(formula.null,data=dataMerge_sort, family=binomial)
    print(fit0)
    obj.noK = SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula.null, data = dataMerge_sort)


    if(!skipModelFitting){
      system.time(modglmm<-glmmkin.ai_PCG_Rcpp_Binary(plinkFile, fit0, tau = c(0,0), fixtau = c(0,0), maxiter =maxiter, tol = tol, verbose = TRUE, Is.Trace.New=TRUE, nrun=30, tolPCG = tolPCG, maxiterPCG = maxiterPCG, subPheno = dataMerge_sort))
      save(modglmm, file = modelOut)
    }else{
      setgeno(plinkFile, dataMerge_sort$IndexGeno)
      load(modelOut)
    }
    scoreTest_SPAGMMAT_forVarianceRatio_binaryTrait(obj.glmm.null = modglmm,
                                                    obj.glm.null = fit0,
                                                    obj.noK = obj.noK,
                                                    Cutoff = Cutoff,
                                                    maxiterPCG = maxiterPCG,
                                                    tolPCG = tolPCG,
                                                    numMarkers = numMarkers,
                                                    varRatioOutFile = varRatioFile,
                                                    testOut = SPAGMMATOut)
    closeGenoFile_plink()

  }else if(traitType == "quantitative"){

    cat(phenoCol, " is a quantitative trait\n")
    if(invNormalize){
      cat("Perform the inverse nomalization for ", phenoCol, "\n")
      invPheno = qnorm((rank(dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)], na.last="keep")-0.5)/sum(!is.na(dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)])))
      dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)] = invPheno
    }

    obj.noK = ScoreTest_wSaddleApprox_NULL_Model_q(formula.null, dataMerge_sort)
    fit0 = glm(formula.null, data=dataMerge_sort,family=gaussian(link = "identity"))
 
    if(!skipModelFitting){

      system.time(modglmm<-glmmkin.ai_PCG_Rcpp_Quantitative(plinkFile,fit0, tau = c(0,0), fixtau = c(0,0), maxiter =maxiter, tol = tol, verbose = TRUE, Is.Trace.New=TRUE, nrun=30, tolPCG = tolPCG, maxiterPCG = maxiterPCG, subPheno = dataMerge_sort))
      save(modglmm, file = modelOut)
      print("step2")
    }else{
      setgeno(plinkFile, dataMerge_sort$IndexGeno)
      load(modelOut)
    }
 
    scoreTest_SPAGMMAT_forVarianceRatio_quantitativeTrait(obj.glmm.null = modglmm,
                                                    obj.glm.null = fit0,
                                                    obj.noK = obj.noK,
                                                    Cutoff = Cutoff,
                                                    maxiterPCG = maxiterPCG,
                                                    tolPCG = tolPCG,
                                                    numMarkers = numMarkers,
                                                    varRatioOutFile = varRatioFile,
                                                    testOut = SPAGMMATOut)
    closeGenoFile_plink()
  }
}


scoreTest_SPAGMMAT_forVarianceRatio_binaryTrait = function(obj.glmm.null,
                                                    obj.glm.null,
                                                    obj.noK,
                                                    Cutoff = 2,
                                                    maxiterPCG = 500,
                                                    tolPCG = 0.01,
                                                    numMarkers,
                                                    varRatioOutFile,
                                                    testOut){

  if(file.exists(testOut)){file.remove(testOut)}
  #resultHeader = c("markerIndex","p.value", "p.value.NA", "Is.converge","var1","var2", "N", "NCase", "NCtrl", "AC", "AC.Case", "AC.Ctrl", "AF", "AF.Case", "AF.Ctrl")
  resultHeader = c("markerIndex","p.value", "p.value.NA", "Is.converge","varT","varTstar", "N", "AC", "AF")
  write(resultHeader,file = testOut, ncolumns = length(resultHeader))

  if(Cutoff < 10^-2){
    Cutoff = 10^-2
  }

  family = obj.glm.null$family
  print(family)


  print(names(obj.glmm.null))
  eta = obj.glmm.null$linear.predictors
  mu = obj.glmm.null$fitted.values

  mu.eta = family$mu.eta(eta)
  sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
  W = sqrtW^2
  tauVecNew = obj.glmm.null$theta
  X1 = obj.noK$X1
  Sigma_iX = getSigma_X(W, tauVecNew, X1, maxiterPCG, tolPCG)
  y = obj.glm.null$y

  ##randomize the marker orders to be tested
  mMarkers = gettotalMarker()
  listOfMarkersForVarRatio = sample(c(1:mMarkers), size = mMarkers, replace = FALSE)
  freqVec = getAlleleFreqVec()
  Nnomissing = length(mu)

  OUTtotal = NULL
  OUT = NULL

  indexInMarkerList = 1
  numTestedMarker = 0
  while(numTestedMarker < numMarkers){
    i = listOfMarkersForVarRatio[indexInMarkerList]
    cat("i is ", i, "\n")
    G0 = Get_OneSNP_Geno(i-1)
    #cat("G0", G0[1:10], "\n")
   
    if(sum(G0)/(2*Nnomissing) > 0.5){
      G0 = 2-G0
    }
    NAset = which(G0==0)
    AC = sum(G0)

    #if (AC <= 20 | AC >= (2*Nnomissing - 20)){
    if (AC <= 20){
      indexInMarkerList = indexInMarkerList + 1
    }else{
      AF = AC/(2*Nnomissing)
      #NCase = sum(y == 1)
      #NCtrl = sum(y == 0)
      #AC.Case = sum(G0[which(y == 1)])
      #AC.Ctrl = sum(G0[which(y == 0)])
      #AF.Case = AC.Case/(2*NCase)
      #AF.Ctrl = AC.Ctrl/(2*NCtrl)
      G = G0  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G0) # G1 is X adjusted 
      g = G/sqrt(AC)
      q = innerProduct(g,y)
      Sigma_iG = getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)
      var1a = t(G)%*%Sigma_iG - t(G)%*%Sigma_iX%*%(solve(t(X1)%*%Sigma_iX))%*%t(X1)%*%Sigma_iG
      var1 = var1a/AC
      m1 = innerProduct(mu,g)
      var2 = innerProduct(mu*(1-mu), g*g)
      qtilde = (q-m1)/sqrt(var1) * sqrt(var2) + m1

      if(length(NAset)/length(G) < 0.5){
        out1 = SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
      }else {
        out1 = SPAtest:::Saddle_Prob_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8)
      }

      #out1 = SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
      #OUT = rbind(OUT, c(i, out1$p.value, out1$p.value.NA, out1$Is.converge, var1, var2, Nnomissing, NCase, NCtrl, AC, AC.Case, AC.Ctrl,AF, AF.Case, AF.Ctrl))
      OUT = rbind(OUT, c(i, out1$p.value, out1$p.value.NA, out1$Is.converge, var1, var2, Nnomissing, AC, AF))
      indexInMarkerList = indexInMarkerList + 1
      numTestedMarker = numTestedMarker + 1
      if(numTestedMarker %% 10 == 0 | numTestedMarker == numMarkers){
        OUT = as.data.frame(OUT)
        OUTtotal = rbind(OUTtotal, OUT)
        write.table(OUT, testOut, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
        OUT = NULL
      }
    }
  }

  OUTtotal = as.data.frame(OUTtotal)
  colnames(OUTtotal) = resultHeader

  OUT1 = OUTtotal
  varRatio = mean(OUT1$varT/OUT1$varTstar)
  cat("varRatio", varRatio, "\n")
  write(varRatio, varRatioOutFile)
  print(varRatio)
}


scoreTest_SPAGMMAT_forVarianceRatio_quantitativeTrait = function(obj.glmm.null,
                                                    obj.glm.null,
                                                    obj.noK,
                                                    Cutoff = 2,
                                                    maxiterPCG = 500,
                                                    tolPCG = 0.01,
                                                    numMarkers,
                                                    varRatioOutFile,
                                                    testOut){

  if(file.exists(testOut)){file.remove(testOut)}

  #resultHeader = c("markerIndex","p.value", "p.value.NA","var1","var2","Tv1", "Tv2", "p.value.Tv2","N", "AC", "AF")
  resultHeader = c("markerIndex","p.value", "p.value.NA","var1","var2","Tv1","N", "AC", "AF")

  write(resultHeader,file = testOut, ncolumns = length(resultHeader))

  if(Cutoff < 10^-2){
    Cutoff = 10^-2
  }

  family = obj.glm.null$family
  print(family)

  eta = obj.glmm.null$linear.predictors
  mu = obj.glmm.null$fitted.values
  mu.eta = family$mu.eta(eta)
  sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
  W = sqrtW^2
  tauVecNew = obj.glmm.null$theta

  X1 = obj.noK$X1
  Sigma_iX = getSigma_X(W, tauVecNew, X1, maxiterPCG, tolPCG)
  y = obj.glm.null$y

  ##randomize the marker orders to be tested
  mMarkers = gettotalMarker()
  listOfMarkersForVarRatio = sample(c(1:mMarkers), size = mMarkers, replace = FALSE)
 # listOfMarkersForVarRatio = c(1:mMarkers)
  freqVec = getAlleleFreqVec()
  Nnomissing = length(mu)

  OUTtotal = NULL
  OUT = NULL

  indexInMarkerList = 1
  numTestedMarker = 0
  while(numTestedMarker < numMarkers){
    i = listOfMarkersForVarRatio[indexInMarkerList]
    cat("i is ", i, "\n")
    G0 = Get_OneSNP_Geno(i-1)
    cat("G0", G0[1:10], "\n")
    AC = sum(G0)
    if (AC <= 20 | AC >= (2*Nnomissing - 20)){
      indexInMarkerList = indexInMarkerList + 1
    }else{
      AF = AC/(2*Nnomissing)
      G = G0  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G0) # G1 is X adjusted 
      g = G/sqrt(AC)
      q = innerProduct(g,y)
 #     print(g[1:20])
 #     print(y[1:20])
 #     print(q)
      Sigma_iG = getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)
      var1a = t(G)%*%Sigma_iG - t(G)%*%Sigma_iX%*%(solve(t(X1)%*%Sigma_iX))%*%t(X1)%*%Sigma_iG
      ###var1 = g'Pg, var2 = g'g
      var1 = var1a/AC
      m1 = innerProduct(mu,g)
#      m1 = sum(mu * g)
      #mu2 = 1-mu
      #innermumu2= innerProduct(mu, mu2)
      #innerg = innerProduct(g, g)
      #var2 = innerProduct(innermumu2, innerg)
      var2 = innerProduct(g, g)
      #qtilde = (q-m1)/sqrt(var1) * sqrt(var2) + m1
      #out1 = SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
      Tv1 = (q-m1)/tauVecNew[1]
#      cat("q ",q, " m1 ", m1, " tauVecNew[1] ", tauVecNew[1], " Tv1 ", Tv1, "\n")
      p.value = pchisq(Tv1^2/var1, lower.tail = FALSE, df=1)
      p.value.NA = pchisq(Tv1^2/var2, lower.tail = FALSE, df=1)
      #####Tversion2 = g'Pytilde
      #Tv2 = t(g)%*%Pytilde
      #p.value.Tv2 = pchisq(Tv2^2/var1, lower.tail = FALSE, df=1)

      #OUT = rbind(OUT, c(mth,p.value, p.value.NA, var1, var2, Tv1, Tv2, p.value.Tv2, Nnomissing, AC, AF))
      OUT = rbind(OUT, c(i, p.value, p.value.NA, var1, var2, Tv1, Nnomissing, AC, AF))

      indexInMarkerList = indexInMarkerList + 1
      numTestedMarker = numTestedMarker + 1
      if(numTestedMarker %% 10 == 0 | numTestedMarker == numMarkers){
        OUT = as.data.frame(OUT)
        OUTtotal = rbind(OUTtotal, OUT)
        write.table(OUT, testOut, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
        OUT = NULL
      }
    }
  }

  OUTtotal = as.data.frame(OUTtotal)
  colnames(OUTtotal) = resultHeader
  OUT1 = OUTtotal
  varRatio = mean(OUT1$var1/OUT1$var2)
  cat("varRatio", varRatio, "\n")
  write(varRatio, varRatioOutFile)
  print(varRatio)
}
