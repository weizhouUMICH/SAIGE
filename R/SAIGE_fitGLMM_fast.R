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



# Functon to get working vector and fixed & random coefficients when LOCO is TRUE
# getCoefficients_LOCO is used
# Run iterations to get converged alpha and eta
Get_Coef_LOCO = function(y, X, tau, family, alpha0, eta0,  offset, maxiterPCG, tolPCG, maxiter, verbose=FALSE){
  tol.coef = 0.1
  mu = family$linkinv(eta0)
  mu.eta = family$mu.eta(eta0)
  Y = eta0 - offset + (y - mu)/mu.eta
  sqrtW = mu.eta/sqrt(family$variance(mu))
  W = sqrtW^2

  for(i in 1:maxiter){
    cat("iGet_Coef: ", i, "\n")
    re.coef = getCoefficients_LOCO(Y, X, W, tau, maxiter=maxiterPCG, tol=tolPCG)
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


#Fits the null glmm for binary traits
glmmkin.ai_PCG_Rcpp_Binary = function(genofile, fit0, tau=c(0,0), fixtau = c(0,0), maxiter =20, tol = 0.02, verbose = TRUE, nrun=30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne) {
  #Fits the null generalized linear mixed model for a binary trait
  #Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix 
  #  fit0: glm model. Logistic model output (with no sample relatedness accounted for) 
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the glmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model  
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates  
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome  
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  #Returns:
  #  model output for the null glmm

  subSampleInGeno = subPheno$IndexGeno
  if(verbose){
    print("Start reading genotype plink file here")
  }

  re1 = system.time({setgeno(genofile, subSampleInGeno, memoryChunk, isDiagofKinSetAsOne)})

  if(verbose){
    print("Genotype reading is done")
  }

  y = fit0$y
  n = length(y)
  X = model.matrix(fit0)

  offset = fit0$offset
  if(is.null(offset)){
    offset = rep(0, n)
  }

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

  if(tauInit[fixtau == 0] == 0){
  tau[fixtau == 0] = 0.5
  }else{
    tau[fixtau == 0] = tauInit[fixtau == 0]
  }
  cat("inital tau is ", tau,"\n")
  tau0=tau

  re.coef = Get_Coef(y, X, tau, family, alpha0, eta0,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
  re = getAIScore(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG,tolPCG = tolPCG, traceCVcutoff = traceCVcutoff)


  tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)

  if(verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  for (i in seq_len(maxiter)) {

    if(verbose) cat("\nIteration ", i, tau, ":\n")
    alpha0 = re.coef$alpha
    tau0 = tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 = eta

    # use Get_Coef before getAIScore        
    re.coef = Get_Coef(y, X, tau, family, alpha0, eta0,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
    fit = fitglmmaiRPCG(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol, traceCVcutoff = traceCVcutoff)

    tau = as.numeric(fit$tau)
    cov = re.coef$cov
    alpha = re.coef$alpha
    eta = re.coef$eta
    Y = re.coef$Y
    mu = re.coef$mu

     print(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol))
      cat("tau: ", tau, "\n")
      cat("tau0: ", tau0, "\n")


    if(tau[2] == 0) break
      # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
      if(max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
      #print(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol))	
      #cat("tau: ", tau, "\n")
      #cat("tau0: ", tau0, "\n")

      if(max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      	i = maxiter
      	break
      }
  }

  if(verbose) cat("\nFinal " ,tau, ":\n")

    #added these steps after tau is estimated 04-14-2018
  re.coef = Get_Coef(y, X, tau, family, alpha, eta,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
  cov = re.coef$cov
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu

  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu
  if(isCovariateTransform){
  coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  }else{
    coef.alpha = alpha
  }
  glmmResult = list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = subPheno$IID, obj.noK=obj.noK, obj.glm.null=fit0, traitType="binary")

  #LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme  

  glmmResult$LOCO = LOCO
  if(LOCO){
    glmmResult$LOCOResult = list()
     
    for (j in 1:22){
      startIndex = chromosomeStartIndexVec[j]
      endIndex = chromosomeEndIndexVec[j]
      if(!is.na(startIndex) && !is.na(endIndex)){
        setStartEndIndex(startIndex, endIndex)
        re.coef_LOCO = Get_Coef_LOCO(y, X, tau, family, alpha, eta,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter)
        cov = re.coef_LOCO$cov
        alpha = re.coef_LOCO$alpha
        eta = re.coef_LOCO$eta
        Y = re.coef_LOCO$Y
        mu = re.coef_LOCO$mu
        res = y - mu
        if(isCovariateTransform){
        coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
	}else{
	coef.alpha = alpha
	}
        glmmResult$LOCOResult[[j]] = list(isLOCO = TRUE, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov)
      }else{
        glmmResult$LOCOResult[[j]] = list(isLOCO = FALSE)
      }
    }
  }

  return(glmmResult)
}



#Fits the null glmm for a quantitative trait
glmmkin.ai_PCG_Rcpp_Quantitative = function(genofile, fit0, tau = c(0,0), fixtau = c(0,0), maxiter = 20, tol = 0.02, verbose = TRUE, nrun=30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne){
  #Fits the null linear mixed model for a quantitative trait
  #Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. linear model output (with no sample relatedness accounted for), family=gaussian(link = "identity")
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the lmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  #Returns:
  #  model output for the null lmm


  subSampleInGeno = subPheno$IndexGeno
  if(verbose){
    print("Start reading genotype plink file here")
  }

  re1 = system.time({setgeno(genofile, subSampleInGeno, memoryChunk, isDiagofKinSetAsOne)})

  #test time
#  for(a in 1:100){
#  time1 = proc.time()
#  Get_OneSNP_StdGeno(a-1)
#  time2 = proc.time()
#  timediff = time2 - time1
#  print(timediff)
#  }
#  break

  #test time

  if(verbose){
    print("Genotype reading is done")
  }

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
#  cat("sqrtW: ",sqrtW,"\n")
  W = sqrtW^2
#  cat("fit0\n")
#  print(fit0)
  X = model.matrix(fit0)
#  cat("X\n")
#  print(X)


  X1 = SPAtest:::ScoreTest_wSaddleApprox_Get_X1(X)

#  print("X")
#  print(X)
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
  if(sum(tauInit[fixtau == 0]) == 0){
    tau[fixtau == 0] = var(Y)/(q+1)
  }else{
    tau[fixtau == 0] = tauInit[fixtau == 0]
  }

  tau0 = tau
  cat("initial tau is ", tau,"\n")
  #bvtest = rep(3.08474, n)
#  A = NULL
#  for(i in c(1:n)){
#  	bvtest = rep(0, n)
#  	bvtest[i] = 1
#  	a = getCrossprodMatAndKin(bvtest)
#  	A = c(A, a[i])
#  }
#  A = matrix(A, ncol=1)
#  write.table(A, "/net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation_08_2018/jobs/SAIGE_SKATO/step1/jobs/diagOfKin.txt", quote=F, row.names=F, col.names=F)

  #bvtest = rep(0.05, n)
  #a = getCrossprodMatAndKin(bvtest)
  #print(a)

  #bvtest = rep(5, n)
  #a = getCrossprodMatAndKin(bvtest)
  #print(a)


  #bvtest = rep(1, n)
  #a = getCrossprodMatAndKin(bvtest)
  #print(a)


  re = getAIScore_q(Y, X, W, tau, nrun, maxiterPCG, tolPCG, traceCVcutoff)
#  cat(names(re))
#  cat("X\n")
#  print(X)
#  cat("Sigma_iX:\n")
#  print(re$Sigma_iX[1:20,])
#  cat("PY\n")
#  print(re$PY)
#  cat("Trace\n")
#  print(re$Trace)
#  cat("YPAPY:\n")
#  print(re$YPAPY)
#  cat("YPA0PY:\n")
#  print(re$YPA0PY)
  #print(sum((re$PY/W)^2))
  tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace[2])/n)
  tau[1] = max(0, tau0[1] + tau0[1]^2 * (re$YPA0PY - re$Trace[1])/n)
  #tau[1] = max(0, tau0[1] + tau0[1]^2 * (sum((re$PY/W)^2) - re$Trace[1])/n) #try 
  #testVec=rep(1,100)
  #testVec[50] = 1
  #testVecResult=getCrossprodMatAndKin(testVec)
  #cat("testVecResult\n")

  #print(testVecResult)
  #cat("tauv3 ",tau,"\n")


  if(verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }


  for (i in seq_len(maxiter)) {
    W = sqrtW^2

    if(verbose) cat("\nIteration ", i, ":\n")
    alpha0 = alpha
    tau0 = tau
#    cat("tau0: ", tau0,"\n")
    fit = fitglmmaiRPCG_q(Y, X, W, tau, nrun, maxiterPCG, tolPCG, tol, traceCVcutoff)
#    cat("tau0_after_fit: ", tau0,"\n")
#    print(fit)
    tau = as.numeric(fit$tau)
    cov = as.matrix(fit$cov)
    cat("cov: ", cov, "\n")
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
if(FALSE){
    cat("abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol)\n")
    print(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol))	
    cat("tau: ", tau,"\n")
    cat("tau0: ", tau0,"\n")
    cat("abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)\n")
    print(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol))
    cat("tol: ")
    print(tol)

    print(2*max(max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol)), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)))
    print(2*max(max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol)), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol)
}

    if(tau[2] == 0) break
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
  #cat("Sigma_iX: ", Sigma_iX, "\n")

  if(isCovariateTransform){
    coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  }else{
    coef.alpha = alpha
  }

  #coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)


  lmmResult = list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged, sampleID = subPheno$IID, Sigma_iy = Sigma_iy, Sigma_iX = Sigma_iX, obj.noK=obj.noK, obj.glm.null=fit0, traitType="quantitative")

  #LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme
  lmmResult$LOCO = LOCO  
  if(LOCO){
    lmmResult$LOCOResult = list()

    for (j in 1:22){
      startIndex = chromosomeStartIndexVec[j]
      endIndex = chromosomeEndIndexVec[j]
      if(!is.na(startIndex) && !is.na(endIndex)){
        setStartEndIndex(startIndex, endIndex)

        re.coef_LOCO=fitglmmaiRPCG_q_LOCO(Y, X, W, tau, nrun, maxiterPCG, tolPCG, tol, traceCVcutoff)

	cov = as.matrix(re.coef_LOCO$cov)
    	alpha = as.numeric(re.coef_LOCO$alpha)
    	eta = as.numeric(re.coef_LOCO$eta) + offset
    	mu = family$linkinv(eta)
    	mu.eta = family$mu.eta(eta)
    	Y = eta - offset + (y - mu)/mu.eta
        res = y - mu
        if(isCovariateTransform){
        coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
        }else{
        coef.alpha = alpha
        }

        #coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
        lmmResult$LOCOResult[[j]] = list(isLOCO = TRUE, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov)
      }else{
        lmmResult$LOCOResult[[j]] = list(isLOCO = FALSE)
      }
    }
  }

  return(lmmResult)

  #return(list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged, sampleID = subPheno$IID, Sigma_iy = Sigma_iy, Sigma_iX = Sigma_iX, obj.noK=obj.noK, obj.glm.null=fit0, traitType="quantitative"))
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
  re = list(y = glmfit$y, mu = mu, res = res, V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv)
  class(re) = "SA_NULL"
  return(re)
}


#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants 
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM). Genetic markers are also randomly selected from the plink file to estimate the variance ratios
#' @param phenoFile character. Path to the phenotype file. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Addiitonal columns can be included in the phenotype file for covariates in the null GLMM. Please note covariates to be used in the NULL GLMM need to specified using the argument covarColList.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary" or "quantitative". By default, "binary"
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. E.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null GLM model e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the glm model (NOT work yet)
#' @param sampleIDColinphenoFile character.  Column name for the sample IDs in the phenotype file e.g. "IID".  
#' @param tol numeric.The tolerance for fitting the null GLMMM to converge. By default, 0.02.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param tolPCG numeric. The tolerance for PCG to converge. By default, 1e-5.
#' @param maxiterPCG integer. The maximum number of iterations for PCG. By default, 500. 
#' @param nThreads integer. Number of threads to be used. By default, 1 
#' @param SPAcutoff numeric. The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2.
#' @param numMarkers integer (>0). Minimum number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed 
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param tauInit vector of numbers. e.g. c(1,1), Unitial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is not specified, the second element will be 0.5 for binary traits.  
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. 
#' @param traceCVcutoff numeric. The threshold for coefficient of variantion (CV) for the trace estimator to increase nrun. By default, 0.0025
#' @param ratioCVcutoff numeric. The threshold for coefficient of variantion (CV) for the variance ratio estimate. If ratioCV > ratioCVcutoff. numMarkers will be increased by 10. By default, 0.001 
#' @param outputPrefix character. Path to the output files with prefix.
#' @param outputPrefix_varRatio character. Path to the output variance ratio file with prefix. variace ratios will be output to outputPrefix_varRatio.varianceRatio.txt. If outputPrefix_varRatio is not specified, outputPrefix_varRatio will be the same as the outputPrefix
#' @param IsSparseKin logical. Whether to exploit the sparsity of GRM to estimate the variance ratio. By default, TRUE
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. 
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples for sparse GRM. By default, 1000
#' @param isCateVarianceRatio logical. Whether to estimate variance ratio based on different MAC categories. If yes, variance ratio will be estiamted for multiple MAC categories corresponding to cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude. Currently, if isCateVarianceRatio=TRUE, then LOCO=FALSE. By default=FALSE 
#' @param relatednessCutoff float. The threshold to treat two samples as unrelated if IsSparseKin is TRUE. By default, 0.125
#' @param cateVarRatioIndexVec vector of integer 0 or 1. The length of cateVarRatioIndexVec is the number of MAC categories for variance ratio estimation. 1 indicates variance ratio in the MAC category is to be estimated, otherwise 0. By default, NULL. If NULL, variance ratios corresponding to all specified MAC categories will be estimated. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(1.5,2.5,3.5,4.5,5.5,10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param isDiagofKinSetAsOne logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM = function(plinkFile = "", 
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
		isDiagofKinSetAsOne = FALSE){


  
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
    ydat = data.table:::fread(phenoFile, header=T, stringsAsFactors=FALSE)
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

#  data.new = data.frame(cbind(out.transform$Y, out.transform$X1))
#  colnames(data.new) = c("Y",out.transform$Param.transform$X_name)
#  cat("colnames(data.new) is ", colnames(data.new), "\n")
#  cat("out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr), "\n")


  if(traitType == "binary"){
    cat(phenoCol, " is a binary trait\n")
  #  uniqPheno = sort(unique(dataMerge_sort[,which(colnames(dataMerge_sort) == phenoCol)]))
    uniqPheno = sort(unique(out.transform$Y))
    if (uniqPheno[1] != 0 | uniqPheno[2] != 1){
      stop("ERROR! phenotype value needs to be 0 or 1 \n")
    }
    #fit0 = glm(formula.null,data=dataMerge_sort, family=binomial)
    #fit0 = glm(out.transform$Y ~ out.transform$X1,family=binomial)
    fit0 = glm(formula.new, data=data.new, family=binomial)
    cat("glm:\n")
    print(fit0)
    
    #obj.noK = SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula.null, data = dataMerge_sort)
    obj.noK = SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula.new, data = data.new)


    if(!skipModelFitting){
      cat("Start fitting the NULL GLMM\n")
      t_begin = proc.time()
      print(t_begin)

      system.time(modglmm<-glmmkin.ai_PCG_Rcpp_Binary(plinkFile, fit0, tau = c(0,0), fixtau = c(0,0), maxiter =maxiter, tol = tol, verbose = TRUE, nrun=30, tolPCG = tolPCG, maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, obj.noK = obj.noK, out.transform = out.transform, tauInit=tauInit, memoryChunk=memoryChunk, LOCO=LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec, chromosomeEndIndexVec = chromosomeEndIndexVec, traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform, isDiagofKinSetAsOne = isDiagofKinSetAsOne))
      save(modglmm, file = modelOut)

      t_end = proc.time()
      print(t_end)
      cat("t_end - t_begin, fitting the NULL model took\n")
      print(t_end - t_begin)

    }else{
      cat("Skip fitting the NULL GLMM\n")
      load(modelOut)
      if(is.null(modglmm$LOCO)){modglmm$LOCO = FALSE}
      setgeno(plinkFile, dataMerge_sort$IndexGeno, memoryChunk, isDiagofKinSetAsOne)	
      set.seed(98765)
n <- 5e3
# 5000 x 5000 matrices, 99% sparse
a <- rsparsematrix(n, n, 0.01, rand.x=function(n) rpois(n, 1) + 1)
b <- rsparsematrix(n, n, 0.01, rand.x=function(n) rpois(n, 1) + 1)

d=a %*% b
print("print d")
print(dim(d))

m1 <- mult_sp_sp_to_sp(a, b)
print("print m1")
print(dim(m1))

m2 = gen_sp(a)
print("print m2")
print(dim(m2))

sparseGRMtest = Matrix:::readMM(sparseGRMFile)

m3 = gen_sp_v2(a)
print("print m3")
print(dim(m3))

m4 = gen_sp_v2(sparseGRMtest)
print("print m4")
print(dim(m4))
A = summary(m4)

locationMatinR = rbind(A$i-1, A$j-1)
valueVecinR = A$x
setupSparseGRM(length(A$x), locationMatinR, valueVecinR)
B = gen_sp_GRM()
print("print B")
print(dim(B))

x = gen_spsolve_v3()
print(x[1:30])
    }
    cat("Start estimating variance ratios\n")
    scoreTest_SPAGMMAT_forVarianceRatio_binaryTrait(obj.glmm.null = modglmm,
                                                    obj.glm.null = fit0,
                                                    obj.noK = obj.noK,
                                                    Cutoff = SPAcutoff,
                                                    maxiterPCG = maxiterPCG,
                                                    tolPCG = tolPCG,
                                                    numMarkers = numMarkers,
                                                    varRatioOutFile = varRatioFile,
						    ratioCVcutoff = ratioCVcutoff,
                                                    testOut = SPAGMMATOut,
						    plinkFile = plinkFile,
						    chromosomeStartIndexVec = chromosomeStartIndexVec, 
						    chromosomeEndIndexVec = chromosomeEndIndexVec,
						    isCateVarianceRatio = isCateVarianceRatio,
						    cateVarRatioIndexVec = cateVarRatioIndexVec,
                                                    IsSparseKin = IsSparseKin,
                                                    sparseGRMFile = sparseGRMFile,
                                                    sparseGRMSampleIDFile = sparseGRMSampleIDFile,
                                                    numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
                                                    relatednessCutoff = relatednessCutoff,
                                                    nThreads = nThreads)
    closeGenoFile_plink()

  }else if(traitType == "quantitative"){

    cat(phenoCol, " is a quantitative trait\n")
 
    obj.noK = ScoreTest_wSaddleApprox_NULL_Model_q(formula.new, data.new)
    fit0 = glm(formula.new, data=data.new,family=gaussian(link = "identity"))
    cat("glm:\n")
    print(fit0)


    if(!skipModelFitting){
      cat("Start fitting the NULL GLMM\n")
      t_begin = proc.time()
      print(t_begin)

      system.time(modglmm<-glmmkin.ai_PCG_Rcpp_Quantitative(plinkFile,fit0, tau = c(0,0), fixtau = c(0,0), maxiter =maxiter, tol = tol, verbose = TRUE, nrun=30, tolPCG = tolPCG, maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, obj.noK=obj.noK, out.transform=out.transform, tauInit=tauInit, memoryChunk = memoryChunk, LOCO=LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec, chromosomeEndIndexVec = chromosomeEndIndexVec, traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform, isDiagofKinSetAsOne = isDiagofKinSetAsOne))
      save(modglmm, file = modelOut)

      t_end = proc.time()
      print(t_end)
      cat("t_end - t_begin, fitting the NULL model took\n")
      print(t_end - t_begin)
      print("step2")

    }else{

      cat("Skip fitting the NULL GLMM\n")
      load(modelOut)
      if(is.null(modglmm$LOCO)){modglmm$LOCO = FALSE}
      setgeno(plinkFile, dataMerge_sort$IndexGeno, memoryChunk, isDiagofKinSetAsOne)

      #test time
#	btest = rnorm(nrow(data.new))
#	for(i in 1:10){
#		tTimeVec = testTime(i, btest)	
#	}
#      if(is.null(modglmm$LOCO)){modglmm$LOCO = FALSE}
    }
#    cat("dataMerge_sort$IndexGeno: ", dataMerge_sort$IndexGeno, "\n") 

    cat("Start estimating variance ratios\n")
    scoreTest_SPAGMMAT_forVarianceRatio_quantitativeTrait(obj.glmm.null = modglmm,
                                                    obj.glm.null = fit0,
                                                    obj.noK = obj.noK,
                                                    Cutoff = SPAcutoff,
                                                    maxiterPCG = maxiterPCG,
                                                    tolPCG = tolPCG,
                                                    numMarkers = numMarkers,
                                                    varRatioOutFile = varRatioFile,
						    ratioCVcutoff = ratioCVcutoff,
                                                    testOut = SPAGMMATOut,
						    plinkFile = plinkFile,
                                                    chromosomeStartIndexVec = chromosomeStartIndexVec,
                                                    chromosomeEndIndexVec = chromosomeEndIndexVec,
						    isCateVarianceRatio = isCateVarianceRatio,
						    cateVarRatioIndexVec = cateVarRatioIndexVec,
						    IsSparseKin = IsSparseKin,
						    sparseGRMFile = sparseGRMFile,
                				    sparseGRMSampleIDFile = sparseGRMSampleIDFile,	
						    numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
						    relatednessCutoff = relatednessCutoff,
						    nThreads = nThreads)
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
						    ratioCVcutoff,
                                                    testOut,
                                                    plinkFile,
						    chromosomeStartIndexVec, 
						    chromosomeEndIndexVec,
						    isCateVarianceRatio,
						    cateVarRatioIndexVec,
                                                    IsSparseKin,
                                                    sparseGRMFile,
                                                    sparseGRMSampleIDFile,
                                                    numRandomMarkerforSparseKin,
                                                    relatednessCutoff,
                                                    nThreads){


  if(file.exists(testOut)){file.remove(testOut)}
  


  resultHeader = c("CHR","SNPID","POS","A1","A2","p.value", "p.value.NA", "Is.converge","var1","var2", "N", "AC", "AF")
  write(resultHeader,file = testOut, ncolumns = length(resultHeader))
  bimPlink = data.frame(data.table:::fread(paste0(plinkFile,".bim"), header=F))


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
  Sigma_iX_noLOCO = getSigma_X(W, tauVecNew, X1, maxiterPCG, tolPCG)
  y = obj.glm.null$y
  ##randomize the marker orders to be tested

    #####sparse Kin

  if(IsSparseKin){
    sparseSigma = getSparseSigma(outputPrefix=varRatioOutFile,
                sparseGRMFile=sparseGRMFile,
                sparseGRMSampleIDFile=sparseGRMSampleIDFile,
                numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
                relatednessCutoff = relatednessCutoff,
                obj.glmm.null = obj.glmm.null,
                W=W, tauVecNew=tauVecNew)
  }

  mMarkers = gettotalMarker()
#  listOfMarkersForVarRatio = sample(c(1:mMarkers), size = mMarkers, replace = FALSE)
  listOfMarkersForVarRatio = list()
  MACvector = getMACVec()
#  freqVec = getAlleleFreqVec()
#  Nnomissing = length(mu)

#  OUTtotal = NULL
#  OUT = NULL
#  indexInMarkerList = 1
#  numTestedMarker = 0
#  ratioCV = ratioCVcutoff + 0.1

  if(!isCateVarianceRatio){
    cat("Only one variance ratio will be estimated using randomly selected markers with MAC >= 20\n")
    MACindex = which(MACvector >= 20)
    listOfMarkersForVarRatio[[1]] = sample(MACindex, size = length(MACindex), replace = FALSE)
    cateVarRatioIndexVec=c(1)
  }else{
    cat("Categorical variance ratios will be estimated\n")

    if(is.null(cateVarRatioIndexVec)){cateVarRatioIndexVec = rep(1, length(cateVarRatioMinMACVecExclude))}
    numCate = length(cateVarRatioIndexVec)
    for(i in 1:(numCate-1)){
       #print("i 1:(numCate-1)")
       #print(i)
       #print(cateVarRatioMinMACVecExclude[i])
       #print(cateVarRatioMaxMACVecInclude[i])
       #print(length(MACvector))
       #print(MACvector[1:10])
       #print(min(MACvector))

      MACindex = which(MACvector > cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
      #print(length(MACindex))
      #tempindex = which(MACvector > 0.5 & MACvector <= 1.5)
      #print(length(tempindex))

      listOfMarkersForVarRatio[[i]] = sample(MACindex, size = length(MACindex), replace = FALSE)

    }
    if(length(cateVarRatioMaxMACVecInclude) == (numCate-1)){
      MACindex = which(MACvector > cateVarRatioMinMACVecExclude[numCate])
    }else{
      MACindex = which(MACvector > cateVarRatioMinMACVecExclude[numCate] & MACvector <= cateVarRatioMaxMACVecInclude[numCate])
    }
    listOfMarkersForVarRatio[[numCate]] = sample(MACindex, size = length(MACindex), replace = FALSE)

    for(k in 1:length(cateVarRatioIndexVec)){
      if(k <= length(cateVarRatioIndexVec)-1){
        if(cateVarRatioIndexVec[k] == 1){
          cat(cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k],"\n")
          if(length(listOfMarkersForVarRatio[[k]]) < numMarkers){
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], " is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }else{
        if(cateVarRatioIndexVec[k] == 1){
          cat(cateVarRatioMinMACVecExclude[k], "< MAC\n")
          if(length(listOfMarkersForVarRatio[[k]]) < numMarkers){
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC  is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }
    }


  }# if(!isCateVarianceRatio){

  freqVec = getAlleleFreqVec()



  Nnomissing = length(mu)
  varRatioTable = NULL
  Sigma_iX_noLOCO = getSigma_X(W, tauVecNew, X1, maxiterPCG, tolPCG)


  for(k in 1:length(listOfMarkersForVarRatio)){
    #if(length(listOfMarkersForVarRatio[[k]]) == 0){
    #  cateVarRatioIndexVec[k] = 0
    #  cat("no marker is found in the MAC category ", k, "\n")
    #}
    if(cateVarRatioIndexVec[k] == 1){

      numMarkers0 = numMarkers
      OUTtotal = NULL
      OUT = NULL
      indexInMarkerList = 1
      numTestedMarker = 0
      ratioCV = ratioCVcutoff + 0.1

      while(ratioCV > ratioCVcutoff){
        while(numTestedMarker < numMarkers0){
          i = listOfMarkersForVarRatio[[k]][indexInMarkerList]
          cat(i, "th marker\n")
          G0 = Get_OneSNP_Geno(i-1)
          cat("G0", G0[1:10], "\n")
          #AC = sum(G0)
          CHR = bimPlink[i,1]

	  if(sum(G0)/(2*Nnomissing) > 0.5){
            G0 = 2-G0
          }
          NAset = which(G0==0)
          AC = sum(G0)

         if (CHR < 1 | CHR > 22){
           indexInMarkerList = indexInMarkerList + 1
         }else{
          AF = AC/(2*Nnomissing)
          G = G0  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G0) # G1 is X adjusted
          g = G/sqrt(AC)
          q = innerProduct(g,y)
 #     print(g[1:20])
 #     print(y[1:20])
 #     print(q)
          if(!obj.glmm.null$LOCO){
            Sigma_iG = getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)
            Sigma_iX = Sigma_iX_noLOCO
          }else if(!(obj.glmm.null$LOCOResult[[CHR]]$isLOCO)){
             eta = obj.glmm.null$linear.predictors
             mu = obj.glmm.null$fitted.values
             mu.eta = family$mu.eta(eta)
             sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
             W = sqrtW^2
             Sigma_iG = getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)
             Sigma_iX = Sigma_iX_noLOCO
          }else{
             eta = obj.glmm.null$LOCOResult[[CHR]]$linear.predictors
             mu = obj.glmm.null$LOCOResult[[CHR]]$fitted.values
             mu.eta = family$mu.eta(eta)
             sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
             W = sqrtW^2
             startIndex = chromosomeStartIndexVec[CHR]
             endIndex = chromosomeEndIndexVec[CHR]
             setStartEndIndex(startIndex, endIndex)
             Sigma_iG = getSigma_G_LOCO(W, tauVecNew, G, maxiterPCG, tolPCG)
             Sigma_iX = getSigma_X_LOCO(W, tauVecNew, X1, maxiterPCG, tolPCG)
          }

          var1a = t(G)%*%Sigma_iG - t(G)%*%Sigma_iX%*%(solve(t(X1)%*%Sigma_iX))%*%t(X1)%*%Sigma_iG
      ###var1 = g'Pg, var2 = g'g

      #cat("Sigma_iG: \n")
      #print(Sigma_iG/AC)

          var1 = var1a/AC
          m1 = innerProduct(mu,g)


          if(IsSparseKin){
            t1 = proc.time()
            cat("t1\n")
             cat("t1again\n")
#       pcginvSigma = getPCG1ofSparseSigmaAndVector(sparseSigma, g)
#       pcginvSigma = pcgSparse(sparseSigma, g)
             #pcginvSigma = pcg(sparseSigma, g)
             pcginvSigma = solve(sparseSigma, g, sparse=T)
             t2 = proc.time()
             cat("t2-t1\n")
             print(t2-t1)
             var2_a = t(g) %*% pcginvSigma
             var2 = var2_a[1,1]
        #cat("qrinvSigma: \n")
        #print(qrinvSigma)
        }else{
          var2 = innerProduct(mu*(1-mu), g*g)
        }

      var2q = innerProduct(mu*(1-mu), g*g)
      qtilde = (q-m1)/sqrt(var1) * sqrt(var2q) + m1

      if(length(NAset)/length(G) < 0.5){
        out1 = SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
      }else {
        out1 = SPAtest:::Saddle_Prob_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="p")
      }


      OUT = rbind(OUT, c(bimPlink[i,1], bimPlink[i,2], bimPlink[i,4], bimPlink[i,5], bimPlink[i,6], out1$p.value, out1$p.value.NA, out1$Is.converge, var1, var2, Nnomissing, AC, AF))	
#        OUT = rbind(OUT, c(i, p.value, p.value.NA, var1, var2, Tv1, Nnomissing, AC, AF))

      indexInMarkerList = indexInMarkerList + 1
      numTestedMarker = numTestedMarker + 1

      if(numTestedMarker %% 10 == 0 | numTestedMarker == numMarkers | indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]]) ){
          OUT = as.data.frame(OUT)
          print("OK")
          OUTtotal = rbind(OUTtotal, OUT)
          print("OK1")
          write.table(OUT, testOut, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
          OUT = NULL
        }
      }

      if(indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]])){
        numTestedMarker = numMarkers0
      }
    }#end of while(numTestedMarker < numMarkers)

    print("OK2")
    #OUTtotal = as.data.frame(OUTtotal)
    #colnames(OUTtotal) = resultHeader
    OUT1 = OUTtotal
    OUT1 = as.data.frame(OUT1)
    colnames(OUT1) = resultHeader
    ratioVec = as.numeric(OUT1$var1)/as.numeric(OUT1$var2)
    ratioCV = calCV(ratioVec)

    if(ratioCV > ratioCVcutoff){
      cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " > ", ratioCVcutoff, "\n")
      numMarkers0 = numMarkers0 + 10
      cat("try ", numMarkers0, " markers\n")
    }else{
      cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " < ", ratioCVcutoff, "\n")
    }

    if(indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]])){
      ratioCV = ratioCVcutoff
      cat("no more markers are available in the MAC category ", k, "\n")
      print(indexInMarkerList-1)
    }

  }#end of while(ratioCV > ratioCVcutoff)

  #OUTtotal = as.data.frame(OUTtotal)
  #colnames(OUTtotal) = resultHeader
  OUT1 = OUTtotal
  OUT1 = as.data.frame(OUT1)
  colnames(OUT1) = resultHeader
  varRatio = mean(as.numeric(OUT1$var1)/as.numeric(OUT1$var2))
  cat("varRatio", varRatio, "\n")
  varRatioTable = rbind(varRatioTable, c(varRatio))
#  write(varRatio, varRatioOutFile)
  print(varRatio)
  print(varRatioTable)

  }else{# if(cateVarRatioVec[k] == 1)
    varRatioTable = rbind(varRatioTable, c(1))
  }

} #for(k in 1:length(listOfMarkersForVarRatio)){


  print(varRatioTable)
  print(varRatioOutFile)
  write.table(varRatioTable, varRatioOutFile, quote=F, col.names=F, row.names=F)
  data = read.table(varRatioOutFile, header=F)
  print(data)

}




scoreTest_SPAGMMAT_forVarianceRatio_quantitativeTrait = function(obj.glmm.null,
                                                    obj.glm.null,
                                                    obj.noK,
                                                    Cutoff = 2,
                                                    maxiterPCG = 500,
                                                    tolPCG = 0.01,
                                                    numMarkers,
                                                    varRatioOutFile,
						    ratioCVcutoff,
                                                    testOut,
						    plinkFile,
                                                    chromosomeStartIndexVec,
                                                    chromosomeEndIndexVec,
						    isCateVarianceRatio,
						    cateVarRatioIndexVec,
						    IsSparseKin,
						    sparseGRMFile,
                                                    sparseGRMSampleIDFile,
						    numRandomMarkerforSparseKin,
                                                    relatednessCutoff,
						    nThreads){	


  if(file.exists(testOut)){file.remove(testOut)}

#  if(nThreads > 1){
#    RcppParallel:::setThreadOptions(numThreads = nThreads)
#    cat(nThreads, " threads are set to be used ", "\n")
#  }
    
  #resultHeader = c("markerIndex","p.value", "p.value.NA","var1","var2","Tv1", "Tv2", "p.value.Tv2","N", "AC", "AF")
  resultHeader = c("markerIndex","p.value", "p.value.NA","var1","var2","Tv1","N", "AC", "AF")
  write(resultHeader,file = testOut, ncolumns = length(resultHeader))

  bimPlink = data.frame(data.table:::fread(paste0(plinkFile,".bim"), header=F))

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
  y = obj.glm.null$y

    #####sparse Kin

  if(IsSparseKin){
    sparseSigma = getSparseSigma(outputPrefix=varRatioOutFile,
                sparseGRMFile=sparseGRMFile,
                sparseGRMSampleIDFile=sparseGRMSampleIDFile,
                numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
                relatednessCutoff = relatednessCutoff,
                obj.glmm.null = obj.glmm.null,
                W=W, tauVecNew=tauVecNew)
  }


  ##randomize the marker orders to be tested
  mMarkers = gettotalMarker()


  listOfMarkersForVarRatio = list()	
  MACvector = getMACVec()

  if(!isCateVarianceRatio){
    cat("Only one variance ratio will be estimated using randomly selected markers with MAC >= 20\n")
    MACindex = which(MACvector >= 20)
    listOfMarkersForVarRatio[[1]] = sample(MACindex, size = length(MACindex), replace = FALSE)
    cateVarRatioIndexVec=c(1)
  }else{
    cat("Categorical variance ratios will be estimated\n")
    if(is.null(cateVarRatioIndexVec)){cateVarRatioIndexVec = rep(1, length(cateVarRatioMinMACVecExclude))}
    numCate = length(cateVarRatioIndexVec)
    for(i in 1:(numCate-1)){
       #print("i 1:(numCate-1)")
       #print(i)
       #print(cateVarRatioMinMACVecExclude[i])
       #print(cateVarRatioMaxMACVecInclude[i])
       #print(length(MACvector))	
       #print(MACvector[1:10])	
       #print(min(MACvector))

      MACindex = which(MACvector > cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
      #print(length(MACindex))	
      #tempindex = which(MACvector > 0.5 & MACvector <= 1.5)	
      #print(length(tempindex))

      listOfMarkersForVarRatio[[i]] = sample(MACindex, size = length(MACindex), replace = FALSE)

    }
    if(length(cateVarRatioMaxMACVecInclude) == (numCate-1)){
      MACindex = which(MACvector > cateVarRatioMinMACVecExclude[numCate])
    }else{
      MACindex = which(MACvector > cateVarRatioMinMACVecExclude[numCate] & MACvector <= cateVarRatioMaxMACVecInclude[numCate])
    }
    listOfMarkersForVarRatio[[numCate]] = sample(MACindex, size = length(MACindex), replace = FALSE)

    for(k in 1:length(cateVarRatioIndexVec)){
      if(k <= length(cateVarRatioIndexVec)-1){
        if(cateVarRatioIndexVec[k] == 1){
          cat(cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k],"\n")
          if(length(listOfMarkersForVarRatio[[k]]) < numMarkers){
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], " is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }else{
        if(cateVarRatioIndexVec[k] == 1){	
          cat(cateVarRatioMinMACVecExclude[k], "< MAC\n")
          if(length(listOfMarkersForVarRatio[[k]]) < numMarkers){
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC  is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }
    }


  }

  freqVec = getAlleleFreqVec()


  Nnomissing = length(mu)
  varRatioTable = NULL

  Sigma_iX_noLOCO = getSigma_X(W, tauVecNew, X1, maxiterPCG, tolPCG)


  for(k in 1:length(listOfMarkersForVarRatio)){
    #if(length(listOfMarkersForVarRatio[[k]]) == 0){
    #  cateVarRatioIndexVec[k] = 0
    #  cat("no marker is found in the MAC category ", k, "\n")
    #}
    if(cateVarRatioIndexVec[k] == 1){

      numMarkers0 = numMarkers
      OUTtotal = NULL
      OUT = NULL
      indexInMarkerList = 1
      numTestedMarker = 0
      ratioCV = ratioCVcutoff + 0.1

      while(ratioCV > ratioCVcutoff){  
        while(numTestedMarker < numMarkers0){
          i = listOfMarkersForVarRatio[[k]][indexInMarkerList]
          cat(i, "th marker\n")
          G0 = Get_OneSNP_Geno(i-1)
          cat("G0", G0[1:10], "\n")
          AC = sum(G0)
          CHR = bimPlink[i,1]
         if (CHR < 1 | CHR > 22){
           indexInMarkerList = indexInMarkerList + 1
         }else{
          AF = AC/(2*Nnomissing)
          G = G0  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G0) # G1 is X adjusted 
          g = G/sqrt(AC)
          q = innerProduct(g,y)
 #     print(g[1:20])
 #     print(y[1:20])
 #     print(q)
          if(!obj.glmm.null$LOCO){          
            Sigma_iG = getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)
            Sigma_iX = Sigma_iX_noLOCO
          }else if(!(obj.glmm.null$LOCOResult[[CHR]]$isLOCO)){
             eta = obj.glmm.null$linear.predictors
             mu = obj.glmm.null$fitted.values
             mu.eta = family$mu.eta(eta)
             sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
             W = sqrtW^2
             Sigma_iG = getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)
             Sigma_iX = Sigma_iX_noLOCO
          }else{
             eta = obj.glmm.null$LOCOResult[[CHR]]$linear.predictors
             mu = obj.glmm.null$LOCOResult[[CHR]]$fitted.values
             mu.eta = family$mu.eta(eta)
             sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
             W = sqrtW^2
             startIndex = chromosomeStartIndexVec[CHR]
             endIndex = chromosomeEndIndexVec[CHR]
             setStartEndIndex(startIndex, endIndex)
             Sigma_iG = getSigma_G_LOCO(W, tauVecNew, G, maxiterPCG, tolPCG)
             Sigma_iX = getSigma_X_LOCO(W, tauVecNew, X1, maxiterPCG, tolPCG)
          }

          var1a = t(G)%*%Sigma_iG - t(G)%*%Sigma_iX%*%(solve(t(X1)%*%Sigma_iX))%*%t(X1)%*%Sigma_iG
      ###var1 = g'Pg, var2 = g'g

      #cat("Sigma_iG: \n")
      #print(Sigma_iG/AC)	

          var1 = var1a/AC
          m1 = innerProduct(mu,g)


          if(IsSparseKin){
	    t1 = proc.time()
	    cat("t1\n")
#	pcginvSigma = getPCG1ofSparseSigmaAndVector(sparseSigma, g)
#	pcginvSigma = pcgSparse(sparseSigma, g)
	     #pcginvSigma = pcg(sparseSigma, g)
	     pcginvSigma = solve(sparseSigma, g, sparse=T)
	#print(class(sparseSigma))
	#print(dim(pcginvSigma))
	#print(class(pcginvSigma))
	#require(Matrix)	
	#a1<-methods:::as(sparseSigma, "dsTMatrix")
#	pcginvSigma = pcg(sparseSigma, g)
	#pcginvSigma = pcg(a1, g)
	     t2 = proc.time()
             cat("t2-t1\n")
	     print(t2-t1)
	     var2_a = t(g) %*% pcginvSigma
	     var2 = var2_a[1,1]
	#cat("qrinvSigma: \n")
	#print(qrinvSigma)
          }else{
             var2 = innerProduct(g, g)
          }

        Tv1 = (q-m1)/tauVecNew[1]
        p.value = pchisq(Tv1^2/var1, lower.tail = FALSE, df=1)
        p.value.NA = pchisq(Tv1^2/var2, lower.tail = FALSE, df=1)
        OUT = rbind(OUT, c(i, p.value, p.value.NA, var1, var2, Tv1, Nnomissing, AC, AF))

        indexInMarkerList = indexInMarkerList + 1
        numTestedMarker = numTestedMarker + 1
      	
	
        if(numTestedMarker %% 10 == 0 | numTestedMarker == numMarkers | indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]]) ){
          OUT = as.data.frame(OUT)
	  print("OK")
          OUTtotal = rbind(OUTtotal, OUT)
	  print("OK1")
          write.table(OUT, testOut, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
          OUT = NULL
        }
      }

      if(indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]])){
        numTestedMarker = numMarkers0
      }
    }#end of while(numTestedMarker < numMarkers)

    print("OK2")
    #OUTtotal = as.data.frame(OUTtotal)
    #colnames(OUTtotal) = resultHeader
    OUT1 = OUTtotal
    OUT1 = as.data.frame(OUT1)
    colnames(OUT1) = resultHeader
    ratioVec = as.numeric(OUT1$var1)/as.numeric(OUT1$var2)
    ratioCV = calCV(ratioVec)

    if(ratioCV > ratioCVcutoff){
      cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " > ", ratioCVcutoff, "\n")
      numMarkers0 = numMarkers0 + 10
      cat("try ", numMarkers0, " markers\n")
    }else{
      cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " < ", ratioCVcutoff, "\n")
    }

    if(indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]])){
      ratioCV = ratioCVcutoff
      cat("no more markers are available in the MAC category ", k, "\n")
      print(indexInMarkerList-1)	
    }

  }#end of while(ratioCV > ratioCVcutoff)

  #OUTtotal = as.data.frame(OUTtotal)
  #colnames(OUTtotal) = resultHeader
  OUT1 = OUTtotal
  OUT1 = as.data.frame(OUT1)
  colnames(OUT1) = resultHeader
  varRatio = mean(as.numeric(OUT1$var1)/as.numeric(OUT1$var2))
  cat("varRatio", varRatio, "\n")
  varRatioTable = rbind(varRatioTable, c(varRatio))
#  write(varRatio, varRatioOutFile)
  print(varRatio)
  print(varRatioTable)

  }else{# if(cateVarRatioVec[k] == 1)
    varRatioTable = rbind(varRatioTable, c(1))
  }

} #for(k in 1:length(listOfMarkersForVarRatio)){


  print(varRatioTable)
  print(varRatioOutFile)
  write.table(varRatioTable, varRatioOutFile, quote=F, col.names=F, row.names=F)
  data = read.table(varRatioOutFile, header=F)
  print(data)

}


##suggested by Shawn 01-19-2018
Covariate_Transform<-function(formula, data){
  X1<-model.matrix(formula,data=data)
#  X1=X1[,c(2:ncol(X1))] #remove intercept
  formula.frame<-model.frame(formula,data=data)
  Y = model.response(formula.frame, type = "any")
  X_name = colnames(X1)
		
  # First run linear regression to identify multi collinearity 
  out.lm<-lm(Y ~ X1 - 1, data=data)
#  out.lm<-lm(Y ~ X1, data=data)
  idx.na<-which(is.na(out.lm$coef))
  if(length(idx.na)> 0){
	X1<-X1[, -idx.na]
	X_name = X_name[-idx.na]		
        cat("Warning: multi collinearity is detected in covariates! ", X_name[idx.na], " will be excluded in the model\n")
  }
  if(!(1 %in% idx.na)){
    X_name[1] = "minus1"
  }

	
 # QR decomposition
  Xqr = qr(X1)
  X1_Q = qr.Q(Xqr)
  qrr = qr.R(Xqr)
	
  N<-nrow(X1)
	
  # Make square summation=N (so mean=1)
  X1_new<-X1_Q * sqrt(N)	
  Param.transform<-list(qrr=qrr, N=N, X_name = X_name, idx.na=idx.na)
  re<-list(Y =Y, X1 = X1_new, Param.transform=Param.transform)
}

# In case to recover original scale coefficients
# X \beta = Q R \beta = (Q \sqrt(N)) ( R \beta / \sqrt(N))
# So coefficient from fit.new is the same as R \beta / \sqrt(N)
Covariate_Transform_Back<-function(coef, Param.transform){	
	#coef<-fit.new$coef; Param.transform=out.transform$Param.transform
	coef1<-coef * sqrt(Param.transform$N)
	coef.org<-solve(Param.transform$qrr, coef1)
	
	names(coef.org)<-Param.transform$X_name
	return(coef.org)
}


pcg<-function (A, b, M=NULL, maxiter = 1e+05, tol = 1e-06){
  
  # A<-a; b<-c1[,1]; M<-NULL;maxiter = 1e+05; tol = 1e-06
  if (is.null(M)) {
    dA <- diag(A)
    dA[which(dA == 0)] = 1e-04
    #print("dA")
    #print(dA)
    Minv = 1/dA
  } else Minv = solve(M)
  x = rep(0, length(b))
  r = b 
  if(is.null(M)){
    z = Minv *r
  } else {
    z= Minv %*% r
  }
  p = z
  iter = 0
  sumr2 = sum(r^2)
  while (sumr2 > tol & iter < maxiter) {
    iter = iter + 1
#    cat("iter is ", iter, "\n")
    Ap = crossprod(p, A)[1,]
    a = as.numeric((t(r) %*% z)/(t(p) %*% Ap))
    x = x + a * p
    r1 = r - a * Ap
    
    if(is.null(M)){
      z1 = Minv * r1
    } else {
      z1 = Minv %*% r1
    }
    
    
    bet = as.numeric((t(z1) %*% r1)/(t(z) %*% r))

    p = z1 + bet * p

    z = z1
    r = r1
    sumr2 = sum(r^2)
  }
  if (iter >= maxiter) 
    x = "pcg did not converge. You may increase maxiter number."
  return(x)
}


pcgSparse<-function (A, b, M=NULL, maxiter = 1e+05, tol = 1e-06){
  # A<-a; b<-c1[,1]; M<-NULL;maxiter = 1e+05; tol = 1e-06
  if (is.null(M)) {
    dA <- diag(A)
    dA[which(dA == 0)] = 1e-04
    #print("dA")
    #print(dA)
    Minv = 1/dA
  } else Minv = solve(M)
  x = rep(0, length(b))
  r = b
  if(is.null(M)){
    z = Minv *r
  } else {
    z= Minv %*% r
  }
  p = z
  iter = 0
  sumr2 = sum(r^2)
  print("psparse0")
  psparse =  Matrix:::sparseMatrix(i = rep(1, length(p)), j = c(1:length(p)), x = as.vector(p))
  cat("nrow(psparse) ", nrow(psparse), "\n")
  cat("ncol(psparse) ", ncol(psparse), "\n")
  print(class(psparse))

  while (sumr2 > tol & iter < maxiter) {
    iter = iter + 1
    #Ap = crossprod(p, A)[1,]
    print(class(psparse)[1])	
    print(class(A)[1])	
    print(nrow(psparse))	
    print(ncol(psparse))	
    print(nrow(A))	
    print(ncol(A))	

    Ap = sparse_row_idx_mult(psparse, A)
    #Ap = psparse%*%A
    print("psparse1")
    a = as.numeric((t(r) %*% z)/(t(p) %*% Ap))
    x = x + a * p
    r1 = r - a * Ap

    if(is.null(M)){
      z1 = Minv * r1
    } else {
      z1 = Minv %*% r1
    }


    bet = as.numeric((t(z1) %*% r1)/(t(z) %*% r))

    p = z1 + bet * p

    z = z1
    r = r1
    sumr2 = sum(r^2)
  }
  if (iter >= maxiter)
    x = "pcg did not converge. You may increase maxiter number."
  return(x)
}


#https://gist.github.com/bobthecat/5024079
bigGRMPar = function(x, nblocks = 10, verbose = TRUE, ncore= 1, relatednessCutoff = 0){
#  library(foreach)
#  library(doParallel)
  #register cores
  doParallel:::registerDoParallel(ncore)

  NCOL <- ncol(x)
  NROW <- nrow(x)

  ## test if ncol(x) %% nblocks gives remainder 0
#  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}

#  if(NCOL %% nblocks == 0){
  ## split column numbers into 'nblocks' groups
#  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  SPLIT <- split(1:NCOL, ceiling(seq_along(1:NCOL)/(NCOL%/%nblocks)))
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
#  }else{
#    NCOLNEW = NCOL - (NCOL%%nblocks + NCOL%/%nblocks)
#    SPLIT <- split(1:NCOLNEW, rep(1:(nblocks-1), each = NCOLNEW/(nblocks-1)))
#  }
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  `%dopar%` <- foreach::`%dopar%`
  results <- foreach:::foreach(i = 1:nrow(COMBS), .combine='rbind')%dopar%{
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        GRM <- t(x[, G1])%*%(x[, G2])

        if(sum(G1 != G2) == 0){
                GRM[lower.tri(GRM, diag = FALSE)] = 0
        }
        GRM <- GRM/NROW
        indice <- which(GRM >= relatednessCutoff, arr.ind=T)
        cbind(G1[indice[,1]], G2[indice[,2]])
#       resultsIVec =c(resultsIVec, G1[indice[,1]])
#        resultsJVec =c(resultsJVec, G2[indice[,2]])
        #corMAT[G1, G2] <- COR
        #corMAT[G2, G1] <- t(COR)
#       GRM <- NULL
#       indice = NULL
}
  gc()
  return(results)
}



bigGRMPar_new = function(nblocks = 10, verbose = TRUE, ncore= 1, relatednessCutoff = 0){
  doParallel:::registerDoParallel(ncore)

  NCOL <- getNColStdGenoMultiMarkersMat()
  NROW <- getNRowStdGenoMultiMarkersMat()
  cat("NCOL: ", NCOL, "\n")
  cat("NROW: ", NROW, "\n")

  nblocks <- NCOL%/%10

  ## test if ncol(x) %% nblocks gives remainder 0
#  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}

#  if(NCOL %% nblocks == 0){
  ## split column numbers into 'nblocks' groups
#  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

  SPLIT <- split(1:NCOL, ceiling(seq_along(1:NCOL)/(NCOL%/%nblocks)))
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
#  }else{
#    NCOLNEW = NCOL - (NCOL%%nblocks + NCOL%/%nblocks)
#    SPLIT <- split(1:NCOLNEW, rep(1:(nblocks-1), each = NCOLNEW/(nblocks-1)))
#  }
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  `%dopar%` <- foreach::`%dopar%`
  results <- foreach:::foreach(i = 1:nrow(COMBS), .combine='rbind')%dopar%{
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]] - 1
        G2 <- SPLIT[[COMB[2]]] - 1
        #if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        G1M = getColfromStdGenoMultiMarkersMat(G1)
        G2M = getColfromStdGenoMultiMarkersMat(G2)
#	cat("G1M ", dim(G1M), "\n")
#	cat("G2M ", dim(G2M), "\n")
        GRM <- t(G1M)%*%(G2M)
#	cat("G1M: ", G1M, "\n")
#	cat("G2M: ", G2M, "\n")
#	cat("GRM: ", GRM, "\n")
        if(sum(G1 != G2) == 0){
                GRM[lower.tri(GRM, diag = FALSE)] = 0
        }
        GRM <- GRM/NROW
        indice <- which(GRM >= relatednessCutoff, arr.ind=T)
#	cat("GRM[1:10,1:10]: ", GRM[1:10,1:10], "\n")
#	cat("relatednessCutoff: ", relatednessCutoff, "\n")
#	cat("indice: ", indice, "\n")
        cbind(G1[indice[,1]], G2[indice[,2]])
#       resultsIVec =c(resultsIVec, G1[indice[,1]])
#        resultsJVec =c(resultsJVec, G2[indice[,2]])
        #corMAT[G1, G2] <- COR
        #corMAT[G2, G1] <- t(COR)
#       GRM <- NULL
#       indice = NULL
}
  gc()
  return(results)
}


#refineKinPar = function(iMat, relatednessCutoff, W, tauVecNew, nblocks = 10, verbose = TRUE, ncore= 1){
refineKinPar = function(relatednessCutoff, W, tauVecNew, nblocks = 10, verbose = TRUE, ncore= 1){

	#chunk iMat
  doParallel:::registerDoParallel(ncore)
  NROW <- nrow(iMat)
  SPLIT <- split(1:NROW, ceiling(seq_along(1:NROW)/(NROW%/%nblocks)))
  GRMvec = rep(0, NROW)
  print(length(GRMvec))
  `%dopar%` <- foreach::`%dopar%`
  mMarkers = gettotalMarker()
  for(j in 1:mMarkers){
    cat("j is ", j, "\n")
    stdGeno = Get_OneSNP_StdGeno(j-1)
    results <- foreach:::foreach(i = 1:length(SPLIT))%dopar%{
        G1 <- SPLIT[[i]]
	#print(G1)
        if (verbose) cat("Block", i , "\n")
        flush.console()
	for(m in G1){
		print(m)
 		print("OK2")
		print(iMat[m,1])	
		print(iMat[m,2])	
		print(stdGeno[iMat[m,1]])	
		print(stdGeno[iMat[m,2]])	
		print(GRMvec[m])	
		GRMvec[m] = GRMvec[m] + stdGeno[iMat[m,1]]*stdGeno[iMat[m,2]]/mMarkers
		print(GRMvec[m])	
		print("OK3")	
	}
  } 
  print("OK4")
 }
 print("OK1")
 return(GRMvec)		

}


#createSparseKinParallel = function(markerIndexVec, nblocks, ncore, relatednessCutoff, W, tauVecNew){
createSparseKinParallel = function(nblocks, ncore, relatednessCutoff){
  #get MAT
  #MAT = Get_MultiMarkersBySample_StdGeno_Mat(markerIndexVec)  
  #MAT = Get_MultiMarkersBySample_StdGeno_Mat()  
  setRelatednessCutoff(relatednessCutoff)
  Get_MultiMarkersBySample_StdGeno_Mat()  
#  cat("dim(MAT) is ", dim(MAT), "\n")
  tp0 = proc.time()
#  indexVec = bigGRMPar(MAT, nblocks = nblocks, verbose = FALSE, ncore = nblocks, relatednessCutoff = relatednessCutoff)
#  indexVec = bigGRMPar_new(nblocks = nblocks, verbose = TRUE, ncore = nblocks, relatednessCutoff = relatednessCutoff)

  printComb(3)
  #indexVec = findIndiceRelatedSample()
  findIndiceRelatedSample()

  #print(indexVec)
  tp1 = proc.time()
  cat("tp1 - tp0: ", tp1-tp0, "\n")
#  cat(indexVec)
  #sparseKinList = refineKin(indexVec-1, relatednessCutoff, W, tauVecNew)
  sparseKinList = refineKin(relatednessCutoff)

#  sparseKinList$kinValue = sparseKinList$kinValue * tauVecNew[2]

  Nval = getNnomissingOut()
	
  sparseKinList$iIndex = c(sparseKinList$iIndex, seq(1:Nval))
  sparseKinList$jIndex = c(sparseKinList$jIndex, seq(1:Nval))
#  diagKin = getDiagOfSigma(W, tauVecNew)
  diagKin = rep(1, Nval)
  sparseKinList$kinValue = c(sparseKinList$kinValue, diagKin)
#  rm(diagKin)

  #sparseKinList = refineKin(indexVec, relatednessCutoff, W, tauVecNew)
  #GRMvec = refineKinPar(indexVec, relatednessCutoff = relatednessCutoff, W = W, tauVecNew = tauVecNew, nblocks = nblocks, verbose = TRUE, ncore= nblocks) 
  #sparseKinList = shortenList(indexVec-1, GRMvec, relatednessCutoff, W, tauVecNew)
 tp2 = proc.time()
  cat("tp2 - tp1: ", tp2-tp1, "\n")
  return(sparseKinList)
}



getSparseSigma = function(outputPrefix="",
                sparseGRMFile=NULL,
                sparseGRMSampleIDFile="",
                numRandomMarkerforSparseKin = 500,
                relatednessCutoff = 0.125,
		obj.glmm.null,
                W, tauVecNew){

  cat("sparse GRM will be used\n")
#  sparseGRMFile = paste0(outputPrefix, ".sparseGRM.mtx")
  if(is.null(sparseGRMFile)){
    freqVec = getAlleleFreqVec()
    MAFindex = which(freqVec >= 0.01 & freqVec <= 0.99)
    cat(numRandomMarkerforSparseKin, "genetic markers are randomly selected to decide which samples are related\n")
    if(length(MAFindex) < numRandomMarkerforSparseKin){
      stop("ERROR! not enough genetic markers with MAC >= 1% to detect which samples are related\n","Try include at least ", numRandomMarkerforSparseKin, " genetic markers with MAC >= 1% in the plink file\n")
    }

    markerIndexforSparseM = sample(MAFindex, size = numRandomMarkerforSparseKin, replace=FALSE)

    cat("Start detecting related samples for the sparse GRM\n")
    ta = proc.time()
    setSubMarkerIndex(markerIndexforSparseM -1)
    tb = proc.time()
    cat("tb-ta\n")
    print(tb-ta)


    cat("Start creating sparse GRM\n")
    ta = proc.time()
    sparseMList = createSparseKinParallel(nblocks = nThreads, ncore = nThreads, relatednessCutoff)
    tb = proc.time()
    cat("tb-ta\n")
    print(tb-ta)



    cat("length(sparseMList$iIndex): ", length(sparseMList$iIndex), "\n")
    print(sparseMList$iIndex[1:102])
    cat("length(sparseMList$jIndex): ", length(sparseMList$jIndex), "\n")
    print(sparseMList$jIndex[1:102])
    cat("length(sparseMList$kinValue): ", length(sparseMList$kinValue), "\n")
    print(sparseMList$kinValue[1:102])
    sparseGRM = Matrix:::sparseMatrix(i = as.vector(sparseMList$iIndex), j = as.vector(sparseMList$jIndex), x = as.vector(sparseMList$kinValue), symmetric = TRUE)
    cat("nrow(sparseGRM): ", nrow(sparseGRM), "\n")
    cat("ncol(sparseGRM): ", ncol(sparseGRM), "\n")
    cat("ncol(sparseGRM): ", sum(sparseGRM != 0), "\n")

    tc = proc.time()
    cat("tc-tb\n")
    print(tc-tb)

#    cat("td-tc\n")
#    print(td-tc)
    #cat("OK3", "\n")
  }else{ # if(sparseGRMFile=="")

       cat("sparse GRM has been specified\n")
       cat("read in sparse GRM from ",sparseGRMFile,"\n")

    sparseGRMLarge = Matrix:::readMM(sparseGRMFile)
    #cat("sparseSigmaFile: ", sparseSigmaFile, "\n")
    if(sparseGRMSampleIDFile != ""){
      if(!file.exists(sparseGRMSampleIDFile)){
        stop("ERROR! sparseSigmaSampleIDFile ", sparseGRMSampleIDFile, " does not exsit\n")
      }else{
        sparseGRMSampleID = data.frame(data.table:::fread(sparseGRMSampleIDFile, header=F, stringsAsFactors=FALSE))
        colnames(sparseGRMSampleID) = c("sampleID")
        sparseGRMSampleID$IndexGRM = seq(1,nrow(sparseGRMSampleID), by=1)
        sampleInModel = NULL
        sampleInModel$IID = obj.glmm.null$sampleID
        sampleInModel = data.frame(sampleInModel)
        sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
        cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
        mergeID = merge(sampleInModel, sparseGRMSampleID, by.x="IID", by.y = "sampleID")
        mergeID = mergeID[with(mergeID, order(IndexInModel)), ]
        indexIDofGRM=mergeID$IndexGRM
        #cat("Subset sparse GRM to be ", indexIDofSigma," by ", indexIDofSigma, "\n")
        sparseGRM = sparseGRMLarge[indexIDofGRM, indexIDofGRM]
        rm(sparseGRMLarge)
      }
    }else{#end of if(sparseSigmaSampleIDFile != "")
      stop("ERROR! sparseSigmaSampleIDFile is not specified\n")
    }

  #cat("sparse GRM has been specified\n")
  #cat("read in sparse GRM from ",sparseSigmaOutFile,"\n")
  #sparseSigma = Matrix:::readMM(sparseSigmaOutFile)
 }
  sparseGRMFile = paste0(outputPrefix,"_relatednessCutoff_",relatednessCutoff, ".sparseGRM.mtx")
  cat("write sparse GRM to ", sparseGRMFile ,"\n")
  Matrix:::writeMM(sparseGRM, sparseGRMFile)
  Nval = length(W)


  sparseSigma = sparseGRM * tauVecNew[2]
  diag(sparseSigma) = getDiagOfSigma(W, tauVecNew)

  sparseSigmaFile = paste0(outputPrefix, "_relatednessCutoff_",relatednessCutoff, ".sparseSigma.mtx")
  cat("write sparse Sigma to ", sparseSigmaFile ,"\n")
  Matrix:::writeMM(sparseSigma, sparseSigmaFile)

#    td = proc.time()
#    cat("td-tc\n")
#    print(td-tc)


  return(sparseSigma)
}

