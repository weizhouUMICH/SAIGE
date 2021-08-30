##################################################################
## HELPER FUNCTION:  fast.logistf.control
## This function provides the convergence control parameters
##################################################################
fast.logistf.control <- function (maxit = 50, maxhs = 15, maxstep = 15, lconv = 1e-05, 
    gconv = 1e-05, xconv = 1e-05) 
{
    list(maxit = maxit, maxhs = maxhs, maxstep = maxstep, lconv = lconv, 
        gconv = gconv, xconv = xconv)
}


##################################################################
## HELPER FUNCTION:  fast.logDet
##################################################################
fast.logDet <- function (x) {
    my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (2 * sum(log(diag(my.chol))))
    }
}   

##################################################################
## HELPER FUNCTION:  fast.invFisher
##################################################################
fast.invFisher <- function(x) {
  my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (chol2inv(my.chol))
    }
  #ifelse(is.na(my.chol), NA, chol2inv(my.chol))
}


fast.poisf.fit <- function (x, y, weight = NULL, Lambda = NULL, offset = NULL, firth = TRUE, col.fit = NULL, 
    init = NULL, control) {
    n <- nrow(x)
    k <- ncol(x)
    if (is.null(Lambda)) 		#These will be baseline hazards Lambda0
        Lambda = rep(1, n)
    if (is.null(init)) 
        init = rep(0, k)
    if (is.null(col.fit)) 
        col.fit = 1:k
    if (is.null(offset)) 
        offset = rep(0, n)
    if (is.null(weight)) 
        weight = rep(1, n)
    if (col.fit[1] == 0) 
        maxit <- 0
    if (missing(control)) 
        control <- fast.logistf.control()
    maxit <- control$maxit
    maxstep <- control$maxstep
    maxhs <- control$maxhs
    lconv <- control$lconv
    gconv <- control$gconv
    xconv <- control$xconv
    beta <- init
    iter <- 0
    pi <- as.vector(Lambda*exp(x %*% beta + offset))	#Changed it to Poisson mean with Lambda(t) multiplied
    evals <- 1
    repeat {
        beta.old <- beta
        XW2 <- t(x * (weight * pi)^0.5)		#Changed to Poisson variance
        myQR <- qr(t(XW2))
        Q <- qr.Q(myQR)
        h <- (Q*Q) %*% rep(1, ncol(Q))        
        if (firth) 
            U.star <- crossprod(x, weight * (y - pi) + h * 0.5)	#Changed this to reflect derivative of log|I(beta)|
        else U.star <- crossprod(x, weight * (y - pi))
        XX.covs <- matrix(0, k, k)
        if (col.fit[1] != 0) {
            XX.XW2 <- t(x[, col.fit, drop=FALSE] * (weight * pi)^0.5)	#Changed to reflect Poisson variance
            XX.Fisher <- crossprod(t(XX.XW2))
            XX.covs[col.fit, col.fit] <- fast.invFisher(XX.Fisher)   
        }
        if(all(is.na(XX.covs)) == T) {
            break
        }  
        delta <- as.vector(XX.covs %*% U.star)
        delta[is.na(delta)] <- 0
        mx <- max(abs(delta))/maxstep
        if (mx > 1) 
            delta <- delta/mx
        evals <- evals + 1
        if (maxit > 0) {
            iter <- iter + 1
            beta <- beta + delta
                pi <- as.vector(Lambda * exp(x %*% beta + offset))		#Changed this

        }
        if (iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star[col.fit]) < 
            gconv)))) 
            break
    }
    # Error catching (if chol(x) not positive definite)
    if(all(is.na(XX.covs))==T) {
        var <- XX.covs			#This is not supposed to be the final variance of the estimates.
        list(beta = NA, var = var, pi = NA, hat.diag = NA, 
  iter = NA, evals = NA, conv = c(NA, 
            NA, NA))
    } else {
        var <- XX.covs
        list(beta = beta, var = var, pi = pi, hat.diag = h, 
            iter = iter, evals = evals, conv = c(max(abs(U.star)),
		 max(abs(delta))))
    }
}
