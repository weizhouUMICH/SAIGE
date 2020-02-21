##saddlepoint approxmation for sum of weighted Poisson distribution
Korg_Poi<-function(t, mu, g)
{
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-mu*(exp(g*t1) - g*t1 - 1)
    out[i]<-sum(temp)
  }
  return(out)
}

K1_Poi<-function(t, mu, g)
{
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-mu * g * exp(g*t1) - mu * g
    out[i]<-sum(temp)
  }
  return(out)
}


K1_adj_Poi<-function(t, mu, g, q)
{
  n.t<-length(t)	
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-mu * g * exp(g*t1) - mu * g
    out[i]<-sum(temp)-q
  }
  return(out)
}


K2_Poi<-function(t, mu, g)
{
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-mu * g^2 * exp(g*t1)
    out[i]<-sum(temp, na.rm=TRUE)
  }
  return(out)
}


getroot_K1_uniroot_Poi = function(init,mu,g,q,m1,tol=.Machine$double.eps^0.25,maxiter=1000){

	newroot = uniroot(function(x) K1_adj_Poi(x,mu,g,q), c(-1000,1000),extendInt="yes", check.conv=TRUE, maxiter=maxiter) 
	if(newroot$iter==maxiter)
        {
        	conv<-FALSE
        }else{
		conv<-TRUE
	}
	
	return(list(root=newroot$root,n.iter=newroot$iter,Is.converge=conv))
}

getroot_K1_Poi<-function(init,mu,g,q,m1,tol=.Machine$double.eps^0.25,maxiter=1000)
{
    t<-init
    K1_eval<-K1_adj_Poi(t,mu,g,q)
    prevJump<- Inf
    rep<-1
    repeat
    {
      K2_eval<-K2_Poi(t,mu,g)
      tnew<-t-K1_eval/K2_eval
      if(is.na(tnew))
      {
        conv=FALSE
        break
      }
      if(abs(tnew-t)<tol)
      {
        conv<-TRUE
        break
      }
      if(rep==maxiter)
      {
        conv<-FALSE
        break
      }
      
      newK1<-K1_adj_Poi(tnew,mu,g,q)
      if(sign(K1_eval)!=sign(newK1))
      {
        if(abs(tnew-t)>prevJump-tol)
        {
          tnew<-t+sign(newK1-K1_eval)*prevJump/2
          newK1<-K1_adj_Poi(tnew,mu,g,q)
          prevJump<-prevJump/2
        } else {
          prevJump<-abs(tnew-t)
        }
      }
      
      rep<-rep+1
      t<-tnew
      K1_eval<-newK1
    } 
    return(list(root=t,n.iter=rep,Is.converge=conv))
 # }
}


Get_Saddle_Prob_Poi<-function(zeta, mu, g, q) 
{
  k1<-Korg_Poi(zeta, mu, g)
  #cat("k1 is ", k1, "\n")
  k2<-K2_Poi(zeta, mu, g)
  #cat("k2 is ", k2, "\n")
  if(is.finite(k1) && is.finite(k2))
  {
    temp1<-zeta * q - k1
    
    
    w<-sign(zeta) * (2 *temp1)^{1/2}
    v<- zeta * (k2)^{1/2}
    
    Z.test<-w + 1/w * log(v/w)	
    
    if(Z.test > 0){
      pval<-pnorm(Z.test, lower.tail = FALSE)
    } else {
      pval= -pnorm(Z.test, lower.tail = TRUE)
    }	
  } else {
    pval<-0
  }
  
  return(pval)
}
