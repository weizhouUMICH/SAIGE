##saddlepoint approxmation for sum of weighted Poisson distribution
Korg_Poi_fast<-function(t, mu, g, gNA,gNB,muNA,muNB,NAsigma)
{
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-muNB*(exp(gNB*t1) - gNB*t1 - 1)
    #out[i]<-sum(temp)+NAmu*t1+0.5*NAsigma*t1^2
    out[i]<-sum(temp)+0.5*NAsigma*t1^2 
 }
  return(out)
}

#K1_Poi<-function(t, mu, g)
#{
#  n.t<-length(t)
#  out<-rep(0,n.t)
  
#  for(i in 1:n.t){
#    t1<-t[i]
#    temp<-mu * g * exp(g*t1) - mu * g
#    out[i]<-sum(temp)
#  }
#  return(out)
#}


K1_adj_Poi_fast<-function(t, mu, g, q, gNA,gNB,muNA,muNB,NAsigma)
{
  n.t<-length(t)	
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-muNB * gNB * exp(gNB*t1) - muNB * gNB
    #temp2<-NAmu+NAsigma*t1
    temp2<-NAsigma*t1
    out[i]<-sum(temp)-q + temp2
  }
  return(out)
}


K2_Poi_fast<-function(t, mu, g, gNA,gNB,muNA,muNB,NAsigma)
{
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-muNB * gNB^2 * exp(gNB*t1)
    out[i]<-sum(temp, na.rm=TRUE) + NAsigma
  }
  return(out)
}


getroot_K1_uniroot_Poi_fast = function(init,mu,g,q,m1,gNA,gNB,muNA,muNB,NAsigma,tol=.Machine$double.eps^0.25,maxiter=1000){

	newroot = uniroot(function(x) K1_adj_Poi_fast(x,mu,g,q,gNA,gNB,muNA,muNB,NAsigma), c(-1000,1000),extendInt="yes", check.conv=TRUE, maxiter=maxiter) 
	if(newroot$iter==maxiter)
        {
        	conv<-FALSE
        }else{
		conv<-TRUE
	}
	
	return(list(root=newroot$root,n.iter=newroot$iter,Is.converge=conv))
}

getroot_K1_Poi_fast<-function(init,mu,g,q,m1,gNA,gNB,muNA,muNB,NAsigma,tol=.Machine$double.eps^0.25,maxiter=1000)
{
    t<-init
    K1_eval<-K1_adj_Poi_fast(t,mu,g,q,gNA,gNB,muNA,muNB,NAsigma)
    #cat("K1_eval: ", K1_eval, "\n")
    prevJump<- Inf
    rep<-1
    repeat
    {
      K2_eval<-K2_Poi_fast(t,mu,g,gNA,gNB,muNA,muNB,NAsigma)
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
      
      newK1<-K1_adj_Poi_fast(tnew,mu,g,q,gNA,gNB,muNA,muNB,NAsigma)
      if(sign(K1_eval)!=sign(newK1))
      {
        if(abs(tnew-t)>prevJump-tol)
        {
          tnew<-t+sign(newK1-K1_eval)*prevJump/2
          newK1<-K1_adj_Poi_fast(tnew,mu,g,q,gNA,gNB,muNA,muNB,NAsigma)
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


Get_Saddle_Prob_Poi_fast<-function(zeta, mu, g, q,gNA,gNB,muNA,muNB,NAsigma) 
{
  k1<-Korg_Poi_fast(zeta, mu, g,gNA,gNB,muNA,muNB,NAsigma)
  #cat("k1 is ", k1, "\n")
  k2<-K2_Poi_fast(zeta, mu, g,gNA,gNB,muNA,muNB,NAsigma)
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
