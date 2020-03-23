scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast=function(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y,varRatio, Cutoff, rowHeader, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }
  Run1=TRUE

if(!isCondition){
  if(IsSparse==TRUE){
    if(MAF < 0.05){
       out.score<-Score_Test_Sparse(obj.noK, G0, mu.a, mu2.a, varRatio );
    }else{
       out.score<-Score_Test(obj.noK, G0,mu.a, mu2.a, varRatio );
    }
    #if(out.score["pval.noadj"] > 0.05){
    if(abs(as.numeric(unlist(out.score["Tstat"])[1])/sqrt(as.numeric(unlist(out.score["var1"])[1]))) < Cutoff){
       if(AF > 0.5){
         out.score$BETA = (-1)*out.score$BETA
         out.score$Tstat = (-1)*out.score$Tstat
       }
       outVec = list(BETA = out.score$BETA, SE = out.score$SE, Tstat = out.score$Tstat, p.value = out.score$pval.noadj, p.value.NA = out.score$pval.noadj, Is.converge = 1, var1 = out.score$var1, var2 = out.score$var2)
       Run1=FALSE
     }
  }
}

  #cat("Run1: ", Run1, "\n")
  if(Run1){
    G0 = matrix(G0, ncol = 1)
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
    g = G
    NAset = which(G0==0)
    #cat("length(NAset): ", length(NAset), "\n")	
    #if(length(NAset)/length(G)<0.5){
    #out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma(g, AC2, AC,NAset, y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv)
    #}else{
    out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast(g, AC2, AC,NAset, y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv)
    #print(out1)
    #}

    if(isCondition){
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"], Tstat_c = out1["Tstat_c"], p.value.c = out1["p.value.c"], var1_c = out1["var1_c"], BETA_c = out1["BETA_c"], SE_c = out1["SE_c"])

    }else{
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"])
     #outVec = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
   }
  }
  return(outVec)
}



scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast=function(g, AC, AC_true, NAset, y, mu, varRatio, Cutoff, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

  #g = G/sqrt(AC)
  q = innerProduct(g, y)
  m1 = innerProduct(g, mu)
  Tstat = q-m1
  var2 = innerProduct(mu, g*g)
  var1 = var2 * varRatio

  if(!is.null(sparseSigma)){
    #pcginvSigma<-pcg(sparseSigma, g)
    pcginvSigma<-solve(sparseSigma, g, sparse=T)
    var2b = as.matrix(t(g) %*% pcginvSigma)
    var1 = var2b * varRatio
  }

  if(isCondition){
    T2stat = OUT_cond[,2]
    G1tilde_P_G2tilde = matrix(G1tilde_P_G2tilde,nrow=1)
    Tstat_c = Tstat - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% T2stat
    var1_c = var1 - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% t(G1tilde_P_G2tilde)
  }

  AF = AC_true/(2*length(y))
  if(AF > 0.5){
    Tstat = (-1)*Tstat
    if(isCondition){
      Tstat_c = (-1)*Tstat_c
    }
  }

  qtilde = Tstat/sqrt(var1) * sqrt(var2) + m1

  if(length(NAset)/length(g) < 0.5){
	#print("Saddle_Prob_Poisson")
    out1 = Saddle_Prob_Poisson(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
  }else{
	#print("Saddle_Prob_Poisson_fast")
    out1 = Saddle_Prob_Poisson_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8)
  }

  #out1 = c(out1, var1 = var1)
  #out1 = c(out1, var2 = var2)
  out1$var1 = var1
  out1$var2 = var2

  #01-27-2019
  #as g is not divided by sqrt(AC), the sqrt(AC) is removed from the denominator
  #logOR = (Tstat/var1)/sqrt(AC)
  logOR = Tstat/var1
  SE = abs(logOR/qnorm(out1$p.value/2))
#  out1 = c(out1, BETA = logOR, SE = SE, Tstat = Tstat)
  out1$BETA=logOR
  out1$SE=SE
  out1$Tstat = Tstat

  if(isCondition){
    if(var1_c <= (.Machine$double.xmin)^2){
      out1 = c(out1, var1_c = var1_c,BETA_c = NA, SE_c = NA, Tstat_c = Tstat_c, p.value.c = 1, p.value.NA.c = 1)
    }else{

      qtilde_c = Tstat_c/sqrt(var1_c) * sqrt(var2) + m1
      if(length(NAset)/length(g) < 0.5){
        out1_c = Saddle_Prob_Poisson(q=qtilde_c, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
      }else{
        out1_c = Saddle_Prob_Poisson_fast(q=qtilde_c,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8)
      #  out1_c = SPAtest:::Saddle_Prob_fast(q=qtilde_c,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="p")
      }
    #01-27-2019
    #logOR_c = (Tstat_c/var1_c)/sqrt(AC)
    logOR_c = Tstat_c/var1_c
    SE_c = abs(logOR_c/qnorm(out1_c$p.value/2))
    out1 = c(out1, var1_c = var1_c,BETA_c = logOR_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = out1_c$p.value, p.value.NA.c = out1_c$p.value.NA)
    }

  }

  #print("out1")
  #print(out1)

  return(out1)
}


Saddle_Prob_Poisson_fast=function (q, mu, g, gNA,gNB,muNA,muNB, Cutoff = 2, alpha = 5*10^-8){
    m1 <- sum(mu * g)
    var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL

    #NAmu= m1-sum(gNB*muNB)
    NAsigma=var1-sum(muNB*gNB^2)

    #NAsigma = sum(muNA*gNA^2)
    Score <- q - m1
    #qinv = -sign(q - m1) * abs(q - m1) + m1
    pval.noadj <- pchisq((q - m1)^2/var1, lower.tail = FALSE,
        df = 1)
    Is.converge = TRUE

    if (abs(q - m1)/sqrt(var1) < Cutoff) {
        pval = pval.noadj
    }else {
	#print("Saddle_Prob_Poisson_fast >= Cutoff")

	#Korg_Poi_result = Korg_Poi(t=0.1, mu, g)
	#Korg_Poi_fast_result = Korg_Poi_fast(t=0.1, mu, g, gNA,gNB,muNA,muNB,NAmu,NAsigma)	
	#cat("Korg_Poi_result: ", Korg_Poi_result, "\n")
	#cat("Korg_Poi_fast_result: ", Korg_Poi_fast_result, "\n")
        out.uni1 <- getroot_K1_Poi_fast(0, mu = mu, g = g, q = Score, gNA=gNA,gNB=gNB,muNA=muNA,muNB=muNB,NAsigma=NAsigma)
        out.uni2 <- getroot_K1_Poi_fast(0, mu = mu, g = g, q = (-1)*Score, gNA=gNA,gNB=gNB,muNA=muNA,muNB=muNB,NAsigma=NAsigma)
         #cat("out.uni1 out.uni2: ", out.uni1$root, " ", out.uni2$root, "\n")

        if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
		p1<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni1$root, mu, g, q=Score,gNA,gNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
		p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu, g, q=(-1)*Score,gNA,gNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})	
		#cat("p1 p2: ", p1, " ", p2, "\n")

            pval = abs(p1) + abs(p2)
            Is.converge = TRUE
        }else {
            print("Error_Converge")
            pval <- pval.noadj
            Is.converge = FALSE
        }
    }
    return(list(p.value = pval, p.value.NA = pval.noadj,
            Is.converge = Is.converge, Score = Score))
}



