scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast=function(resq, q0, Wq, qW1, XWq, G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y,varRatio, Cutoff, rowHeader, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }
  Run1=TRUE
  idx_no0<-which(G0>0)
  G0 = matrix(G0, ncol = 1)
  tp0 = proc.time()
  mug = mean(G0)
  tp1 = proc.time()
  print("tp1-tp0")
  print(tp1-tp0)
  ##########old way to get g_mc
  #G0_mc = matrix(G0-mean(G0), ncol = 1)
  #XVG0_mc = eigenMapMatMult(obj.noK$XV, G0_mc)
  #g_mc = G0_mc - eigenMapMatMult(obj.noK$XXVX_inv, XVG0_mc)
  #########
if(!isCondition){
  if(IsSparse==TRUE){
    if(MAF < 0.05){
       out.score<-Score_Test_Sparse_Survival(obj.noK, G0, mug, mu.a, mu2.a, varRatio, resq, Wq, qW1, XWq);
        tp2a = proc.time()
	print("tp2a-tp1")
 	print(tp2a-tp1)
    }else{
       out.score<-Score_Test_Survival(obj.noK, G0, mug, mu.a, mu2.a, varRatio, resq, Wq, qW1, XWq);
        tp2b = proc.time()
	print("tp2b-tp1")
  	print(tp2b-tp1)
    }
    #cat("out.score$Tstat: ", out.score$Tstat, "\n")
    #cat("out.score1$Tstat: ", out.score1$Tstat, "\n")
    #cat("out.score$var2: ", out.score$var2, "\n")
    #cat("out.score1$var2: ", out.score1$var2, "\n")
    ##if(out.score["pval.noadj"] > 0.05){
    if(abs(as.numeric(unlist(out.score["Tstat"])[1])/sqrt(as.numeric(unlist(out.score["var1"])[1]))) < Cutoff){
       if(AF > 0.5){
         out.score$BETA = (-1)*out.score$BETA
         out.score$Tstat = (-1)*out.score$Tstat
       }
       outVec = list(BETA = out.score$BETA, SE = out.score$SE, Tstat = out.score$Tstat, p.value = out.score$pval.noadj, p.value.NA = out.score$pval.noadj, Is.converge = 1, var1 = out.score$var1, var2 = out.score$var2)
       Run1=FALSE
   }else{
       if(MAF < 0.05){

	if(FALSE){
         X1 = obj.noK$X1
       ##cat("dim(X1) ", dim(X1), "\n")
       ##cat("dim(out.score$Z) ", dim(out.score$Z), "\n")
         g = rep(0, nrow(X1))
         g[idx_no0] = out.score$g_tilde1 + (out.score$B)
         g = g - X1%*%(out.score$Z)
         g = g - mug*q0
	}
         XVG0 = eigenMapMatMult(obj.noK$XV, G0)
         g = G0 - eigenMapMatMult(obj.noK$XXVX_inv, XVG0)

        tp3a = proc.time()
	print("tp3a-tp2a")
  	print(tp3a-tp2a)
      #print(g[1:20])
      #print(g_mc[1:20])
       }else{
        g = out.score$g_tilde	
        tp3b = proc.time()
	print("tp3b-tp2b")
  	print(tp3b-tp2b)
      }
   }
 }
}

  #cat("Run1: ", Run1, "\n")
  if(Run1){
    #G0 = matrix(G0, ncol = 1)
    #XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    #G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
    #g = G
    NAset = which(G0==0)
    tp4 = proc.time()
    out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast(g, Score = out.score$Tstat, pval.noadj = out.score$pval.noadj, var1_a = out.score$var1, var2_a = out.score$var2, Wq, qW1,mug, AC2, AC,NAset, y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv)
    tp5 = proc.time()
    print("tp5-tp4")
        print(tp5-tp4)	
    if(isCondition){
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"], Tstat_c = out1["Tstat_c"], p.value.c = out1["p.value.c"], var1_c = out1["var1_c"], BETA_c = out1["BETA_c"], SE_c = out1["SE_c"])

    }else{
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"])
     #outVec = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
   }
 }
  #cat("p.value: ", as.numeric(outVec$p.value), "\n")
  #cat("p.value.NA: ", as.numeric(outVec$p.value.NA), "\n")
  return(outVec)
}



scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast=function(g, Score, pval.noadj, var1_a, var2_a, Wq, qW1, mug, AC, AC_true, NAset, y, mu, varRatio, Cutoff, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

  #g = G/sqrt(AC)
  #q = innerProduct(g, y)
  m1 = innerProduct(g, mu)
  #Tstat = q-m1
  #Tstat = Score
  #var2 = innerProduct(mu, g*g)
  #var2c_old = innerProduct(mu, g_mc*g_mc)
  #var2c = var2 - 2*mug*innerProduct(g,Wq) + mug^2*qW1
  
  #cat("var2c_old ", var2c_old, "\n")
  #cat("var2c ", var2c, "\n")
  #var1 = var2c * varRatio
  var1 = var1_a
  var2 = var2_a

  if(!is.null(sparseSigma)){
    #pcginvSigma<-pcg(sparseSigma, g)
    pcginvSigma<-solve(sparseSigma, g, sparse=T)
    var2b = as.matrix(t(g) %*% pcginvSigma)
    var1 = var2b * varRatio
  }

  if(isCondition){
    T2stat = OUT_cond[,2]
    G1tilde_P_G2tilde = matrix(G1tilde_P_G2tilde,nrow=1)
    Tstat_c = Score - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% T2stat
    var1_c = var1 - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% t(G1tilde_P_G2tilde)
  }

  AF = AC_true/(2*length(y))
  #if(AF > 0.5){
  #  Tstat = (-1)*Tstat
  #  if(isCondition){
  #    Tstat_c = (-1)*Tstat_c
  #  }
  #}

  #qtilde = Tstat/sqrt(var1) * sqrt(var2) + m1
  Score2 = Score/sqrt(var1) * sqrt(var2)
  if(length(NAset)/length(g) < 0.5){
    #print("Saddle_Prob_Poisson")
    out1 = Saddle_Prob_Poisson(Score=Score2, pval.noadj=pval.noadj, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8, m1=m1, var1=var2)
  }else{
    #print("Saddle_Prob_Poisson_fast")
    out1 = Saddle_Prob_Poisson_fast(Score=Score2, pval.noadj=pval.noadj, g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, m1=m1, var1=var2)
  }

  out1$var1 = var1
  out1$var2 = var2

  logOR = Score/var1
  SE = abs(logOR/qnorm(out1$p.value/2))
  out1$BETA=logOR
  out1$SE=SE
  out1$Tstat = Score

  if(isCondition){
    if(var1_c <= (.Machine$double.xmin)^2){
      out1 = c(out1, var1_c = var1_c,BETA_c = NA, SE_c = NA, Tstat_c = Tstat_c, p.value.c = 1, p.value.NA.c = 1)
    }else{

      #qtilde_c = Tstat_c/sqrt(var1_c) * sqrt(var2) + m1
      pval.noadj_c<-pchisq((Tstat_c)^2/(var1_c), lower.tail = FALSE, df=1)
      if(length(NAset)/length(g) < 0.5){
######To improve
        out1_c = Saddle_Prob_Poisson(Score=Tstat_c, pval.noadj=pval.noadj_c, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8, m1=m1, var1=var1_c)
      }else{
######To improve
        out1_c = Saddle_Prob_Poisson_fast(Score=Tstat_c, pval.noadj=pval.noadj_c, g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, m1=m1, var1=var1_c)
      }
      logOR_c = Tstat_c/var1_c
      SE_c = abs(logOR_c/qnorm(out1_c$p.value/2))
      out1 = c(out1, var1_c = var1_c,BETA_c = logOR_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = out1_c$p.value, p.value.NA.c = out1_c$p.value.NA)
    }

  }

  return(out1)
}


Saddle_Prob_Poisson_fast=function (Score, pval.noadj, mu, g, gNA,gNB,muNA,muNB, Cutoff = 2, alpha = 5*10^-8, m1, var1){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL

    #NAmu= m1-sum(gNB*muNB)
    NAsigma=var1-sum(muNB*gNB^2)

    #cat("Score is ", Score, "\n")
    #cat("NAsigma is ", NAsigma, "\n")
    #print(mu[1:20])
    #print(g[1:20])

    #NAsigma = sum(muNA*gNA^2)
    #Score <- q - m1
    #qinv = -sign(q - m1) * abs(q - m1) + m1
    #pval.noadj <- pchisq((q - m1)^2/var1, lower.tail = FALSE,
    #    df = 1)
    #Is.converge = TRUE

    #if (abs(q - m1)/sqrt(var1) < Cutoff) {
    #    pval = pval.noadj
    #}else {
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
#		cat("p1 p2: ", p1, " ", p2, "\n")

            pval = abs(p1) + abs(p2)
            Is.converge = TRUE
        }else {
            print("Error_Converge")
            pval <- pval.noadj
            Is.converge = FALSE
        }
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
            Is.converge = Is.converge, Score = Score))
}



