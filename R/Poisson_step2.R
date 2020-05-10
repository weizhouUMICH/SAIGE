scoreTest_SAIGE_survivalTrait_cond_sparseSigma=function(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y,varRatio, Cutoff, rowHeader, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }
  Run1=TRUE
    G0 = matrix(G0, ncol = 1)
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
    g = G
    G0_mc = matrix(G0-mean(G0), ncol = 1)
    XVG0_mc = eigenMapMatMult(obj.noK$XV, G0_mc)
    g_mc = G0_mc - eigenMapMatMult(obj.noK$XXVX_inv, XVG0_mc)	

if(!isCondition){
  if(IsSparse==TRUE){
    #if(MAF < 0.05){
    #   out.score<-Score_Test_Sparse(obj.noK, G0, mu.a, mu2.a, varRatio );
    #}else{
       #out.score<-Score_Test(obj.noK, G0,mu.a, mu2.a, varRatio );
       out.score<-Score_Test_Survival(obj.noK, g,  mu.a, mu2.a, varRatio );
    #}
#    if(out.score["pval.noadj"] > 0.05){
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

    NAset = which(G0==0)
    out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma(g, g_mc, AC2, AC,NAset, y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv)

    if(isCondition){
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"], Tstat_c = out1["Tstat_c"], p.value.c = out1["p.value.c"], var1_c = out1["var1_c"], BETA_c = out1["BETA_c"], SE_c = out1["SE_c"])

    }else{
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"])
     #outVec = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
   }
  }
  return(outVec)
}




Saddle_Prob_Poisson=function (Score, pval.noadj, mu, g, Cutoff = 2, alpha = 5*10^-8, m1, var1){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL
    #cat("Score is ", Score, "\n")
    #print(g[1:20])


    #Score <- q - m1
    #qinv = -sign(q - m1) * abs(q - m1) + m1
    #pval.noadj <- pchisq((q - m1)^2/var1, lower.tail = FALSE,
    #    df = 1)
    Is.converge = TRUE

    #if (abs(q - m1)/sqrt(var1) < Cutoff) {
    #    pval = pval.noadj
    #}else {
#        print("Saddle_Prob_Poisson")
        out.uni1 <- getroot_K1_Poi(0, mu = mu, g = g, q = Score)
        out.uni2 <- getroot_K1_Poi(0, mu = mu, g = g, q = (-1)*Score)
        if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
	   p1 <- tryCatch(Get_Saddle_Prob_Poi(out.uni1$root, mu, g, q=Score), error=function(e) {return(pval.noadj/2)})	
	   p2 <- tryCatch(Get_Saddle_Prob_Poi(out.uni2$root, mu, g, q = (-1)*Score), error=function(e) {return(pval.noadj/2)})	
            #p1 <- Get_Saddle_Prob_Poi(out.uni1$root, mu, g, q)
            #p2 <- Get_Saddle_Prob_Poi(out.uni2$root, mu, g, qinv)
#	    cat("p1 p2: ", p1, " ", p2, "\n")	

            pval = abs(p1) + abs(p2)
            Is.converge = TRUE
        }
        else {
            print("Error_Converge")
            pval <- pval.noadj
            Is.converge = FALSE
        }
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
            Is.converge = Is.converge, Score = Score))
}



