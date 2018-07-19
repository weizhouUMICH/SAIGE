  if(traitType == "binary"){
    cat("It is a binary trait\n")
    resultHeader = c(dosageFilecolnamesSkip, "N", "BETA", "SE", "Tstat", "p.value", "p.value.NA", "Is.SPA.converge","varT","varTstar")

    if(IsOutputAFinCaseCtrl){
      resultHeader = c(resultHeader, "AF.Cases", "AF.Controls")
    }

    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
    if(Cutoff < 10^-2){
        Cutoff=10^-2
    }
 
    y = obj.glm.null$y
    OUT = NULL
    numPassMarker = 0
    NSparse=0
    mth = 0
    y1Index = which(y == 1)
    NCase = length(y1Index)
    y0Index = which(y == 0)
    NCtrl = length(y0Index)

    cat("Analyzing ", NCase, " cases and ",NCtrl, " controls \n")
   
    N = length(y)
    obj.noK$XVX_inv_XV = obj.noK$XXVX_inv * obj.noK$V
    indChromCheck = FALSE
    if(!obj.glmm.null$LOCO){
      mu = obj.glmm.null$fitted.values
      mu.a<-as.vector(mu)
      mu2.a<-mu.a *(1-mu.a)
      obj.noK$XVX = t(obj.noK$X1) %*% (obj.noK$X1 * mu2.a)      
      obj.noK$S_a = colSums(obj.noK$X1 * (y - mu.a))
    }else if(chrom != ""){
      chrom_v2 = as.character(chrom)
      chrom_v3 = as.numeric(gsub("[^0-9.]", "", chrom_v2))
      if(obj.glmm.null$LOCOResult[[chrom_v3]]$isLOCO){
        mu = obj.glmm.null$LOCOResult[[chrom_v3]]$fitted.values
        mu.a<-as.vector(mu)
        mu2.a<-mu.a *(1-mu.a)
      }else{
        mu = obj.glmm.null$fitted.values
        mu.a<-as.vector(mu)
        mu2.a<-mu.a *(1-mu.a)
      }
      obj.noK$XVX = t(obj.noK$X1) %*% (obj.noK$X1 * mu2.a)
      obj.noK$S_a = colSums(obj.noK$X1 * (y - mu.a))
    }else{
      cat("LOCO will be used, but chromosome for the dosage file is not specified. Will check each marker for its chromosome for LOCO!\n")
      indChromCheck = TRUE
    }
#####


  }else if(traitType == "quantitative"){
    cat("It is a quantitative trait\n")
    if(!isCondition){
    	resultHeader = c(dosageFilecolnamesSkip,  "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar")
    }else{
	resultHeader = c(dosageFilecolnamesSkip,  "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar","Tstat_cond", "p.value_cond", "varT_cond", "BETA_cond", "SE_cond" )
    }	


    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
    OUT = NULL
    numPassMarker = 0
    mth = 0
    sampleIndex = sampleIndex - 1
    y = obj.glm.null$y
    N = length(y)
    tauVec = obj.glmm.null$theta
    obj.noK$XVX = t(obj.noK$X1) %*% (obj.noK$X1)
    obj.noK$XVX_inv_XV = obj.noK$XXVX_inv * obj.noK$V

    indChromCheck = FALSE
    cat("obj.glmm.null$LOCO ", obj.glmm.null$LOCO, "\n")
   if(!obj.glmm.null$LOCO){
      mu = obj.glmm.null$fitted.values
      mu.a<-as.vector(mu)
      obj.noK$S_a = colSums(obj.noK$X1 * (y - mu.a))

    }else if(chrom != ""){
      chrom_v2 = as.character(chrom)
      chrom_v3 = as.numeric(gsub("[^0-9.]", "", chrom_v2))
      if(obj.glmm.null$LOCOResult[[chrom_v3]]$isLOCO){
        mu = obj.glmm.null$LOCOResult[[chrom_v3]]$fitted.values
        mu.a<-as.vector(mu)
      }else{
        mu = obj.glmm.null$fitted.values
        mu.a<-as.vector(mu)
      }
      obj.noK$S_a = colSums(obj.noK$X1 * (y - mu.a))      
#      if(!is.null(sparseSigma) & traitType == "binary" ){
#        sparseSigma = sparseSigma * diag(mu2.a)
#      }

    }else{
      cat("LOCO will be used, but chromosome for the dosage file is not specified. Will check each marker for its chromosome for LOCO!\n")
      indChromCheck = TRUE
    }

  }else{
    stop("ERROR! The type of the trait has to be either binary or quantitative\n")
  }


if(isCondition){
  OUT_cond = NULL
    for(i in 1:ncol(dosage_cond)){
      G0  = dosage_cond[,i]
      AC = sum(G0)
      N  = length(G0)
      AF = AC/(2*N)
      MAF = AF
      MAC = AC
      if(AF > 0.5){
        MAF = 1-AF
        MAC = 2*N - MAC
      }
      varRatio = getvarRatio(MAC, ratioVec)
    
      rowHeader = paste0("condMarker",i)

      if(traitType == "binary"){
	out1 = scoreTest_SAIGE_binaryTrait(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y, varRatio, Cutoff, rowHeader)
        OUT_cond = rbind(OUT_cond, c(as.numeric(out1[3]), as.numeric(out1[5]), as.numeric(out1[9])))
      }else if(traitType == "quantitative"){
        mu = obj.glmm.null$fitted.values
        mu.a<-as.vector(mu)
        obj.noK$S_a = colSums(obj.noK$X1 * (y - mu.a))

        out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.noK, AC, AF, y, mu, varRatio, tauVec, sparseSigma=sparseSigma)
        #out1 = scoreTest_SAIGE_quantitativeTrait(G0, obj.noK, AC, AF, y, mu, varRatio, tauVec)
        OUT_cond = rbind(OUT_cond, c(as.numeric(out1$BETA), as.numeric(out1$Tstat), as.numeric(out1$var1)))
      }
   OUT_cond = as.matrix(OUT_cond)

   } # end of for(i in 1:ncol(dosage_cond)){

  #covM for conditioning markers
  cat("dim(obj.noK$XVX_inv_XV): ", dim(obj.noK$XVX_inv_XV), "\n")
      cat("dim(dosage_cond): ", dim(dosage_cond), "\n")
      dosage_cond_tilde<-dosage_cond  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% dosage_cond)	
#      dosage_cond_tilde = dosage_cond - (obj.noK$XVX_inv_XV)%*%dosage_cond
      Mcond = ncol(dosage_cond_tilde)	
      #G0_tilde = G0 - (obj.noK$XVX_inv_XV)%*%G0
      #dosage_cond_tilde = cbind(G0_tilde, dosage_cond_tilde)
      covM = matrix(0,nrow=Mcond+1, ncol = Mcond+1)
      covMsub = getcovM(dosage_cond_tilde, dosage_cond_tilde, sparseSigma)
      covM[2:(Mcond+1), 2:(Mcond+1)] = covMsub

      GratioMatrix_cond = getGratioMatrix(dosage_cond, ratioVec)
      G2tilde_P_G2tilde_inv = solve(covMsub %*% GratioMatrix_cond)

}else{# end of if(isCondition)	
  OUT_cond = NULL
  G2tilde_P_G2tilde_inv = NULL
}
