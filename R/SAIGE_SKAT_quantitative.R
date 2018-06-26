options(stringsAsFactors=F)
#' Run single variant score tests with SPA based on the logistic mixed model.
#'
#' @param dosageFile character. Path to dosage file. Each line contains dosages for a marker to be tested
#' @param dosageFileNrowSkip integer(>=0). Number of lines to be skiped in the dosage file.
#' @param dosageFileNcolSkip integer(>=0). Number of columns to be skiped in the dosage file
#' @param dosageFilecolnamesSkip vector of characters. The column names of the skipped columns. Default: c("SNPID", "CHR", "POS", "Allele0", "Allele1")
#' @param dosageFileChrCol string. The column name for the chromosome column. Must be in the dosageFilecolnamesSkip. Required If LOCO = TRUE and chrom ="".  
#' @param bgenFile character. Path to bgen file. Currently version 1.2 with 8 bit compression is supported
#' @param bgenFileIndex character. Path to the .bgi file (index of the bgen file)
#' @param vcfFile character. Path to vcf file
#' @param vcfFileIndex character. Path to index for vcf file by tabix, ".tbi" by "tabix -p vcf file.vcf.gz"
#' @param vcfField character. genotype field in vcf file to use. "DS" for dosages or "GT" for genotypes. By default, "DS".
#' @param savFile character. Path to sav file
#' @param savFileIndex character. Path to index for sav file .s1r
#' @param idstoExcludeFile character. Path to the file containing variant ids to be excluded from the bgen or vcf file
#' @param idstoIncludeFile character. Path to the file containing variant ids to be included from the bgen or vcf file
#' @param rangestoExcludeFile character. Path to the file containing genome regions to be excluded from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header 
#' @param rangestoIncludeFile character. Path to the file containing genome regions to be included from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header 
#' @param chrom character. string for the chromosome to include from vcf file. Required for vcf file. If LOCO is specified, providing chrom will save computation cost
#' @param start numeric. start genome position to include from vcf file. 
#' @param end numeric. end genome position to include from vcf file. 
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 1. The higher threshold between minMAC and minMAF will be used
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used
#' @param maxMAF numeric. Maximum minor allele frequency of markers to test. By default 0.5. 
#' @param minInfo numeric. Minimum imputation info of markers to test (in bgen file)
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the dosage, vcf, sav, or bgen file with NO header
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param Cutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param IsSparse logical. Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for dichotomous traits. By default, TRUE 
#' @param numLinesOutput numeric. Output results for how many marker each time.    
#' @param SAIGEOutputFile character. Path to the output file containing the SPAGMMAT test results
#' @param IsOutputAFinCaseCtrl logical. Whether to output allele frequency in cases and controls. By default, FALSE
#' @param groupFile character. Path to the group file containing one column "GeneID", and one column for ID of the tested genetic markers delimited by comma. This file is space-delimited can only work with the bgen,sav,and vcf format 
#' @return SAIGEOutputFile
#' @export
SKATtest = function(dosageFile = "",
                 dosageFileNrowSkip = 0, 
                 dosageFileNcolSkip = 0,
                 dosageFilecolnamesSkip = c("SNPID", "CHR", "POS", "Allele0", "Allele1"),
		 dosageFileChrCol = "CHR",   ##for LOCO
		 bgenFile = "",
		 bgenFileIndex = "", 
		 vcfFile = "",
                 vcfFileIndex = "",
		 vcfField = "DS",
		 savFile = "",
		 savFileIndex = "",
		 sampleFile = "", 
		 idstoExcludeFile = "",
		 idstoIncludeFile = "",
		 rangestoExcludeFile = "",
		 rangestoIncludeFile = "",
		 chrom = "",
		 start = 1,
		 end = 250000000,	
		 minMAC = 1, 
                 minMAF = 0,
                 maxMAF = 0.5,
        	 minInfo = 0,
                 GMMATmodelFile = "", 
                 varianceRatioFile = "", 
                 Cutoff=2, 
                 SAIGEOutputFile = "",
		 numLinesOutput = 10000, 
		 IsSparse=TRUE,
		 IsSparseSigma = TRUE,
		 sparseSigmaFile = "", 
		 IsOutputAFinCaseCtrl=FALSE,
		 groupFile=""
){


  #check and read files
  #output file
  if(!file.exists(SAIGEOutputFile)){
    file.create(SAIGEOutputFile, showWarnings = TRUE)
  }

  if(IsSparseSigma){
    if(!file.exists(sparseSigmaFile)){
      stop("ERROR! sparseSigmaFile ", sparseSigmaFile, " does not exsit\n")	
    }else{
      sparseSigma = Matrix:::readMM(sparseSigmaFile)
    }
  }else{

    sparseSigma = NULL
  }
  #file for the glmm null model
  if(!file.exists(GMMATmodelFile)){
    stop("ERROR! GMMATmodelFile ", GMMATmodelFile, " does not exsit\n")
  }else{
    load(GMMATmodelFile)
    obj.glmm.null = modglmm
    sampleInModel = NULL
    sampleInModel$IID = obj.glmm.null$sampleID
    sampleInModel = data.frame(sampleInModel)
    sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
    cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
    #print(sampleInModel$IID[1:10])
    obj.glm.null = obj.glmm.null$obj.glm.null
    obj.noK = obj.glmm.null$obj.noK   
    traitType = obj.glmm.null$traitType 
  }

  if(!file.exists(varianceRatioFile)){
    stop("ERROR! varianceRatioFile ", varianceRatioFile, " does not exsit\n")
  }else{
    varRatio = data.frame(data.table:::fread(varianceRatioFile, header=F, stringsAsFactors=FALSE))
    if(nrow(varRatio) == 1){
      ratioVec = rep(varRatio[1,1],6)
    }else{
      ratioVec = varRatio[,1]
    }
    #cat("variance Ratio is ", varRatio, "\n")
    cat("variance Ratio is ", ratioVec, "\n")
  }


  #sample file
  if(!file.exists(sampleFile)){
    stop("ERROR! sampleFile ", sampleFile, " does not exsit\n")
  }else{
    sampleListinDosage = data.frame(data.table:::fread(sampleFile, header=F, stringsAsFactors=FALSE))
    sampleListinDosage$IndexDose = seq(1,nrow(sampleListinDosage), by=1)
    cat(nrow(sampleListinDosage), " sample IDs are found in sample file\n")
    colnames(sampleListinDosage)[1] = "IIDDose"

    dataMerge = merge(sampleInModel, sampleListinDosage, by.x="IID", by.y = "IIDDose")
    dataMerge_sort = dataMerge[with(dataMerge, order(IndexInModel)), ]
    if(nrow(dataMerge_sort) < nrow(sampleInModel)){
      stop("ERROR!", nrow(sampleInModel) - nrow(dataMerge_sort), " samples used in glmm model fit do not have dosages\n")
    }else{
      #0909 modified by WZ
      dataMerge_v2 = merge(dataMerge_sort, sampleListinDosage, by.x="IID", by.y = "IIDDose", all.y = TRUE)
      print(dim(dataMerge_v2))
      print(colnames(dataMerge_v2))
      dataMerge_v2_sort = dataMerge_v2[with(dataMerge_v2, order(IndexDose.y)), ]
      sampleIndex = dataMerge_v2_sort$IndexInModel
      N = sum(!is.na(sampleIndex))
      cat(N, " samples were used in fitting the NULL glmm model and are found in sample file\n")
      sampleIndex[is.na(sampleIndex)] = -10  ##with a negative number
    }
  }


  ##Needs to check the number of columns and the number of samples in sample file
  if(dosageFile != ""){
    if(!file.exists(dosageFile)){
      stop("ERROR! dosageFile ", dosageFile, " does not exsit\n")
    }else{
      if(dosageFileNrowSkip < 0 | dosageFileNcolSkip < 0){
        stop("ERROR! dosageFileNrowSkip or dosageFileNcolSkip can't be less than zero\n")
      }
    }
    dosageFileType = "plain"

  }else if(bgenFile != ""){ 
    if(!file.exists(bgenFile)){
      stop("ERROR! bgenFile ", bgenFile, " does not exsit\n")
    }
    if(!file.exists(bgenFileIndex)){
      stop("ERROR! bgenFileIndex ", bgenFileIndex, " does not exsit\n")
    }
    dosageFileType = "bgen"

  }else if(vcfFile != ""){
    if(!file.exists(vcfFile)){
      stop("ERROR! vcfFile ", vcfFile, " does not exsit\n")
    }
    if(!file.exists(vcfFileIndex)){
      stop("ERROR! vcfFileIndex ", vcfFileIndex, " does not exsit\n")
    }
    dosageFileType = "vcf"

    ###chrom needs to be specified 
    if(chrom == ""){stop("ERROR! chrom needs to be specified for the vcf file\n")}

  }else if(savFile != ""){
    if(!file.exists(savFile)){
      stop("ERROR! savFile ", savFile, " does not exist\n")
    }else{
      vcfFile = savFile	
    }

    if(!file.exists(savFileIndex)){
      stop("ERROR! savFileIndex ", savFileIndex, " does not exist\n")
    }else{
      vcfFileIndex = savFileIndex
    }	
    dosageFileType = "vcf"
  }

  if(dosageFileType != "plain"){
    if(!file.exists(groupFile)){
      stop("ERROR! groupFile ", groupFile, " does not exist\n")
    }
  }

  #determine minimum MAF for markers to be tested
  if(minMAC < 1){minMAC = 1} ##01-19-2018
  cat("minMAC: ",minMAC,"\n")
  cat("minMAF: ",minMAF,"\n")
  cat("maxMAF: ",maxMAF,"\n")

  minMAFBasedOnMAC = minMAC/(2*N) 
  testMinMAF = max(minMAFBasedOnMAC, minMAF) 
  cat("Minimum MAF of markers to be testd is ", testMinMAF, "\n")
  
  if(dosageFileType == "vcf"){ setMAFcutoffs(testMinMAF, maxMAF) }


  ##############START TEST########################
  startTime = as.numeric(Sys.time())  # start time of the SPAGMMAT tests
  cat("Analysis started at ", startTime, "Seconds\n")

  if(file.exists(SAIGEOutputFile)){file.remove(SAIGEOutputFile)}

  sampleIndex = sampleIndex - 1

isVariant = TRUE
if(dosageFileType == "plain"){

}else if (dosageFileType == "bgen"){
  SetSampleIdx(sampleIndex, N)

}else if(dosageFileType == "vcf"){
  isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
  SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N) 
}
#cat("sampleIndex: ", sampleIndex, "\n")

  if(traitType == "quantitative"){
    OUT = NULL
    cat("It is a quantitative trait\n")
    mth = 0
    resultHeader = c("Gene", "Pvalue", "N_MAC1","N_MAC2","N_MAC3","N_MAC4","N_MAC5","N_MACgt5","markerIDs","markerAFs")
    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))

    gf = file(groupFile, "r")
    while ( TRUE ) {
      marker_group_line = readLines(gf, n = 1)
      print(marker_group_line)
      geneID = strsplit(marker_group_line, split="\t")[[1]][1]
      if(length(line) == 0 ){
        break
      }else{
        if(dosageFileType == "vcf"){
          Gx = getGenoOfGene_vcf(marker_group_line)
        }else if(dosageFileType == "bgen"){
          Gx = getGenoOfGene_bgen(bgenFile,bgenFileIndex,marker_group_line, testMinMAF, maxMAF)          
        }

        G0 = Gx$dosages
        cntMarker = Gx$cnt
#	cat("markerIDs: ", Gx$markerIDs, "\n")
#	cat("G0: ", G0, "\n")
        if(cntMarker > 0){
          Gmat = matrix(G0, byrow=F, ncol = cntMarker)
#	 cat("Gmat[,1]: ", Gmat[,1], "\n")	
#	 cat("Gmat[,2]: ", Gmat[,2], "\n")	
	#  saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec)

	  saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec, isSparseSigma = TRUE, sparseSigma = sparseSigma)
	
          OUT = rbind(OUT, c(geneID, saigeskatTest$p.value, saigeskatTest$markerNumbyMAC, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))
          mth = mth + 1
          if(mth %% numLinesOutput == 0){
	    ptm <- proc.time()
            print(ptm)
            print(mth)
            OUT = as.data.frame(OUT)
            write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
            OUT = NULL
          }
        }    
      }#end of else for if(length(line) == 0 )
    } # end of while ( TRUE ) {
    if(!is.null(OUT)){
      OUT = as.data.frame(OUT)
      write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
      OUT = NULL
    }
  }else{
    stop("ERROR! The type of the trait has to be quantitative\n")
  }

  #close the dosage file after tests
  if(dosageFileType == "plain"){
    closetestGenoFile_plainDosage()  
  }else if (dosageFileType == "bgen"){
    closetestGenoFile_bgenDosage()
  }else if(dosageFileType == "vcf"){
    closetestGenoFile_vcfDosage()
  }


  endTime = as.numeric(Sys.time()) #end time of the SPAGMMAT tests
  cat("Analysis ended at ", endTime, "Seconds\n")
  tookTime = endTime - startTime
  cat("Analysis took ", tookTime, "Seconds\n")
  
}


#No Sparsity
scoreTest_SAIGE_quantitativeTrait_old=function(G0, obj.noK, AC, y, mu, varRatio, tauVec){
  XVG0 = eigenMapMatMult(obj.noK$XV, G0)
  G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G1 is X adjusted 
  g = G/sqrt(AC)
  var2 = innerProduct(g, g)
  q = innerProduct(g, y)
  m1 = innerProduct(mu, g)
  var1 = var2 * varRatio
  Tstat = (q-m1)/tauVec[1]
  p.value = pchisq(Tstat^2/var1, lower.tail = FALSE, df=1)
  BETA = (Tstat/var1)/sqrt(AC)
  SE = abs(BETA/qnorm(p.value/2))
  out1 = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
  return(out1)
}

#Use Sparsity trick for rare variants
scoreTest_SAIGE_quantitativeTrait=function(G0, obj.noK, AC, AF, y, mu, varRatio, tauVec){
  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }
  maf = min(AF, 1-AF)
  if(maf < 0.05){
    idx_no0<-which(G0>0)
    #cat("length(idx_no0): ", length(idx_no0), "\n")
    #cat("maf: ", maf, "\n")
    g1<-G0[idx_no0]/sqrt(AC2)
    A1<-obj.noK$XVX_inv_XV[idx_no0,]
    X1<-obj.noK$X1[idx_no0,]
    mu1<-mu[idx_no0]
    y1<-obj.noK$y[idx_no0]


## V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv
    if(length(idx_no0) > 1){
      Z = t(A1) %*% g1
      B<-X1 %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      S_a2 = obj.noK$S_a - colSums(X1 * (y1 - mu1))
      S2 = -S_a2 %*% Z
    }else{
      Z = A1 * g1    
      B<-X1 %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      S_a2 = obj.noK$S_a - X1 * (y1 - mu1)
      S2 = -S_a2 %*% Z
    }
    S<- S1+S2
    Tstat = S/tauVec[1]
  }else{
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G1 is X adjusted
    g = G/sqrt(AC2)
    q = innerProduct(g, y)
    m1 = innerProduct(mu, g)
    var2 = innerProduct(g, g)
    var1 = var2 * varRatio
    Tstat = (q-m1)/tauVec[1]
  }
  p.value = pchisq(Tstat^2/var1, lower.tail = FALSE, df=1)
  BETA = (Tstat/var1)/sqrt(AC2)
  SE = abs(BETA/qnorm(p.value/2))
  out1 = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
  return(out1)
}


Score_Test_Sparse<-function(obj.null, G, mu, mu2, varRatio ){
  # mu=mu.a; mu2= mu2.a; G=G0; obj.null=obj.noK
  idx_no0<-which(G>0)
  #print(G[1:10])
  #print(idx_no0)
  g1<-G[idx_no0]
  A1<-obj.null$XVX_inv_XV[idx_no0,]
  X1<-obj.null$X1[idx_no0,]
  mu21<-mu2[idx_no0]
  mu1<-mu[idx_no0]
  y1<-obj.null$y[idx_no0]

  if(length(idx_no0) > 1){
    #cat("idx_no0 ", idx_no0, "\n")
    #cat("dim(X1) ", X1, "\n")
    Z = t(A1) %*% g1
    B<-X1 %*% Z
    #cat("dim(Z) ", Z, "\n")
    g_tilde1 = g1 - B
    var2 = t(Z) %*% obj.null$XVX %*% Z - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
    var1 = var2 * varRatio
    S1 = crossprod(y1-mu1, g_tilde1)
    S_a2 = obj.null$S_a - colSums(X1 * (y1 - mu1))
    S2 = -S_a2 %*% Z
  }else{
    #cat("idx_no0 ", idx_no0, "\n")
    #cat("dim(X1) ", X1, "\n")
    Z = A1 * g1
    #cat("dim(Z) ", Z, "\n")
    #cat("dim(Z) here ", Z, "\n")
    B<-X1 %*% Z
    g_tilde1 = g1 - B
    var2 = t(Z) %*% obj.null$XVX %*% Z - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
    var1 = var2 * varRatio
    S1 = crossprod(y1-mu1, g_tilde1)
    S_a2 = obj.null$S_a - X1 * (y1 - mu1)
    S2 = -S_a2 %*% Z
  }

  S<- S1+S2
	
  pval.noadj<-pchisq((S)^2/(var1), lower.tail = FALSE, df=1)
  ##add on 10-25-2017
  BETA = S/var1
  SE = abs(BETA/qnorm(pval.noadj/2))
  Tstat = S
  #return(c(BETA, SE, Tstat, pval.noadj, pval.noadj, 1, var1, var2))
  return(list(BETA=BETA, SE=SE, Tstat=Tstat, pval.noadj=pval.noadj, pval.noadj=pval.noadj, is.converge=TRUE, var1=var1, var2=var2))	
}




Score_Test<-function(obj.null, G, mu, mu2, varRatio){
  g<-G  -  obj.null$XXVX_inv %*%  (obj.null$XV %*% G)
  q<-crossprod(g, obj.null$y) 
  m1<-crossprod(mu, g)
  var2<-crossprod(mu2, g^2)
  var1 = var2 * varRatio
  S = q-m1

  pval.noadj<-pchisq((S)^2/var1, lower.tail = FALSE, df=1)

  ##add on 10-25-2017
  BETA = S/var1
  SE = abs(BETA/qnorm(pval.noadj/2))
  #Tstat = S^2
  Tstat = S

  #return(c(BETA, SE, Tstat, pval.noadj, pval.noadj, NA, var1, var2))
  #return(c(pval.noadj, pval.noadj, TRUE, var1, var2))
  return(list(BETA=BETA, SE=SE, Tstat=Tstat, pval.noadj=pval.noadj, pval.noadj=pval.noadj, is.converge=TRUE, var1=var1, var2=var2))
}


####add log(OR), SE, and T estimation on 10-25-2017#######
scoreTest_SPAGMMAT_binaryTrait=function(g, AC, NAset, y, mu, varRatio, Cutoff){
        #g = G/sqrt(AC)
  q = innerProduct(g, y)
  m1 = innerProduct(g, mu)
  var2 = innerProduct(mu*(1-mu), g*g)
  var1 = var2 * varRatio
  Tstat = q-m1
        #cat("Tstat: ", Tstat, "\n")
  qtilde = Tstat/sqrt(var1) * sqrt(var2) + m1
        #cat("var1: ", var1, "\n")


  if(length(NAset)/length(g) < 0.5){
    out1 = SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
  }else{
    out1 = SPAtest:::Saddle_Prob_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8)
  }

  out1 = c(out1, var1 = var1)
  out1 = c(out1, var2 = var2)
  #logOR = Tstat0/var1
  #logOR = Tstat0/(sqrt(var1)*sqrt(var2))
  logOR = (Tstat/var1)/sqrt(AC)
  SE = abs(logOR/qnorm(out1$p.value/2))
  out1 = c(out1, BETA = logOR, SE = SE, Tstat = Tstat)
  return(out1)
}


###add on 10-25-2017###for score test for binary traits for IsSparse 
####add log(OR), SE, and T estimation on 10-25-2017#######
scoreTest_SAIGE_binaryTrait=function(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y,varRatio, Cutoff, rowHeader){
  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }

  ##########################
  ## Added by SLEE 09/06/2017
  Run1=TRUE
  if(IsSparse==TRUE){
    if(MAF < 0.05){
       out.score<-Score_Test_Sparse(obj.noK, G0,mu.a, mu2.a, varRatio );
     }else{
       out.score<-Score_Test(obj.noK, G0,mu.a, mu2.a, varRatio );
     }
     if(out.score["pval.noadj"] > 0.05){
       if(AF > 0.5){
         out.score$BETA = (-1)*out.score$BETA
         out.score$Tstat = (-1)*out.score$Tstat
         #out.score["BETA"][1] = (-1)*out.score["BETA"][1]
         #out.score["Tstat"][1] = (-1)*out.score["Tstat"][1]
       }

       #OUT = rbind(OUT, c(rowHeader, N, unlist(out.score)))
       outVec = c(rowHeader, N, unlist(out.score))
       #NSparse=NSparse+1
       Run1=FALSE
       	
     }
  }

  if(Run1){
    G0 = matrix(G0, ncol = 1)
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
    g = G/sqrt(AC2)
    NAset = which(G0==0)
    out1 = scoreTest_SPAGMMAT_binaryTrait(g, AC2, NAset, y, mu.a, varRatio, Cutoff = Cutoff)
    if(AF > 0.5){
      out1$BETA = (-1)*out1$BETA
      out1$Tstat = (-1)*out1$Tstat
    }
    out1 = unlist(out1)

    #OUT = rbind(OUT, c(rowHeader, N, out1["BETA"], out1["SE"], out1["Tstat"], out1["p.value"], out1["p.value.NA"], out1["Is.converge"], out1["var1"], out1["var2"]))
    outVec = c(rowHeader, N, out1["BETA"], out1["SE"], out1["Tstat"], out1["p.value"], out1["p.value.NA"], out1["Is.converge"], out1["var1"], out1["var2"])
   }
  return(outVec)
}


#       obj for SAIGE_SKAT
#               ratioVec: vector for variance ratio parameter
#               P0: P0 from intermediate sparse Kinship
#               res: residual from the full Kinship
#

SAIGE_SKAT_withRatioVec  = function( Z, obj, ratioVec, kernel= "linear.weighted", method="davies", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, isSparseSigma, sparseSigma, ratioVec1 = NULL){
        m = ncol(Z)
        n = nrow(Z)

        id_include<-1:n
        #if(class(obj) != "SKAT_NULL_Model_SAIGE"){
        #       stop("obj is not the returned object from SAIGE_SKAT_Read_Null_File")
        #}

        # Added by SLEE 4/24/2017
        out.method<-SKAT:::SKAT_Check_Method(method,r.corr, n=n, m=m)
        method=out.method$method
        r.corr=out.method$r.corr
        IsMeta=out.method$IsMeta
        SKAT:::SKAT_Check_RCorr(kernel, r.corr)

        out.z<-SKAT:::SKAT_MAIN_Check_Z(Z, n, id_include, SetID, weights, weights.beta, impute.method, is_check_genotype
        , is_dosage, missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)
        if(out.z$return ==1){
                out.z$param$n.marker<-m
                return(out.z)
        }
        Z = out.z$Z.test
        weights = out.z$weights
        #res = as.numeric(obj$residuls)/(as.numeric(obj$theta[1]))

        ##process variance ratio
        cat("dim(Z) is ", dim(Z), "\n")
        MACvec = colSums(Z)
        MACvec_indVec = MACvec
        MACvec_indVec[which(MACvec <= 1.5)] = 1
        MACvec_indVec[which(MACvec <= 2.5 & MACvec > 1.5)] = 2
        MACvec_indVec[which(MACvec <= 3.5 & MACvec > 2.5)] = 3
        MACvec_indVec[which(MACvec <= 4.5 & MACvec > 3.5)] = 4
        MACvec_indVec[which(MACvec <= 5.5 & MACvec > 4.5)] = 5
        MACvec_indVec[which(MACvec_indVec > 5.5)] = 6

        cat("MACvec_indVec: ", MACvec_indVec, "\n")
        indMatrix = contr.sum(6, contrasts = FALSE)
        GindMatrix = NULL
        for(i in MACvec_indVec){
          GindMatrix = rbind(GindMatrix, indMatrix[i,])
        }

        #mx1 = mx6 %*% 6x1
        GratioVec = GindMatrix %*% matrix(ratioVec, ncol=1)
        #mxm
        GratioMatrix = sqrt(GratioVec) %*% t(sqrt(GratioVec))

        if(!is.null(ratioVec1)){
                GratioVec1 = GindMatrix %*% matrix(ratioVec1, ncol=1)
                #mxm
                GratioMatrix1 = sqrt(GratioVec1) %*% t(sqrt(GratioVec1))
        }

        cat("dim(Z) is ", dim(Z), "\n")
        # If Z is sparse, change it to the sparse matrix
        if(mean(Z) < 0.1){
                Z = as(Z, "sparseMatrix")
        }

        if (kernel == "linear.weighted") {
        Z = t(t(Z) * (weights))
        }

        # Get Score
 #       cat("dim(Z) is ", dim(Z), "\n")
 #       print(dim(t(Z)))
 #       print(dim(obj$residuls))

        Score = as.vector(t(Z) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])


        if(isSparseSigma){
                pcginvSigma<-NULL
                for(i in 1:ncol(Z)){
                        c3<-pcg(sparseSigma, Z[,i])
                        pcginvSigma<-cbind(pcginvSigma, c3)
                }
                Phi = as.matrix(t(Z) %*% pcginvSigma)
        }else{
                # Phi
        if(!is.null(obj$P)){
                Phi = t(Z) %*% (obj$P %*% Z)
        } else {
                XVZ = obj$obj.noK$XV %*% Z
                Z1 = Z  -  (obj$obj.noK$XXVX_inv %*% XVZ) # G1 is X adjusted
#                Score = as.vector(t(Z) %*% obj$residuals)/as.numeric(obj$theta[1])
                Phi = t(Z) %*% Z1
        }

        }
        #Phi = Phi * obj$ratio
#        print(Phi)
        PhiMulti = Phi * GratioMatrix

        cat("Phi[1:10]: ", Phi[1:10], "\n")
        re = SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=PhiMulti, r.corr=r.corr, method=method, Score.Resampling=NULL)
        re$IsMeta=TRUE

        ##summaize the number of markers falling in each MAC category
        markerNumbyMAC = c(sum(MACvec_indVec == 1), sum(MACvec_indVec == 2), sum(MACvec_indVec == 3), sum(MACvec_indVec == 4), sum(MACvec_indVec == 5), sum(MACvec_indVec == 6))
        re$markerNumbyMAC = markerNumbyMAC

        #PhiSingle = Phi * GratioMatrixSingle
        #reSingle = SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=PhiSingle, r.corr=r.corr, method=method, Score.Resampling=NULL)
        #re$p.value.Single = reSingle$p.value
        if(!is.null(ratioVec1)){
                Phi1 = Phi * GratioMatrix1
                re1 = SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi1, r.corr=r.corr, method=method, Score.Resampling=NULL)
                re$p.value.1 = re1$p.value
        }else{
                re$p.value.1 = NA
        }

        print("re")
        print(re)
        return(re)
}

