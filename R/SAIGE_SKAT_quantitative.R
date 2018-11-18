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
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 0. The higher threshold between minMAC and minMAF will be used
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
#' @param condition. For conditional analysis. Genetic marker ids (chr:pos_ref/alt) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A. Note that currently conditional analysis is only for vcf/sav input.
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
		 minMAC = 0, 
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
		 groupFile="",
		 condition=""
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
    #stop("ERROR! varianceRatioFile ", varianceRatioFile, " does not exsit\n")
  }else{
    varRatioData = data.frame(data.table:::fread(varianceRatioFile, header=F, stringsAsFactors=FALSE))
    if(nrow(varRatioData) == 1){
      ratioVec = rep(varRatioData[1,1],6)
    }else{
      ratioVec = varRatioData[,1]
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
      sampleIndex = sampleIndex - 1	
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


  if(condition != ""){
    isCondition = TRUE
  }else{
    isCondition = FALSE
  }


  if(isCondition){
    conditionlist = paste(c("condMarkers",unlist(strsplit(condition,","))),collapse="\t")
    cat("conditionlist is ", conditionlist, "\n")
    if(dosageFileType == "vcf"){
      setMAFcutoffs(0, 0.5)
      isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
      SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N)
      Gx_cond = getGenoOfGene_vcf(conditionlist)
    }else if(dosageFileType == "bgen"){
      Gx_cond = getGenoOfGene_bgen(bgenFile,bgenFileIndex, conditionlist)
    }else{
      cat("WARNING: conditional analysis can only work for dosageFileType vcf, sav or bgen\n")
    }
    #print(Gx_cond)
    cat("conditioning on ", unlist(Gx_cond$markerIDs), "\n")
    #G0 = Gx_cond$dosages
     cntMarker = Gx_cond$cnt

     if(cntMarker > 0){
          dosage_cond = matrix(Gx_cond$dosages, byrow=F, ncol = cntMarker)
     }
    print(dim(dosage_cond))


	


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
    resultHeader = c(dosageFilecolnamesSkip,  "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar")
    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
    OUT = NULL
    numPassMarker = 0
    mth = 0
   # sampleIndex = sampleIndex - 1
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
    }else{
      cat("LOCO will be used, but chromosome for the dosage file is not specified. Will check each marker for its chromosome for LOCO!\n")
      indChromCheck = TRUE
    }

  }else{
    stop("ERROR! The type of the trait has to be either binary or quantitative\n")
  }


    OUT_cond = NULL
    for(i in 1:ncol(dosage_cond)){
    G0  = dosage_cond[,i]
    AC = sum(G0)
    N  = length(G0)
    AF = AC/(2*N)
    MAF = AF
    if(AF > 0.5){
      MAF = 1-AF
      MAC = 2*N - AC
    }else{
      MAC = AC	
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
	cat("G0: ",G0,"\n")	
    #    out1 = scoreTest_SAIGE_quantitativeTrait(G0, obj.noK, AC, AF, y, mu, varRatio, tauVec)
	out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.noK, AC, AF, y, mu, varRatio, tauVec, sparseSigma=sparseSigma)
	print("out1")
	print(out1)
        OUT_cond = rbind(OUT_cond, c(as.numeric(out1$BETA), as.numeric(out1$Tstat), as.numeric(out1$var1)))
    }
  OUT_cond = as.matrix(OUT_cond)

 } # end of for(i in 1:ncol(dosage_cond)){


}else{ # if(isCondition)
  dosage_cond = NULL
  OUT_cond = NULL

}

  #determine minimum MAF for markers to be tested
  #if(minMAC < 1){minMAC = 1} ##01-19-2018
  if(!(minMAC  > 0)){minMAC == 0} ##07-15-2018
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

#  sampleIndex = sampleIndex - 1



if(dosageFileType == "plain"){
  isCondition = FALSE
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
#    resultHeader = c("Gene", "Pvalue", "N_MAC1","N_MAC2","N_MAC3","N_MAC4","N_MAC5","N_MACgt5","markerIDs","markerAFs")
    resultHeader = c("Gene","Pvalue","Pvalue_cond","N_MAC1","N_MAC2","N_MAC3","N_MAC4","N_MAC5","N_MACgt5","markerIDs","markerAFs")
    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))

    gf = file(groupFile, "r")
    while ( TRUE ) {
      marker_group_line = readLines(gf, n = 1)
    #  print(marker_group_line)
      if(length(marker_group_line) == 0 ){
        break
      }else{
      	geneID = strsplit(marker_group_line, split="\t")[[1]][1]
        if(dosageFileType == "vcf"){
          Gx = getGenoOfGene_vcf(marker_group_line)
        }else if(dosageFileType == "bgen"){
          Gx = getGenoOfGene_bgen(bgenFile,bgenFileIndex,marker_group_line, testMinMAF, maxMAF)          
        }

        G0 = Gx$dosages
        cntMarker = Gx$cnt
	cat("cntMarker: ", cntMarker, "\n")
#	cat("markerIDs: ", Gx$markerIDs, "\n")
#	cat("G0: ", G0, "\n")
        if(cntMarker > 0){
         Gmat = matrix(G0, byrow=F, ncol = cntMarker)	  
#	 cat("dim(Gmat): ", dim(Gmat), "\n")	
#	 cat("Gmat[,1]: ", Gmat[,1], "\n")	
#	 cat("ratioVec: ", ratioVec, "\n")
#	 cat("Gmat[,2]: ", Gmat[,2], "\n")	
	#  saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec)
#	 cat("colSums(Gmat): ", colSums(Gmat), "\n")
	 #saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec, Z_cond=dosage_cond, Z_cond_es=OUT_cond[,1], max_maf = maxMAF, sparseSigma = sparseSigma, method = "optimal.adj")
	#print(Gmat[1:20,1:4])
	#print(colMeans(Gmat))
	 saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec, G2_cond=dosage_cond, G2_cond_es=OUT_cond[,1], max_maf = maxMAF, sparseSigma = sparseSigma, method = "optimal.adj")
		


	  cat("saigeskatTest$p.value: ", saigeskatTest$p.value, "\n")


#	  if(isCondition){	
#		saigeskatTest_cond = SAIGE_SKAT_withRatioVec_cond(Gmat, obj.glmm.null, ratioVec, Z_cond=dosage_cond, Z_cond_es=OUT_cond[,1], sparseSigma = sparseSigma )
		
#	  }
	 OUT = rbind(OUT, c(geneID, saigeskatTest$p.value, saigeskatTest$p.value.cond, saigeskatTest$markerNumbyMAC, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))

#	 if(!isCondition){	
#          OUT = rbind(OUT, c(geneID, saigeskatTest$p.value, saigeskatTest$markerNumbyMAC, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))
#	}else{
#	   OUT = rbind(OUT, c(geneID, saigeskatTest$p.value, saigeskatTest_cond$p.value, saigeskatTest$markerNumbyMAC, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))
#	}

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


#       obj for SAIGE_SKAT
#               ratioVec: vector for variance ratio parameter
#               P0: P0 from intermediate sparse Kinship
#               res: residual from the full Kinship
#

SAIGE_SKAT_withRatioVec_old  = function( Z, obj, ratioVec, kernel= "linear.weighted", method="davies", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, sparseSigma){
        m = ncol(Z)
        n = nrow(Z)

        id_include<-1:n

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

  #      cat("MACvec_indVec: ", MACvec_indVec, "\n")
  #      cat("dim(Z) is ", dim(Z), "\n")
	#print(Z)
        # If Z is sparse, change it to the sparse matrix
        if(mean(Z) < 0.1){
                Z = as(Z, "sparseMatrix")
        }

#	obj.noK = obj$obj.noK
#	Z_tilde = Z  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% Z)


        if (kernel == "linear.weighted") {
        Z = t(t(Z) * (weights))
#        Z_tilde = t(t(Z_tilde) * (weights))
        }
	#Z = Z_tilde

	#print(Z)
#        Score = as.vector(t(Z) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
#	cat("Score: ", Score, "\n")


        Z_tilde = Z  -  obj$obj.noK$XXVX_inv %*%  (obj$obj.noK$XV %*% Z)
	Score_tilde = as.vector(t(Z_tilde) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
	cat("Score_tilde: ", Score_tilde, "\n")
	Score = Score_tilde
	Z = Z_tilde

if(is.null(obj$P)){
        if(!is.null(sparseSigma)){
                pcginvSigma<-NULL
                for(i in 1:ncol(Z)){
                        c3<-pcg(sparseSigma, Z[,i])
                        pcginvSigma<-cbind(pcginvSigma, c3)
                }
                Phi = as.matrix(t(Z) %*% pcginvSigma)
		cat("Phi: ", Phi, "\n")

        }else{
		XVZ = obj$obj.noK$XV %*% Z
                Z1 = Z  -  (obj$obj.noK$XXVX_inv %*% XVZ) # G1 is X adjusted
#                Score = as.vector(t(Z) %*% obj$residuals)/as.numeric(obj$theta[1])
                Phi = t(Z) %*% Z1
	
	}

	indMatrix = contr.sum(6, contrasts = FALSE)
        GindMatrix = NULL
        for(i in MACvec_indVec){
          GindMatrix = rbind(GindMatrix, indMatrix[i,])
        }

        #mx1 = mx6 %*% 6x1
        GratioVec = GindMatrix %*% matrix(ratioVec, ncol=1)
        #mxm
        GratioMatrix = sqrt(GratioVec) %*% t(sqrt(GratioVec))

        #cat("GratioMatrix: ", GratioMatrix, "\n")

        if(!is.null(ratioVec1)){
                GratioVec1 = GindMatrix %*% matrix(ratioVec1, ncol=1)
                #mxm
                GratioMatrix1 = sqrt(GratioVec1) %*% t(sqrt(GratioVec1))
        }


	Phi = Phi * GratioMatrix
}else{
                # Phi
        Phi = t(Z) %*% (obj$P %*% Z)

}


        cat("Phi: ", Phi, "\n")
        cat("Score[1:10]: ", Score[1:10], "\n")

        re = SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
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

commentOut=function(a=1){
#SAIGE_SKAT_withRatioVec  = function(Z, obj, ratioVec, Z_cond = NULL, Z_cond_es, kernel= "linear.weighted", method="davies", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
SAIGE_SKAT_withRatioVec  = function(Z, obj, ratioVec, Z_cond = NULL, Z_cond_es, kernel= "linear.weighted", method="optimal.adj", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, sparseSigma = NULL){

	obj.noK = obj$obj.noK
        m = ncol(Z)
        n = nrow(Z)

        id_include<-1:n
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
	print(Z[1,])

	m = ncol(Z)

if(m >  0){
	if(!is.null(Z_cond)){
                m_cond = ncol(Z_cond)
                Zall = cbind(Z, Z_cond)
        }else{
                Zall = Z
        }

        GratioMatrixall = getGratioMatrix(Zall, ratioVec)

        cat("dim(GratioMatrixall) is ", dim(GratioMatrixall), "\n")
        print(GratioMatrixall)

        MACvec_indVec = getMACvec_indVec(Z)

       # cat("MACvec_indVec: ", MACvec_indVec, "\n")
       # cat("m is ", m , "\n")
       # cat("m_cond is ", m_cond , "\n")
        ##summaize the number of markers falling in each MAC category
        markerNumbyMAC = c(sum(MACvec_indVec == 1), sum(MACvec_indVec == 2), sum(MACvec_indVec == 3), sum(MACvec_indVec == 4), sum(MACvec_indVec == 5), sum(MACvec_indVec == 6))




#         If Z is sparse, change it to the sparse matrix
        if(mean(Z) < 0.1){
                Z = as(Z, "sparseMatrix")
        }

        if (kernel == "linear.weighted") {
	        Z = t(t(Z) * (weights))
	        #Z_tilde = t(t(Z_tilde) * (weights))
        }

        Z_tilde = Z  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% Z)

#        Score = as.vector(t(Z) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
	if(!is.null(Z_cond)){
		Z_cond_tilde<- Z_cond  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% Z_cond)

		cat("Z_cond_tilde: ", Z_cond_tilde[1:10], "\n")
		cat("Z_cond_es ", Z_cond_es, "\n")
		cat("Z_cond_tilde%*%Z_cond_es: ", (Z_cond_tilde%*%Z_cond_es)[1:10], "\n")

		Score_cond = as.vector(t(Z) %*% matrix(obj$residuals - Z_cond_tilde%*%Z_cond_es, ncol=1)) / as.numeric(obj$theta[1])
#	Score = as.vector(t(Z_tilde) %*% matrix(obj$residuals - Z_cond_tilde%*%Z_cond_es, ncol=1)) / as.numeric(obj$theta[1])
	}
	Score = as.vector(t(Z) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
	

	if(is.null(obj$P)){

        	if(!is.null(sparseSigma)){
			if(!is.null(Z_cond)){
		#		G1_tilde_Ps_G1_tilde = getcovM(Z_tilde, Z_tilde, sparseSigma)
		#		G2_tilde_Ps_G2_tilde = getcovM(Z_cond_tilde, Z_cond_tilde, sparseSigma)
		#		G1_tilde_Ps_G2_tilde = getcovM(Z_tilde, Z_cond_tilde, sparseSigma)
		#		G2_tilde_Ps_G1_tilde = getcovM(Z_cond_tilde, Z_tilde, sparseSigma)

				G1_tilde_Ps_G1_tilde = t(Z_tilde)%*% solve(sparseSigma) %*% Z_tilde
				G2_tilde_Ps_G2_tilde = t(Z_cond_tilde)%*% solve(sparseSigma) %*% Z_cond_tilde
				G1_tilde_Ps_G2_tilde = t(Z_tilde)%*% solve(sparseSigma) %*% Z_cond_tilde
				G2_tilde_Ps_G1_tilde = t(Z_cond_tilde)%*% solve(sparseSigma) %*% Z_tilde

		#		cat("GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]: ", GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))], "\n")		
				Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - (G1_tilde_Ps_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]))) %*% (G2_tilde_Ps_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m]))

				Phi_cond = as.matrix(Phi_cond)
			}

			Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])
			
		}else{
			if(!is.null(Z_cond)){
				G1_tilde_G1_tilde = t(Z_tilde) %*% Z_tilde
				G2_tilde_G2_tilde = t(Z_cond_tilde) %*% Z_cond_tilde	
				G1_tilde_G2_tilde = t(Z_tilde) %*% Z_cond_tilde
				G2_tilde_G1_tilde = t(Z_cond_tilde) %*% Z_tilde

				Phi_cond = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m]) - (G1_tilde_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]))) %*% (G2_tilde_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])) 

			}

			Phi = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m])
			

		}

        }else{
		if(!is.null(Z_cond)){
                	Phi_cond = t(Z_tilde) %*% (obj$P %*% Z_tilde) - (t(Z_tilde) %*% (obj$P %*% Z_cond_tilde)) %*% solve(t(Z_cond_tilde) %*% (obj$P %*% Z_cond_tilde)) %*% (t(Z_cond_tilde) %*% (obj$P %*% Z_tilde))
		}
                Phi = t(Z_tilde) %*% (obj$P %*% Z_tilde)
		
        }


        cat("diag(Phi_cond): ", diag(Phi_cond), "\n")
        cat("Score_cond: ", Score_cond, "\n")
#	dP = diag(Phi)
#	if(sum(dP < 10^-5) > 0){
#		cat("dP[which(dP < 10^-5)]: ",  dP[which(dP < 10^-5)] , "\n")
#		re = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)
#	}else{
#	}

	if(!is.null(Z_cond)){
		if(sum(diag(Phi_cond)) < 10^-5){
			re_cond = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)
		}else{
        		re_cond = SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL)
		}
	}

	re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)


	if(!is.null(Z_cond)){
		re$p.value.cond = re_cond$p.value
	}else{
		re$p.value.cond = NA
	}

        print(re)

  }else{ #no marker is left for test
    re = list(p.value = NA, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA)
    markerNumbyMAC = c(0,0,0,0,0,0)

  }

  re$IsMeta=TRUE
  re$markerNumbyMAC = markerNumbyMAC
  return(re)

}




getGratioMatrix = function(G, ratioVec){

	MACvec = colSums(G)
        MACvec_indVec = MACvec
        MACvec_indVec[which(MACvec <= 1.5)] = 1
        MACvec_indVec[which(MACvec <= 2.5 & MACvec > 1.5)] = 2
        MACvec_indVec[which(MACvec <= 3.5 & MACvec > 2.5)] = 3
        MACvec_indVec[which(MACvec <= 4.5 & MACvec > 3.5)] = 4
        MACvec_indVec[which(MACvec <= 5.5 & MACvec > 4.5)] = 5
        MACvec_indVec[which(MACvec_indVec > 5.5)] = 6

       # cat("MACvec_indVec: ", MACvec_indVec, "\n")

        indMatrix = contr.sum(6, contrasts = FALSE)

        GindMatrix = NULL
        for(i in MACvec_indVec){
          GindMatrix = rbind(GindMatrix, indMatrix[i,])
        }

        #mx1 = mx6 %*% 6x1
        GratioVec = GindMatrix %*% matrix(ratioVec, ncol=1)
        #mxm
        GratioMatrix = sqrt(GratioVec) %*% t(sqrt(GratioVec))  

	#re = list(GratioMatrix = GratioMatrix
	return(GratioMatrix)
}


getcovM = function(G1, G2, sparseSigma){
   pcginvSigma = NULL
   for(i in 1:ncol(G2)){
     c3<-pcg(sparseSigma, G2[,i])
     pcginvSigma<-cbind(pcginvSigma, c3)
   }
   covM = as.matrix(t(G1) %*% pcginvSigma)
   return(covM)
}

getvarRatio = function(MAC, ratioVec){
	if(MAC == 1){
		i = 1
	}else if(MAC == 2){
		i = 2
	}else if(MAC == 3){
		i = 3
	}else if(MAC == 4){
		i = 4
	}else if(MAC == 5){
		i = 5
	}else if(MAC > 5){
		i = 6
	}
	varRatio = ratioVec[i]
	return(varRatio)
}

getMACvec_indVec = function(Z){

	MACvec = colSums(Z)
        MACvec_indVec = MACvec
        MACvec_indVec[which(MACvec <= 1.5)] = 1
        MACvec_indVec[which(MACvec <= 2.5 & MACvec > 1.5)] = 2
        MACvec_indVec[which(MACvec <= 3.5 & MACvec > 2.5)] = 3
        MACvec_indVec[which(MACvec <= 4.5 & MACvec > 3.5)] = 4
        MACvec_indVec[which(MACvec <= 5.5 & MACvec > 4.5)] = 5
        MACvec_indVec[which(MACvec_indVec > 5.5)] = 6

	return(MACvec_indVec)	
}

}




#obj is the rda. file output from SAIGE step 1
#G1 is genotypes for testing gene, which contains m markers
#G2_cond is G2 in the word document, genotypes for m_cond conditioning marker(s)
#G2_cond_es is beta_2_hat (effect size for the conditioning marker(s))


SAIGE_SKAT_withRatioVec_oldAugust13  = function(G1, obj, ratioVec, G2_cond = NULL, G2_cond_es, kernel= "linear.weighted", method="optimal.adj", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, sparseSigma = NULL){

	#check the input genotype G1 
	obj.noK = obj$obj.noK
        m = ncol(G1)
        n = nrow(G1)
	cat("m =", m, "\n") 
	#cat(colSums(G1),"\n")
	



        id_include<-1:n
        # Added by SLEE 4/24/2017
        out.method<-SKAT:::SKAT_Check_Method(method,r.corr, n=n, m=m)
        method=out.method$method
        r.corr=out.method$r.corr
        IsMeta=out.method$IsMeta
        SKAT:::SKAT_Check_RCorr(kernel, r.corr)

        out.z<-SKAT:::SKAT_MAIN_Check_Z(G1, n, id_include, SetID, weights, weights.beta, impute.method, is_check_genotype
        , is_dosage, missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)


	if(out.z$return ==1){
                out.z$param$n.marker<-m
                #return(out.z)
		m = 0
        }else{

        	G1 = out.z$Z.test
        	weights = out.z$weights
		m = ncol(G1)
	}
	cat("m", m, "\n")
	#if more than 1 marker is left, continue the test
	if(m  >  0){
		#cbind G1 and G2_cond to estimate the variance ratio matrix (m+m_cond) x (m+m_cond) 
		if(!is.null(G2_cond)){
                	m_cond = ncol(G2_cond)
                	Zall = cbind(G1, G2_cond)
        	}else{
                	Zall = G1
        	}

        	GratioMatrixall = getGratioMatrix(Zall, ratioVec)
        	MACvec_indVec = getMACvec_indVec(G1)
        	##summaize the number of markers falling in each MAC category
        	markerNumbyMAC = c(sum(MACvec_indVec == 1), sum(MACvec_indVec == 2), sum(MACvec_indVec == 3), sum(MACvec_indVec == 4), sum(MACvec_indVec == 5), sum(MACvec_indVec == 6))


#	         If G1 is sparse, change it to the sparse matrix
        	if(mean(G1) < 0.1){
                  G1 = as(G1, "sparseMatrix")
        	}

        	if (kernel == "linear.weighted") {
	        	G1 = t(t(G1) * (weights))
	        	#Z_tilde = t(t(Z_tilde) * (weights))
        	}

        	G1_tilde = G1  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G1)
		cat("dim(G1)", dim(G1), "\n")
		Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
		cat("dim(Score)", length(Score), "\n")

		#compute Score test statistics after conditionining		
		if(!is.null(G2_cond)){
			G2_cond_tilde<- G2_cond  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G2_cond)
			T2 = as.vector(t(G2_cond) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
			#Score_cond = as.vector(t(G1) %*% matrix(obj$residuals - G2_cond_tilde%*%G2_cond_es, ncol=1)) / as.numeric(obj$theta[1])
		}

	

		#if no P is provides, use sparseSigma or identity Sigma
		if(is.null(obj$P)){

        		if(!is.null(sparseSigma)){
#				cat("first\n")
				#G1_tilde_Ps_G1_tilde = t(G1_tilde)%*% solve(sparseSigma) %*% G1_tilde
				G1_tilde_Ps_G1_tilde = getcovM(G1_tilde, G1_tilde, sparseSigma)
				if(!is.null(G2_cond)){
#				cat("second\n")
		#		G1_tilde_Ps_G1_tilde = getcovM(Z_tilde, Z_tilde, sparseSigma)
				G2_tilde_Ps_G2_tilde = getcovM(G2_cond_tilde, G2_cond_tilde, sparseSigma)
#				cat("second2\n")
				G1_tilde_Ps_G2_tilde = getcovM(G1_tilde, G2_cond_tilde, sparseSigma)
#				cat("second3\n")
				G2_tilde_Ps_G1_tilde = getcovM(G2_cond_tilde, G1_tilde, sparseSigma)
#				cat("second4\n")

		#		G2_tilde_Ps_G2_tilde = t(G2_cond_tilde)%*% solve(sparseSigma) %*% G2_cond_tilde
		#		G1_tilde_Ps_G2_tilde = t(G1_tilde)%*% solve(sparseSigma) %*% G2_cond_tilde
		#		G2_tilde_Ps_G1_tilde = t(G2_cond_tilde)%*% solve(sparseSigma) %*% G1_tilde

				G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (G1_tilde_Ps_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])))
#				cat("second5\n")
		
				Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2	

		#		cat("GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]: ", GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))], "\n")		
		#		Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - (G1_tilde_Ps_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]))) %*% (G2_tilde_Ps_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])
				Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (G2_tilde_Ps_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m]))
				Phi_cond = as.matrix(Phi_cond)
			}
#			print("OKKKKKK")
			Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])
			cat("dim(Phi)", dim(Phi), "\n")
			
			}else{
				G1_tilde_G1_tilde = t(G1_tilde) %*% G1_tilde

				if(!is.null(G2_cond)){
					G2_tilde_G2_tilde = t(G2_cond_tilde) %*% G2_cond_tilde	
					G1_tilde_G2_tilde = t(G1_tilde) %*% G2_cond_tilde
					G2_tilde_G1_tilde = t(G2_cond_tilde) %*% G1_tilde

					G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (G1_tilde_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])))

					Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2

					Phi_cond = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (G2_tilde_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m]))

#					Phi_cond = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m]) - (G1_tilde_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]))) %*% (G2_tilde_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])) 
					Phi_cond = as.matrix(Phi_cond)
				}

				Phi = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m])
			

			}

        	}else{
			if(!is.null(G2_cond)){
				#G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% solve(t(G2_cond_tilde) %*% (obj$P %*% G2_cond_tilde))
				G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% getcovM(G2_cond_tilde, G2_cond_tilde, obj$P) 

				Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2

				Phi_cond = t(G1_tilde) %*% (obj$P %*% G1_tilde) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (t(G2_cond_tilde) %*% (obj$P %*% G1_tilde))
               # 		Phi_cond = t(G1_tilde) %*% (obj$P %*% G1_tilde) - (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% solve(t(G2_cond_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% (t(G2_cond_tilde) %*% (obj$P %*% G1_tilde))
		}
                	Phi = t(G1_tilde) %*% (obj$P %*% G1_tilde)
		
        	}


		#Perform the SKAT test
		if(!is.null(G2_cond)){
			#if(sum(diag(Phi_cond) < 10^-5) > 0){
			#	re_cond = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)
			#}else{
        		re_cond = SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL)
			#}
		}

#		print("HERE SKAT:::Met_SKAT_Get_Pvalue")
#		print(Phi)
#		print(Score)
#		print(r.corr)
#		print(method)
		re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)

		if(!is.null(G2_cond)){
			re$p.value.cond = re_cond$p.value
		}else{
			re$p.value.cond = NA
		}

 #       print(re)

 	 }else{ 

		#else: no marker is left for test, m = 0
    		re = list(p.value = NA, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA)
    		markerNumbyMAC = c(0,0,0,0,0,0)
  	}

  	re$IsMeta=TRUE
  	re$markerNumbyMAC = markerNumbyMAC
	print(re) 
 	return(re)

}










#obj is the rda. file output from SAIGE step 1
#G1 is genotypes for testing gene, which contains m markers
#G2_cond is G2 in the word document, genotypes for m_cond conditioning marker(s)
#G2_cond_es is beta_2_hat (effect size for the conditioning marker(s))
SAIGE_SKAT_withRatioVec  = function(G1, obj, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, G2_cond = NULL, G2_cond_es, kernel= "linear.weighted", method="optimal.adj", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
, r.corr=0, is_check_genotype=FALSE, is_dosage = TRUE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, sparseSigma = NULL, singleGClambda = 1, mu2 = NULL){
	xt <- proc.time()	
        #check the input genotype G1
        obj.noK = obj$obj.noK
        m = ncol(G1)
        n = nrow(G1)
	MAF = colMeans(G1)/2

	cat("m =", m, "\n")
	#print(MAF)
        #cat(colSums(G1),"\n")
	#print("gc1")
        #gc(verbose = T)
        id_include<-1:n
        # Added by SLEE 4/24/2017
        out.method<-SKAT:::SKAT_Check_Method(method,r.corr, n=n, m=m)
        method=out.method$method
        r.corr=out.method$r.corr
        IsMeta=out.method$IsMeta
        SKAT:::SKAT_Check_RCorr(kernel, r.corr)
#	print(ncol(G1))
#        print(is_check_genotype)
#	print(is_dosage)
#	print(missing_cutoff)
#	print(max_maf)
#	print(estimate_MAF)
	# print("gc2")
        #gc(verbose = T)


#        out.z<-SKAT:::SKAT_MAIN_Check_Z(G1, n, id_include, SetID, weights, weights.beta, impute.method, is_check_genotype, is_dosage, missing_cutoff, max_maf=1, estimate_MAF=estimate_MAF)
	
	if(is.null(weights)){
		weights <- SKAT:::Beta.Weights(MAF, weights.beta)
	}

#	rm(G1)	
	#print("gc3")
        #gc(verbose = T)

#        if(out.z$return ==1){
#                out.z$param$n.marker<-m
                #return(out.z)
#                m = 0
#        }else{

#                G1 = out.z$Z.test
#                weights = out.z$weights
#                m = ncol(G1)
#        }
#	print("gc3b")
#        gc(verbose = T)


#        cat("m", m, "\n")
        #if more than 1 marker is left, continue the test
        if(m  >  0){
		#         If G1 is sparse, change it to the sparse matrix
#                if(mean(G1) < 0.1){
#                  G1 = as(G1, "sparseMatrix")
#                }



                #cbind G1 and G2_cond to estimate the variance ratio matrix (m+m_cond) x (m+m_cond)
                if(!is.null(G2_cond)){
                        m_cond = ncol(G2_cond)
                        Zall = cbind(G1, G2_cond)
                }else{
                        Zall = G1
                }

#		print("gc3c")
#        gc(verbose = T)



		MACvec_indVec_Zall = getCateVarRatio_indVec(Zall, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)

#		 print("gc4")
#        gc(verbose = T)


		rm(Zall)

#		 print("gc5")
#        gc(verbose = T)

		GratioMatrixall = getGratioMatrix(MACvec_indVec_Zall, ratioVec)
		#print(GratioMatrixall)
                #GratioMatrixall = getGratioMatrix(Zall, ratioVec)
                #MACvec_indVec = getMACvec_indVec(G1)
		if(!is.null(G2_cond)){
			MACvec_indVec = MACvec_indVec_Zall[1:m] 
#getCateVarRatio_indVec(G1, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)
		}else{
			MACvec_indVec = MACvec_indVec_Zall
		}

#		 print("gc6")
#        gc(verbose = T)


                ##summaize the number of markers falling in each MAC category
		markerNumbyMAC = NULL
		for(i in 1:length(cateVarRatioMinMACVecExclude)){
			markerNumbyMAC = c(markerNumbyMAC, sum(MACvec_indVec == i))
		}

#		print("gc6a")
#        gc(verbose = T)

#        #         If G1 is sparse, change it to the sparse matrix
#                if(mean(G1) < 0.1){
#                  G1 = as(G1, "sparseMatrix")
#                }

                if (kernel == "linear.weighted") {
                        G1 = t(t(G1) * (weights))
                        #Z_tilde = t(t(Z_tilde) * (weights))
                }

#		print("gc6b")
#        	gc(verbose = T)


#                cat("dim(G1)", dim(G1), "\n")
#                cat("dim(obj.noK$XV)", dim(obj.noK$XV), "\n")
#                cat("dim(obj.noK$XXVX_inv)", dim(obj.noK$XXVX_inv), "\n")
#		print(object.size(G1))
#		print(object.size(obj.noK$XV))
#		print(object.size(obj.noK$XXVX_inv))	

		#G1_tilde = obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G1)
                G1_tilde = G1  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G1)
		#G1_tilde = G1 - G1temp

		#Rcpp_subtractMat_elwise(G1_tilde, G1)


#	        gc()	
#		print("gc6c")
#        	gc(verbose = T)

                Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
#                cat("dim(Score)", length(Score), "\n")
#		print(object.size(G1_tilde))	
#		print(object.size(G1))	
#		print(object.size(obj.noK))	
#		print(object.size(Score))

#		rm(G1)
#	         print("gc7")
#        gc(verbose = T)

                #compute Score test statistics after conditionining
                if(!is.null(G2_cond)){
                        G2_cond_tilde<- G2_cond  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G2_cond)
                        T2 = as.vector(t(G2_cond) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
                        #Score_cond = as.vector(t(G1) %*% matrix(obj$residuals - G2_cond_tilde%*%G2_cond_es, ncol=1)) / as.numeric(obj$theta[1])
                }



                #if no P is provides, use sparseSigma or identity Sigma
                if(is.null(obj$P)){

                        if(!is.null(sparseSigma)){
#                               cat("first\n")
                                #G1_tilde_Ps_G1_tilde = t(G1_tilde)%*% solve(sparseSigma) %*% G1_tilde
#				at <- proc.time()
#				print(at-xt)
#                                print("at-xt")
                                G1_tilde_Ps_G1_tilde = getcovM(G1_tilde, G1_tilde, sparseSigma, mu2 = mu2)
#				bt <- proc.time()
#				print(bt-at)
#				print("bt-at")
                                if(!is.null(G2_cond)){
#                               cat("second\n")
                #               G1_tilde_Ps_G1_tilde = getcovM(Z_tilde, Z_tilde, sparseSigma)
                                G2_tilde_Ps_G2_tilde = getcovM(G2_cond_tilde, G2_cond_tilde, sparseSigma, mu2 = mu2)
#                               cat("second2\n")
                                G1_tilde_Ps_G2_tilde = getcovM(G1_tilde, G2_cond_tilde, sparseSigma, mu2 = mu2)
#                               cat("second3\n")
                                #G2_tilde_Ps_G1_tilde = getcovM(G2_cond_tilde, G1_tilde, sparseSigma, mu2 = mu2)
                                G2_tilde_Ps_G1_tilde = t(G1_tilde_Ps_G2_tilde) 
#                               cat("second4\n")
                #               G2_tilde_Ps_G2_tilde = t(G2_cond_tilde)%*% solve(sparseSigma) %*% G2_cond_tilde
                #               G1_tilde_Ps_G2_tilde = t(G1_tilde)%*% solve(sparseSigma) %*% G2_cond_tilde
                #               G2_tilde_Ps_G1_tilde = t(G2_cond_tilde)%*% solve(sparseSigma) %*% G1_tilde
				#cat("G2_tilde_Ps_G2_tilde: ", G2_tilde_Ps_G2_tilde, "\n")
				#cat("G1_tilde_Ps_G2_tilde: ", G1_tilde_Ps_G2_tilde, "\n")


                                G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (G1_tilde_Ps_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])))
#                               cat("second5\n")

				#cat("G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv: ", G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv, "\n")

                                Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2
				#cat("Score_cond: , ", Score_cond, "\n")
                #               cat("GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]: ", GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))], "\n")
                #               Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - (G1_tilde_Ps_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]))) %*% (G2_tilde_Ps_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])
                                Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (G2_tilde_Ps_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m]))
                                Phi_cond = as.matrix(Phi_cond)
				#cat("Phi_cond, ", Phi_cond, "\n")
                        	}
#                       print("OKKKKKK")
			ct <- proc.time()
                        Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])
                        dt = proc.time()
                        #cat("dim(Phi)", dim(Phi), "\n")
#			cat("time: ", at, " ", bt, " ", ct, " ", dt, "\n")
                        }else{
                                G1_tilde_G1_tilde = t(G1_tilde) %*% G1_tilde

                                if(!is.null(G2_cond)){
                                        G2_tilde_G2_tilde = t(G2_cond_tilde) %*% G2_cond_tilde
                                        G1_tilde_G2_tilde = t(G1_tilde) %*% G2_cond_tilde
                                        G2_tilde_G1_tilde = t(G2_cond_tilde) %*% G1_tilde

                                        G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (G1_tilde_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])))

                                        Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2

                                        Phi_cond = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (G2_tilde_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m]))

#                                       Phi_cond = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m]) - (G1_tilde_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))]))) %*% (G2_tilde_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m]))
                                        Phi_cond = as.matrix(Phi_cond)
                                }

                                Phi = G1_tilde_G1_tilde*(GratioMatrixall[1:m,1:m])
				#print(Phi)
				#cat("dim(Phi)", dim(Phi), "\n")	

                        }

                }else{
                        if(!is.null(G2_cond)){
                                #G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% solve(t(G2_cond_tilde) %*% (obj$P %*% G2_cond_tilde))
                                G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% getcovM(G2_cond_tilde, G2_cond_tilde, obj$P)

                                Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2

                                Phi_cond = t(G1_tilde) %*% (obj$P %*% G1_tilde) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (t(G2_cond_tilde) %*% (obj$P %*% G1_tilde))
               #                Phi_cond = t(G1_tilde) %*% (obj$P %*% G1_tilde) - (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% solve(t(G2_cond_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% (t(G2_cond_tilde) %*% (obj$P %*% G1_tilde))
                }
                        Phi = t(G1_tilde) %*% (obj$P %*% G1_tilde)

                }


                #Perform the SKAT test
                if(!is.null(G2_cond)){
                        if(sum(diag(Phi_cond) < 10^-5) > 0){
                               re_cond = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)
                        }else{
                        re_cond = SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL)
                        }
                }

#               print("HERE SKAT:::Met_SKAT_Get_Pvalue")
#               print(Phi)
#               print(Score)
#               print(r.corr)
#               print(method)

		if(singleGClambda == 1){
                  Phi = Phi * singleGClambda
		  Phi = as.matrix(Phi)
                  re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
		}else{
		  re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
		  Phi = Phi * singleGClambda
		  Phi = as.matrix(Phi)
		  re_GCadj = SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
		  re$P_singlGCadj = re_GCadj$p.value
		  re$GCadjOut = re_GCadj	

		}

                if(!is.null(G2_cond)){
                        re$p.value.cond = re_cond$p.value
			re$condOut = re_cond
                }else{
                        re$p.value.cond = NA
                }
#		 print("gc8")
#        gc(verbose = T)
 #       print(re)

         }else{

                #else: no marker is left for test, m = 0
                re = list(p.value = NA, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA)
                #markerNumbyMAC = c(0,0,0,0,0,0)
                markerNumbyMAC = rep(0, length(cateVarRatioMinMACVecExclude))

        }

        re$IsMeta=TRUE
        re$markerNumbyMAC = markerNumbyMAC
	re$m = m
        print(re)
        return(re)

}












getGratioMatrix = function(MACvec_indVec, ratioVec){

	numCate = length(ratioVec)
        #cat("MACvec_indVec: ", MACvec_indVec, "\n")

        indMatrix = contr.sum(numCate, contrasts = FALSE)

        GindMatrix = NULL
        for(i in MACvec_indVec){
          GindMatrix = rbind(GindMatrix, indMatrix[i,])
        }

        GratioVec = GindMatrix %*% matrix(ratioVec, ncol=1)
        #mxm
        GratioMatrix = sqrt(GratioVec) %*% t(sqrt(GratioVec))

        return(GratioMatrix)
}




getGratioMatrix_old = function(G, ratioVec){

	MACvec = colSums(G)
        MACvec_indVec = MACvec
#        cat("MACvec: ", MACvec, "\n")
        MACvec_indVec[which(MACvec <= 1.5)] = 1
        MACvec_indVec[which(MACvec <= 2.5 & MACvec > 1.5)] = 2
        MACvec_indVec[which(MACvec <= 3.5 & MACvec > 2.5)] = 3
        MACvec_indVec[which(MACvec <= 4.5 & MACvec > 3.5)] = 4
        MACvec_indVec[which(MACvec <= 5.5 & MACvec > 4.5)] = 5
        MACvec_indVec[which(MACvec_indVec > 5.5)] = 6

#        cat("MACvec_indVec: ", MACvec_indVec, "\n")

        indMatrix = contr.sum(6, contrasts = FALSE)

        GindMatrix = NULL
        for(i in MACvec_indVec){
          GindMatrix = rbind(GindMatrix, indMatrix[i,])
        }

        #mx1 = mx6 %*% 6x1
        GratioVec = GindMatrix %*% matrix(ratioVec, ncol=1)
        #mxm
        GratioMatrix = sqrt(GratioVec) %*% t(sqrt(GratioVec))  

	#re = list(GratioMatrix = GratioMatrix
	return(GratioMatrix)
}


getcovM = function(G1, G2, sparseSigma, mu2 = NULL){

  if(!is.null(sparseSigma)){
   pcginvSigma = NULL
   for(i in 1:ncol(G2)){
     c3<-pcg(sparseSigma, G2[,i])
     pcginvSigma<-cbind(pcginvSigma, c3)
   }
   covM = as.matrix(t(G1) %*% pcginvSigma)
  }else{
      if(!is.null(mu2)){
        G2 = G2 * mu2
      }
      covM = as.matrix(t(G1) %*% G2)
  }
   return(covM)
}

getvarRatio = function(MAC, ratioVec){
	if(MAC <= 1.5 ){
		i = 1
	}else if(MAC <= 2.5 & MAC > 1.5){
		i = 2
	}else if(MAC <= 3.5 & MAC > 2.5){
		i = 3
	}else if(MAC <= 4.5 & MAC > 3.5){
		i = 4
	}else if(MAC <= 5.5 & MAC > 4.5){
		i = 5
	}else if(MAC > 5.5){
		i = 6
	}
	varRatio = ratioVec[i]
	return(varRatio)
}

getMACvec_indVec = function(Z){

	MACvec = colSums(Z)
        MACvec_indVec = MACvec
        MACvec_indVec[which(MACvec <= 1.5)] = 1
        MACvec_indVec[which(MACvec <= 2.5 & MACvec > 1.5)] = 2
        MACvec_indVec[which(MACvec <= 3.5 & MACvec > 2.5)] = 3
        MACvec_indVec[which(MACvec <= 4.5 & MACvec > 3.5)] = 4
        MACvec_indVec[which(MACvec <= 5.5 & MACvec > 4.5)] = 5
        MACvec_indVec[which(MACvec_indVec > 5.5)] = 6

	return(MACvec_indVec)	
}




getCateVarRatio_indVec = function(G, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude){

	
	if(ncol(G) > 1){
        	MACvector = colSums(G)
		
	}else{
		MACvector = NULL
		MACvector = c(MACvector, sum(as.vector(G[,1])) )
	}

	MACvector[which(MACvector > nrow(G))] = 2*nrow(G) - MACvector[which(MACvector > nrow(G))]

	#cat("MACvector: ", MACvector, "\n")
        #print(length(MACvector))
        MACvec_indVec = rep(0, length(MACvector))
	#cat("here1 MACvec_indVec: ", MACvec_indVec, "\n")		
	#cat("cateVarRatioMinMACVecExclude: ", cateVarRatioMinMACVecExclude, "\n")
	#cat("cateVarRatioMaxMACVecInclude: ", cateVarRatioMaxMACVecInclude, "\n")
  	numCate = length(cateVarRatioMinMACVecExclude)
#	cat("numCate: ", numCate, "\n")
#	cat("MACvector: ", MACvector, "\n")
    	for(i in 1:(numCate-1)){
		MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
		if(length(MACvecIndex) > 0){
			MACvec_indVec[MACvecIndex] = i
		}

    	}
#	cat("here2 MACvec_indVec: ", MACvec_indVec, "\n")		
#	cat("here2 length(cateVarRatioMaxMACVecInclude): ", length(cateVarRatioMaxMACVecInclude), "\n")		

    	if(length(cateVarRatioMaxMACVecInclude) == (numCate-1)){
		MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate])
    	}else{
		MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate] & MACvector <= cateVarRatioMaxMACVecInclude[numCate])	
    	}

	if(length(MACvecIndex) > 0){
        	MACvec_indVec[MACvecIndex] = numCate
        }
#	cat("here3 MACvec_indVec: ", MACvec_indVec, "\n")		

        return(MACvec_indVec)
}

###
getVarRatio = function(G, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec){
  if(length(ratioVec) == 1 & ncol(as.matrix(G)) == 1){
    return(ratioVec[1])
  }else{
	G = as.matrix(G)
	if(length(ratioVec) == 1){
	  ratioVec = c(ratioVec, rep(ratioVec[1], ncol(as.matrix(G))))
	}
        if(ncol(G) > 1){
                MACvector = colSums(G)

        }else{
                MACvector = NULL
                MACvector = c(MACvector, sum(as.vector(G[,1])) )
        }

        MACvector[which(MACvector > nrow(G))] = 2*nrow(G) - MACvector[which(MACvector > nrow(G))]

        #cat("MACvector: ", MACvector, "\n")
        #print(length(MACvector))
        MACvec_indVec = rep(0, length(MACvector))
        #cat("here1 MACvec_indVec: ", MACvec_indVec, "\n")
        #cat("cateVarRatioMinMACVecExclude: ", cateVarRatioMinMACVecExclude, "\n")
        #cat("cateVarRatioMaxMACVecInclude: ", cateVarRatioMaxMACVecInclude, "\n")
        numCate = length(cateVarRatioMinMACVecExclude)
#       cat("numCate: ", numCate, "\n")
#       cat("MACvector: ", MACvector, "\n")
        for(i in 1:(numCate-1)){
                MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
                if(length(MACvecIndex) > 0){
                        MACvec_indVec[MACvecIndex] = i
                }

        }
#       cat("here2 MACvec_indVec: ", MACvec_indVec, "\n")
#       cat("here2 length(cateVarRatioMaxMACVecInclude): ", length(cateVarRatioMaxMACVecInclude), "\n")

        if(length(cateVarRatioMaxMACVecInclude) == (numCate-1)){
                MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate])
        }else{
                MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate] & MACvector <= cateVarRatioMaxMACVecInclude[numCate])
        }

        if(length(MACvecIndex) > 0){
                MACvec_indVec[MACvecIndex] = numCate
        }
#       cat("here3 MACvec_indVec: ", MACvec_indVec, "\n")


	GratioMat = getGratioMatrix(MACvec_indVec, ratioVec)
        #return(MACvec_indVec)
        return(GratioMat)
  }
}

