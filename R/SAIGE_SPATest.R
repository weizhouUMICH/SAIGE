options(stringsAsFactors=F, digits=4)
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
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 0.5. The higher threshold between minMAC and minMAF will be used
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used
#' @param maxMAFforGroupTest numeric. Maximum minor allele frequency of markers to test in group test. By default 1.
#' @param minInfo numeric. Minimum imputation info of markers to test (in bgen file)
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the dosage, vcf, sav, or bgen file with NO header
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param Cutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param IsSparse logical. Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for dichotomous traits. By default, TRUE 
#' @param numLinesOutput numeric. Output results for how many marker each time.    
#' @param SAIGEOutputFile character. Path to the output file containing the SPAGMMAT test results
#' @param IsOutputAFinCaseCtrl logical. Whether to output allele frequency in cases and controls. By default, FALSE
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out option. By default, FALSE
#' @param condition. For conditional analysis. Genetic marker ids (chr:pos_ref/alt) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A, Note that currently conditional analysis is only for vcf/sav input.
#' @param sparseSigmaFile character. Path to the file containing the sparseSigma from step 1. The suffix of this file is ".mtx" 
#' @return SAIGEOutputFile
#' @export
SPAGMMATtest = function(dosageFile = "",
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
		 minMAC = 0.5, 
                 minMAF = 0,
		 maxMAFforGroupTest = 1,
        	 minInfo = 0,
                 GMMATmodelFile = "", 
                 varianceRatioFile = "", 
                 Cutoff=2, 
                 SAIGEOutputFile = "",
		 numLinesOutput = 10000, 
		 IsSparse=TRUE,
		 IsOutputAFinCaseCtrl=FALSE,
		 LOCO=FALSE,
		 condition="",
		 sparseSigmaFile="",
		 groupFile="",
		 kernel="linear.weighted",
		 method="optimal.adj",
		 weights.beta=c(1,25),
		 r.corr=0,
		 IsSingleVarinGroupTest = TRUE,
		 cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5), 
		 cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5),
		 singleGClambda = 1){


  if(groupFile == ""){
    isGroupTest = FALSE
  }else{
    cat("group-based test will be performed \n")
    if(!file.exists(groupFile)){
      stop("ERROR! groupFile ", groupFile, " does not exsit\n")
    }else{
      isGroupTest = TRUE
    }
  }


####check and read files
  #sparseSigmaFile
 cat("sparseSigmaFile: ", sparseSigmaFile, "\n")
  if(sparseSigmaFile == ""){
    sparseSigma = NULL
  }else{
    cat("sparse kinship matrix is going to be used\n")
    if(!file.exists(sparseSigmaFile)){
      stop("ERROR! sparseSigmaFile ", sparseSigmaFile, " does not exsit\n")
    }else{
      sparseSigma = Matrix:::readMM(sparseSigmaFile)
 cat("sparseSigmaFile: ", sparseSigmaFile, "\n")
    }
  }


  #output file
  if(file.exists(SAIGEOutputFile)){file.remove(SAIGEOutputFile)}

  if(!file.exists(SAIGEOutputFile)){
    file.create(SAIGEOutputFile, showWarnings = TRUE)
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

    if(!LOCO | is.null(obj.glmm.null$LOCO)){
      obj.glmm.null$LOCO = FALSE
      cat("obj.glmm.null$LOCO: ", obj.glmm.null$LOCO, "\n")
      cat("Leave-one-chromosome-out option is not applied\n")
    }
 
  }


  #allowing for categorical variance ratio
  if(!file.exists(varianceRatioFile)){
    stop("ERROR! varianceRatioFile ", varianceRatioFile, " does not exsit\n")
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
    #bgenFileIndex is not necessary if there is no query
    #if(!file.exists(bgenFileIndex)){
    #stop("ERROR! bgenFileIndex ", bgenFileIndex, " does not exsit\n")
    #}
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
      stop("ERROR! savFile ", savFile, " does not exsit\n")
    }else{
      vcfFile = savFile	
    }

    if(!file.exists(savFileIndex)){
      stop("ERROR! savFileIndex ", savFileIndex, " does not exsit\n")
    }else{
      vcfFileIndex = savFileIndex
    }	
    dosageFileType = "vcf"
  }


  ##############START TEST########################
  startTime = as.numeric(Sys.time())  # start time of the SPAGMMAT tests
  cat("Analysis started at ", startTime, "Seconds\n")

#  if(file.exists(SAIGEOutputFile)){file.remove(SAIGEOutputFile)}

  if(!isGroupTest){

    if (dosageFileType == "bgen" | dosageFileType == "vcf"){
      dosageFilecolnamesSkip = c("CHR","POS","SNPID","Allele1","Allele2", "AC_Allele2", "AF_Allele2", "imputationInfo")
    }else{
      dosageFilecolnamesSkip = c(dosageFilecolnamesSkip, "AC", "AF")
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
      Gx_cond = getGenoOfGene_vcf(conditionlist, minInfo)
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
    }else{    #if(cntMarker > 0){
      cat("WARNING: conditioning markers are not found in the provided dosage file \n")
      #print(dim(dosage_cond))
      isCondition = FALSE
      dosage_cond = NULL
    }

  }else{#end of if(isCondition){
    dosage_cond = NULL
  }

if(traitType == "binary"){
  cat("It is a binary trait\n")
  if(!isGroupTest){
    resultHeader = c(dosageFilecolnamesSkip, "N", "BETA", "SE", "Tstat", "p.value", "p.value.NA", "Is.SPA.converge","varT","varTstar")
    if(IsOutputAFinCaseCtrl){
      resultHeader = c(resultHeader, "AF.Cases", "AF.Controls")
    }
    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
  }

  if(Cutoff < 10^-2){
    Cutoff=10^-2
  }

  y = obj.glm.null$y
    #OUT = NULL
    #numPassMarker = 0
    #NSparse=0
    #mth = 0
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
 if(!isGroupTest){
   if(!isCondition){
     resultHeader = c(dosageFilecolnamesSkip,  "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar")
   }else{
     resultHeader = c(dosageFilecolnamesSkip,  "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar","Tstat_cond", "p.value_cond", "varT_cond", "BETA_cond", "SE_cond" )
   }
   write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
 }

    #OUT = NULL
    #numPassMarker = 0
    #mth = 0
    #sampleIndex = sampleIndex - 1
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
#  dosage_cond_tilde = dosage_cond - (obj.noK$XVX_inv_XV)%*%dosage_cond
  Mcond = ncol(dosage_cond_tilde)
      #G0_tilde = G0 - (obj.noK$XVX_inv_XV)%*%G0
      #dosage_cond_tilde = cbind(G0_tilde, dosage_cond_tilde)
  covM = matrix(0,nrow=Mcond+1, ncol = Mcond+1)
  covMsub = getcovM(dosage_cond_tilde, dosage_cond_tilde, sparseSigma)
  covM[2:(Mcond+1), 2:(Mcond+1)] = covMsub


  MACvec_indVec_cond = getCateVarRatio_indVec(dosage_cond, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)
  GratioMatrix_cond = getGratioMatrix(MACvec_indVec_cond, ratioVec)

 print("covMsub")
 print(covMsub)
 print("GratioMatrix_cond")
 print(GratioMatrix_cond)

 print(covMsub%*%GratioMatrix_cond)
#  GratioMatrix_cond = getGratioMatrix(dosage_cond, ratioVec)
  G2tilde_P_G2tilde_inv = solve(covMsub %*% GratioMatrix_cond)

}else{# end of if(isCondition)
  OUT_cond = NULL
  G2tilde_P_G2tilde_inv = NULL
}


#determine minimum MAF for markers to be tested
  if(minMAC == 0){minMAC = 1} ##01-19-2018
  cat("minMAC: ",minMAC,"\n")
  cat("minMAF: ",minMAF,"\n")
  minMAFBasedOnMAC = minMAC/(2*N) 
  testMinMAF = max(minMAFBasedOnMAC, minMAF) 
  cat("Minimum MAF of markers to be testd is ", testMinMAF, "\n")


  ##############START TEST########################
  startTime = as.numeric(Sys.time())  # start time of the SPAGMMAT tests
  cat("Analysis started at ", startTime, "Seconds\n")

if(!isGroupTest){

  isVariant = TRUE

  if(dosageFileType == "plain"){
    Mtest = setgenoTest_plainDosage(dosageFile, dosageFileNrowSkip, dosageFileNcolSkip)
    if(Mtest == 0){isVariant = FALSE}
    SetSampleIdx_plainDosage(sampleIndex, N)

  }else if (dosageFileType == "bgen"){
    if(idstoExcludeFile != ""){
      idsExclude = data.table:::fread(idstoExcludeFile, header=F,sep=" ", stringsAsFactors=FALSE, colClasses=c("character"))
      idsExclude = data.frame(idsExclude)
      ids_to_exclude = as.character(as.vector(idsExclude[,1]))
    }else{
      ids_to_exclude = as.character(vector())
    }	

    if(idstoIncludeFile != ""){
      idsInclude = data.table:::fread(idstoIncludeFile, header=F, sep=" ", stringsAsFactors=FALSE, colClasses=c("character"))
      idsInclude = data.frame(idsInclude)
      ids_to_include = as.character(as.vector(idsInclude[,1]))
    }else{
      ids_to_include = as.character(vector())
    }

    if(rangestoExcludeFile != ""){
      rangesExclude = data.table:::fread(rangestoExcludeFile, header=F)
      ranges_to_exclude = data.frame(rangesExclude)
      colnames(ranges_to_exclude) = c("chromosome","start","end")  
    }else{
      ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
    }

    if(rangestoIncludeFile != ""){
      rangesInclude = data.table:::fread(rangestoIncludeFile, header=F)
      ranges_to_include = data.frame(rangesInclude)
      colnames(ranges_to_include) = c("chromosome","start","end")
    }else{
      ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
    }

    Mtest = setgenoTest_bgenDosage(bgenFile,bgenFileIndex, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
    if(Mtest == 0){isVariant = FALSE}
    isQuery = getQueryStatus()
    SetSampleIdx(sampleIndex, N)
  }else if(dosageFileType == "vcf"){

    setgenoTest_vcfDosage(vcfFile,vcfFileIndex,vcfField,ids_to_exclude_vcf = idstoExcludeFile, ids_to_include_vcf = idstoIncludeFile, chrom, start, end)
    #setTestField(vcfField)
    isVariant = getGenoOfnthVar_vcfDosage_pre()
    SetSampleIdx_vcfDosage(sampleIndex, N)
  }

if(!isVariant){
  cat("isVariant: ", isVariant, "\n") 
  stop("ERROR! opening dosage file failed!!\n")
}

write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
OUT = NULL
numPassMarker = 0
mth = 0

while(isVariant){
    mth = mth + 1
    if(dosageFileType == "plain"){
      G0 = getGenoOfnthVar_plainDosage(mth, dosageFileNrowSkip, dosageFileNcolSkip)
      markerInfo = 1 ##markerInfo is nor provided
      AC = sum(G0)
      AF = AC/(2*N)
      rowHeader=getrowHeaderVec_plainDosage()
      rowHeader = c(rowHeader, AC, AF, markerInfo)
      if(indChromCheck){      
        CHR = rowHeader[which(dosageFilecolnamesSkip == dosageFileChrCol)]
        #cat("CHR is here ", CHR , "\n")
      }

      if(Mtest == mth){isVariant = FALSE}

    }else if (dosageFileType == "bgen"){
      if(isQuery){
        Gx = getDosage_bgen_withquery()
      }else{
        Gx = getDosage_bgen_noquery()
      }
      markerInfo = getMarkerInfo()
      G0 = Gx$dosages
      AC = Gx$variants$AC
      AF = Gx$variants$AF
      Gx$variants$markerInfo = markerInfo
      rowHeader=as.vector(unlist(Gx$variants))
      if(indChromCheck){
	CHR = Gx$variants$chromosome
      }	

      if(Mtest == mth){isVariant = FALSE}
    }else if(dosageFileType == "vcf"){
      markerInfo = 1 ##markerInfo is nor provided
      Gx = getGenoOfnthVar_vcfDosage(mth)
      G0 = Gx$dosages
      AC = Gx$variants$AC
      AF = Gx$variants$AF
      markerInfo = Gx$variants$markerInfo
      rowHeader=as.vector(unlist(Gx$variants))
      if(indChromCheck){
        CHR = Gx$variants$chromosome
      }
      isVariant = getGenoOfnthVar_vcfDosage_pre()
    }

    MAC = AC
    MAF = AF
    if(AF > 0.5){
      MAC = 2*N - AC
      MAF = 1 - AF
    }


#cat("MAF is ", MAF, "\n")
#cat("testMinMAF is ", testMinMAF, "\n")
#cat("markerInfo is ", markerInfo, "\n")
#cat("minInfo is ", minInfo, "\n")
if(MAF >= testMinMAF & markerInfo >= minInfo){
   numPassMarker = numPassMarker + 1


    varRatio = getvarRatio(MAC, ratioVec)

    if(indChromCheck){
      CHR = as.character(CHR)
      CHRv2 = as.numeric(gsub("[^0-9.]", "", CHR))
      #cat("CHRv2 is ", CHRv2, "\n") 
      #cat("CHR is ", CHR, "\n") 
    
      if(obj.glmm.null$LOCOResult[[CHRv2]]$isLOCO){
        mu = obj.glmm.null$LOCOResult[[CHRv2]]$fitted.values
        mu.a<-as.vector(mu)
	if(traitType == "binary"){
          mu2.a<-mu.a *(1-mu.a)
	}
      }else{
        mu = obj.glmm.null$fitted.values
        mu.a<-as.vector(mu)
	if(traitType == "binary"){
          mu2.a<-mu.a *(1-mu.a)
	}
      }
      obj.noK$XVX = t(obj.noK$X1) %*% (obj.noK$X1 * mu2.a)
      obj.noK$XVX_inv_XV = obj.noK$XXVX_inv * obj.noK$V
      obj.noK$S_a = colSums(obj.noK$X1 * (y - mu.a))

    }

  
##conditional analysis
   if(isCondition){	
      #dosage_cond_tilde = dosage_cond - (obj.noK$XVX_inv_XV)%*%dosage_cond
#      G0_tilde = G0 - (obj.noK$XVX_inv_XV)%*%G0
      Gall = cbind(G0, dosage_cond)
if(FALSE){
      G0test = G0	
      if(AF > 0.5){G0test = 2-G0}
      G0condtest=dosage_cond
      for(i in 1:ncol(dosage_cond)){
		AFcond = mean(dosage_cond[,i])/2
			
		if(AFcond > 0.5){
		G0condtest[,i] = 2-dosage_cond[,i]
		MAFcond = 1-AFcond
		G0condtest[,i] = (G0condtest[,i]-MAFcond)/sqrt(MAFcond*2*(1-MAFcond))

		}else{
			MAFcond = AFcond
			G0condtest[,i] = dosage_cond[,i]
			G0condtest[,i] = (G0condtest[,i]-MAFcond)/sqrt(MAFcond*2*(1-MAFcond))

		}
	}

      Gcorr = t((G0test-MAF)/(sqrt(2*MAF*(1-MAF))))%*%G0condtest/N
      print("Gcorr")
      print(Gcorr)
}

      MACvec_indVec_Gall = getCateVarRatio_indVec(Gall, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)
      GratioMatrixall = getGratioMatrix(MACvec_indVec_Gall, ratioVec)
	

      #GratioMatrixall = getGratioMatrix(Gall, ratioVec)
      G0_tilde = G0 - obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G0)
      #dosage_cond_tilde = cbind(G0_tilde, dosage_cond_tilde)	
      #covM = matrix(0,nrow=ncol(dosage_cond_tilde), ncol = ncol(dosage_cond_tilde))	
      #print(dim(G0_tilde))
      #print(dim(dosage_cond_tilde))	
      #covM[1,2:ncol(covM)] = (t(mu.a*(1-mu.a)*G0_tilde)) %*% dosage_cond_tilde			
      covM[1,2:ncol(covM)] = getcovM(G0_tilde, dosage_cond_tilde,sparseSigma)			
      G1tilde_P_G2tilde = covM[1,c(2:ncol(covM))]*(GratioMatrixall[1,c(2:ncol(covM))])

   }else{ #end of if(isCondition)
      G1tilde_P_G2tilde = NULL
      GratioMatrixall = NULL
   }

  if(traitType == "binary"){
    if(!IsOutputAFinCaseCtrl){

      OUT = rbind(OUT, scoreTest_SAIGE_binaryTrait(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y, varRatio, Cutoff, rowHeader) )	
      #OUT = rbind(OUT, scoreTest_SAIGE_binaryTrait_cond(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y, varRatio, Cutoff, rowHeader, covM, OUT_cond, covariateVec))
    }else{
      AFCase = sum(G0[y1Index])/(2*NCase)
      AFCtrl = sum(G0[y0Index])/(2*NCtrl)
      OUT = rbind(OUT, c(scoreTest_SAIGE_binaryTrait(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y, varRatio, Cutoff, rowHeader),AFCase, AFCtrl))	
      #OUT = rbind(OUT, c(scoreTest_SAIGE_binaryTrait_cond(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y, varRatio, Cutoff, rowHeader, covM, OUT_cond, covariateVec),AFCase, AFCtrl))
    }
   }else if(traitType == "quantitative"){

      out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.noK, AC, AF, y, mu, varRatio, tauVec, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv)

      if(!isCondition){
        OUT = rbind(OUT, c(rowHeader, N, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2))
      }else{
        OUT = rbind(OUT, c(rowHeader, N, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2, out1$Tstat_c,  out1$p.value.c, out1$var1_c, out1$BETA_c, out1$SE_c))
      }


   }


 } #end of the if(MAF >= bgenMinMaf & markerInfo >= bgenMinInfo)
      #if(mth %% 100000 == 0 | mth == Mtest){
 if(mth %% numLinesOutput == 0 | !isVariant){
    ptm <- proc.time()
    print(ptm)
    print(mth)
    cat("numPassMarker: ", numPassMarker, "\n")
    OUT = as.data.frame(OUT)
    write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
    OUT = NULL
  }
  #if(mth == 100){break}
} ####end of while(isVariant)

}else{ #end if(!isGroupTest){
  OUT_single = NULL
  if(IsSingleVarinGroupTest){
    SAIGEOutputFile_single = paste0(SAIGEOutputFile, "_single")
    headerline = c("markerID", "AC", "AF", "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar")
    write(headerline,file = SAIGEOutputFile_single, ncolumns = length(headerline))
  }


  if(!isCondition){
    if(dosageFileType == "plain"){
      isCondition = FALSE
    }else if(dosageFileType == "bgen"){
      SetSampleIdx(sampleIndex, N)
    }else if(dosageFileType == "vcf"){
      setMAFcutoffs(testMinMAF, maxMAFforGroupTest)
      isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
      SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N)
    }
  }
  
  if(traitType == "quantitative"){
    OUT = NULL
    cat("It is a quantitative trait\n")
    mth = 0
    if(!isCondition){
      resultHeader = c("Gene", "Pvalue", "N_MAC1","N_MAC2","N_MAC3","N_MAC4","N_MAC5","N_MACgt5","markerIDs","markerAFs")
    }else{
      resultHeader = c("Gene","Pvalue","Pvalue_cond","N_MAC1","N_MAC2","N_MAC3","N_MAC4","N_MAC5","N_MACgt5","markerIDs","markerAFs")
    }
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
          Gx = getGenoOfGene_vcf(marker_group_line, minInfo)
        }else if(dosageFileType == "bgen"){
          Gx = getGenoOfGene_bgen(bgenFile,bgenFileIndex, marker_group_line, testMinMAF, maxMAFforGroupTest, minInfo)
        }

        G0 = Gx$dosages
        cntMarker = Gx$cnt
        cat("cntMarker: ", cntMarker, "\n")
#       cat("markerIDs: ", Gx$markerIDs, "\n")
#       cat("G0: ", G0, "\n")
        if(cntMarker > 0){
         Gmat = matrix(G0, byrow=F, ncol = cntMarker)
#        cat("dim(Gmat): ", dim(Gmat), "\n")
#        cat("Gmat[,1]: ", Gmat[,1], "\n")
#        cat("ratioVec: ", ratioVec, "\n")
#        cat("Gmat[,2]: ", Gmat[,2], "\n")
        #  saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec)
#        cat("colSums(Gmat): ", colSums(Gmat), "\n")
         #saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec, Z_cond=dosage_cond, Z_cond_es=OUT_cond[,1], max_maf = maxMAF, sparseSigma = sparseSigma, method = "optimal.adj")
	 #print(is.null(sparseSigma))
	 #print("is.null(sparseSigma)")
         #saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null, ratioVec, G2_cond=dosage_cond, G2_cond_es=OUT_cond[,1], kernel=kernel, method = method, weights.beta = weights.beta, r.corr = r.corr, max_maf = maxMAFforGroupTest, sparseSigma = sparseSigma)
         saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null,  cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude,ratioVec, G2_cond=dosage_cond, G2_cond_es=OUT_cond[,1], kernel=kernel, method = method, weights.beta = weights.beta, r.corr = r.corr, max_maf = maxMAFforGroupTest, sparseSigma = sparseSigma, singleGClambda = singleGClambda)

          cat("saigeskatTest$p.value: ", saigeskatTest$p.value, "\n")

	if(IsSingleVarinGroupTest){
		for(nc in 1:ncol(Gmat)){
			G0_single = Gmat[,nc]
			AC = sum(G0_single)
			AF = AC/(2*length(G0_single))
			MAC = min(AC, 2*length(G0_single)-AC)
			varRatio = getvarRatio(MAC, ratioVec)
			varRatio = varRatio * singleGClambda			
		      	out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0_single, obj.noK, AC, AF, y, mu, varRatio, tauVec, sparseSigma=sparseSigma)
			OUT_single = rbind(OUT_single, c((Gx$markerIDs)[nc], AC, (Gx$markerAFs)[nc], N, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2))
		}

	}

	if(isCondition){
          OUT = rbind(OUT, c(geneID, saigeskatTest$p.value, saigeskatTest$p.value.cond, saigeskatTest$markerNumbyMAC, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))
	}else{
	  OUT = rbind(OUT, c(geneID, saigeskatTest$p.value, saigeskatTest$markerNumbyMAC, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))
	}
          mth = mth + 1
          if(mth %% numLinesOutput == 0){
            ptm <- proc.time()
            print(ptm)
            print(mth)
            OUT = as.data.frame(OUT)
            write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
            OUT = NULL
		
	    OUT_single = as.data.frame(OUT_single)
	    write.table(OUT_single, SAIGEOutputFile_single, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)	
	    OUT_single = NULL	
          }
        }
      }#end of else for if(length(line) == 0 )
    } # end of while ( TRUE ) {

    if(!is.null(OUT)){
      OUT = as.data.frame(OUT)
      write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
      OUT = NULL

      OUT_single = as.data.frame(OUT_single)
      write.table(OUT_single, SAIGEOutputFile_single, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
      OUT_single = NULL
	

    }
  }else{
    stop("ERROR! The type of the trait has to be quantitative\n")
  }  
  
}#if(groupTest)


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
#  cat("HERE\n")
#  cat("AC: ",AC,"\n")
#  cat("AF: ",AF,"\n")
  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }
  maf = min(AF, 1-AF)
#  cat("HERE2\n")
  if(maf < 0.05){
#  cat("HERE2a\n")
    idx_no0<-which(G0>0)
    #cat("length(idx_no0): ", length(idx_no0), "\n")
    #cat("maf: ", maf, "\n")
    g1<-G0[idx_no0]/sqrt(AC2)
    A1<-obj.noK$XVX_inv_XV[idx_no0,]
    X1<-obj.noK$X1[idx_no0,]
    mu1<-mu[idx_no0]
    y1<-obj.noK$y[idx_no0]
 
    noCov = FALSE
    if(dim(obj.noK$X1)[2] == 1){
     noCov = TRUE
    }

## V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv
    if(length(idx_no0) > 1){
      Z = t(A1) %*% g1
      B<-X1 %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      if(!noCov){
        S_a2 = obj.noK$S_a - colSums(X1 * (y1 - mu1))
      }else{
        S_a2 = obj.noK$S_a - crossprod(X1, y1 - mu1)
      }
      #S_a2 = obj.noK$S_a - colSums(X1 * (y1 - mu1))
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
#  cat("HERE2b\n")
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G1 is X adjusted
    g = G/sqrt(AC2)
    q = innerProduct(g, y)
    m1 = innerProduct(mu, g)
    var2 = innerProduct(g, g)
    var1 = var2 * varRatio
    Tstat = (q-m1)/tauVec[1]
  }
  if(AF > 0.5){
    Tstat = (-1)*Tstat
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

  noCov = FALSE
  if(dim(obj.null$X1)[2] == 1){
   noCov = TRUE 
  }

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

    if(!noCov){
      S_a2 = obj.null$S_a - colSums(X1 * (y1 - mu1))
    }else{
      S_a2 = obj.null$S_a - crossprod(X1, y1 - mu1)
    }

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


scoreTest_SAIGE_binaryTrait_cond=function(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y,varRatio, Cutoff, rowHeader, covM, OUT_cond,covariateVec){
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

  if(Run1){
    G0 = matrix(G0, ncol = 1)
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
    g = G/sqrt(AC2)
    NAset = which(G0==0)
    out1 = scoreTest_SPAGMMAT_binaryTrait_cond(g, AC2, NAset, y, mu.a, varRatio, Cutoff = Cutoff, covM, OUT_cond, covariateVec)
    if(AF > 0.5){
      out1$BETA = (-1)*out1$BETA
      out1$Tstat = (-1)*out1$Tstat
    }
    out1 = unlist(out1)

    #outVec = c(rowHeader, N, out1["BETA"], out1["SE"], out1["Tstat"], out1["p.value"], out1["p.value.NA"], out1["Is.converge"], out1["var1"], out1["var2"])
    outVec = c(rowHeader, N, out1["BETA"], out1["SE"], out1["Tstat"], out1["p.value"], out1["p.value.NA"], out1["Is.converge"], out1["var1"], out1["var2"], out1["p.value1c"], out1["p.value.NA1c"], out1["BETA1c"], out1["SE1c"], out1["Tstat1c"])

   }
  return(outVec)
}


scoreTest_SPAGMMAT_binaryTrait_cond=function(g, AC, NAset, y, mu, varRatio, Cutoff, covM, OUT_cond, covariateVec){
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


  ##condition
  print(OUT_cond[,1])
  print(covM[1,2:ncol(covM)]) 
  Tstat1c = Tstat*sqrt(AC) - innerProduct(as.vector(OUT_cond[,1]), covariateVec)

  if(ncol(covM) > 2){
  G2_tildeG2_tilde1_v1  = t(covM[2:ncol(covM),2:ncol(covM)])
  diag(G2_tildeG2_tilde1_v1) = 0
  G2_tildeG2_tilde = G2_tildeG2_tilde1_v1 + covM[2:ncol(covM),2:ncol(covM)]

  print(covM)
#  print(dim(covM[1,2:ncol(covM)]))
#  print(dim(solve(G2_tildeG2_tilde)))
#  print(dim(t(covM[1,2:ncol(covM)])))

  covMsub = matrix(as.vector(covM[1,2:ncol(covM)]), nrow=1)
#  var1c = var1*AC - varRatio*covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub)
  var1c = varRatio*covM[1,1] - varRatio*covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub)
  }else{
    G2_tildeG2_tilde = covM[2:ncol(covM),2:ncol(covM)]
    covMsub = matrix(as.vector(covM[1,2:ncol(covM)]), nrow=1)
    print(covM)
    cat("covMsub ", covMsub, "\n")
    cat("G2_tildeG2_tilde ", G2_tildeG2_tilde, "\n")
    cat("var1 ", var1, "\n")
    cat("var1*AC ", var1*sqrt(AC), "\n")
    cat("varRatio: ", varRatio, "\n")
    cat("t(covMsub): ", t(covMsub), "\n")
    cat("var1*AC: ", var1*AC, "\n")
    cat("covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub): ", covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub), "\n")
#    var1c = var1*AC - varRatio*(covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub))
    var1c = varRatio*covM[1,1] - varRatio*(covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub))
  }
  #var1c = var1 - varRatio*(covM[1,2:ncol(covM)])%*% solve(G2_tildeG2_tilde) %*% t(covM[1,2:ncol(covM)])
  cat("var1c: ", var1c, "\n")

if(var1c > 10^-5){
  qtilde1c = ((Tstat1c)/sqrt(AC))/sqrt(var1c/AC) * sqrt(var2/AC) + m1
  
  if(length(NAset)/length(g) < 0.5){
    print("OK")
    out1c = SPAtest:::Saddle_Prob(q=qtilde1c, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
  }else{
    print("OK1")
    out1c = SPAtest:::Saddle_Prob_fast(q=qtilde1c, g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8)
  }

  out1c = c(out1c, var1c = var1c)

  logOR1c = (Tstat1c/var1c)/sqrt(AC)
  SE1c = abs(logOR1c/qnorm(out1c$p.value/2))
  out1c = c(out1c, BETA1c = logOR1c, SE1c = SE1c, Tstat1c = Tstat1c)
  out1 = c(out1, p.value1c = out1c$p.value, p.value.NA1c = out1c$p.value.NA, BETA1c = out1c$BETA1c, SE1c=out1c$SE1c, Tstat1c = out1c$Tstat1c)
}else{ #end of if(var1c > 0){
  out1 = c(out1, p.value1c = 1, p.value.NA1c = 1, BETA1c = 0, SE1c=0, Tstat1c = NA)
}

  #out1 = c(out1, p.value1c = out1c$p.value, p.value.NA1c = out1c$p.value.NA, BETA1c = out1c$BETA1c, SE1c=out1c$SE1c, Tstat1c = out1c$Tstat1c)
  return(out1)
}



scoreTest_SAIGE_quantitativeTrait_sparseSigma=function(G0, obj.noK, AC, AF, y, mu, varRatio, tauVec, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

#  cat("HERE\n")
#  cat("AC: ",AC,"\n")
#  cat("AF: ",AF,"\n")
  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }
  maf = min(AF, 1-AF)
#  cat("HERE2\n")


if(maf < 0.05){
#  cat("HERE2a\n")
    idx_no0<-which(G0>0)
    #cat("length(idx_no0): ", length(idx_no0), "\n")
    #cat("maf: ", maf, "\n")
#    g1<-G0[idx_no0]/sqrt(AC2)
    g1<-G0[idx_no0]
    A1<-obj.noK$XVX_inv_XV[idx_no0,]
    X1<-obj.noK$X1[idx_no0,]
    mu1<-mu[idx_no0]
    y1<-obj.noK$y[idx_no0]

    noCov = FALSE
    if(dim(obj.noK$X1)[2] == 1){
     noCov = TRUE
    }

## V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv
    if(length(idx_no0) > 1){
      Z = t(A1) %*% g1
      B<-X1 %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      if(!noCov){
        S_a2 = obj.noK$S_a - colSums(X1 * (y1 - mu1))
      }else{
        S_a2 = obj.noK$S_a - crossprod(X1, y1 - mu1)
      }
      #S_a2 = obj.noK$S_a - colSums(X1 * (y1 - mu1))
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
#    cat("HERE2b\n")
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G1 is X adjusted
#    g = G/sqrt(AC2)
    g = G
    q = innerProduct(g, y)
    m1 = innerProduct(mu, g)
    var2 = innerProduct(g, g)
    var1 = var2 * varRatio
    Tstat = (q-m1)/tauVec[1]
  }

if(!is.null(sparseSigma)){
  XVG0 = eigenMapMatMult(obj.noK$XV, G0)
  g = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G1 is X adjusted
#  print(G0)
#  print(g)
  pcginvSigma<-pcg(sparseSigma, g)
#  print(pcginvSigma)
#  pcginvSigma = solve(sparseSigma) %*% g
#  print(solve(sparseSigma))
  var2 = as.matrix(t(g) %*% pcginvSigma) 
#  cat("var2 is ", var2, "\n")

  var1 = var2 * varRatio 

}
#cat("Tstat is ", Tstat, "\n")

if(isCondition){

  T2stat = OUT_cond[,2]
  #m_all = nrow(GratioMatrixall)


#  G1tilde_P_G2tilde = (covM[1,c(2:m_all)]*(GratioMatrixall[1,c(2:m_all)]))
#  G2tilde_P_G2tilde_inv = solve(covM[c(2:m_all),c(2:m_all)]*(GratioMatrixall[c(2:m_all),c(2:m_all)]))

 
  #Tstat_c = Tstat - covM[1,c(2:m_all)] %*% (solve(covM[c(2:m_all),c(2:m_all)])) %*% T2stat
  Tstat_c = Tstat - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% T2stat
  cat("G1tilde_P_G2tilde: ", G1tilde_P_G2tilde, "\n")
  cat("G2tilde_P_G2tilde_inv: ", G2tilde_P_G2tilde_inv, "\n")
  cat("T2stat: ", T2stat, "\n")

  #var1_c = var1 - (covM[1,c(2:m_all)]*(GratioMatrixall[1,c(2:m_all)])) %*% solve(covM[c(2:m_all),c(2:m_all)]*(GratioMatrixall[c(2:m_all),c(2:m_all)])) %*% (t(covM[1,c(2:m_all)]) * t(GratioMatrixall[1,c(2:m_all)]))
  var1_c = var1 - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% t(G1tilde_P_G2tilde)

#(covM[1,c(2:m_all)]*(GratioMatrixall[1,c(2:m_all)])) %*% solve(covM[c(2:m_all),c(2:m_all)]*(GratioMatrixall[c(2:m_all),c(2:m_all)])) %*% (t(covM[1,c(2:m_all)]) * t(GratioMatrixall[1,c(2:m_all)]))

}


if(AF > 0.5){
    Tstat = (-1)*Tstat
    if(isCondition){
      Tstat_c = (-1)*Tstat_c
    }	
}

if(var1 < 1*10^-5){
  p.value = 1
  BETA = NA
  SE = NA
}else{
  p.value = pchisq(Tstat^2/var1, lower.tail = FALSE, df=1)
#  BETA = (Tstat/var1)/sqrt(AC2)
  BETA = (Tstat/var1)
  SE = abs(BETA/qnorm(p.value/2))
}


  if(isCondition){
    if(var1_c < 1*10^-5){
      p.value.c = 1
      BETA_c = NA
      SE_c = NA
    }else{
      p.value.c = pchisq(Tstat_c^2/var1_c, lower.tail = FALSE, df=1)
#    BETA_c = (Tstat_c/var1_c)/sqrt(AC2)
      BETA_c = (Tstat_c/var1_c)
      SE_c = abs(BETA_c/qnorm(p.value.c/2))
    }
  }
  
  if(isCondition){
    out1 = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2, BETA_c = BETA_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = p.value.c, var1_c = var1_c)
  }else{
    out1 = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
  }


  return(out1)
}


