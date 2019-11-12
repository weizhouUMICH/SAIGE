#' Run single variant score tests with SPA based on the logistic mixed model.
#'
#' @param bgenFile character. Path to bgen file. Currently version 1.2 with 8 bit compression is supported
#' @param bgenFileIndex character. Path to the .bgi file (index of the bgen file)
#' @param vcfFile character. Path to vcf file
#' @param vcfFileIndex character. Path to index for vcf file by tabix, ".tbi" by "tabix -p vcf file.vcf.gz"
#' @param vcfField character. genotype field in vcf file to use. "DS" for dosages or "GT" for genotypes. By default, "DS".
#' @param savFile character. Path to sav file
#' @param savFileIndex character. Path to index for sav file .s1r
#' @param idstoExcludeFile character. Path to the file containing variant ids to be excluded from the bgen file. The file does not have a header and each line is for a marker ID.
#' @param idstoIncludeFile character. Path to the file containing variant ids to be included from the bgen file. The file does not have a header and each line is for a marker ID.
#' @param rangestoExcludeFile character. Path to the file containing genome regions to be excluded from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header 
#' @param rangestoIncludeFile character. Path to the file containing genome regions to be included from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header 
#' @param chrom character. string for the chromosome to include from vcf file. Required for vcf file. Note: the string needs to exactly match the chromosome string in the vcf/sav file. For example, "1" does not match "chr1". If LOCO is specified, providing chrom will save computation cost
#' @param start numeric. start genome position to include from vcf file. By default, 1 
#' @param end numeric. end genome position to include from vcf file. By default, 250000000
#' @param IsDropMissingDosages logical. whether to drop missing dosages (TRUE) or to mean impute missing dosages (FALSE). By default, FALSE. This option only works for bgen, vcf, and sav input.  
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 0.5. The higher threshold between minMAC and minMAF will be used
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used
#' @param maxMAFforGroupTest numeric. Maximum minor allele frequency of markers to test in group test. By default 0.5.
#' @param minInfo numeric. Minimum imputation info of markers to test. By default, 0. This option only works for bgen, vcf, and sav input
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the dosage, vcf, sav, or bgen file with NO header
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param SPAcutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param SAIGEOutputFile character. Path to the output file containing assoc test results
#' @param numLinesOutput numeric. Number of  markers to be output each time. By default, 10000   
#' @param IsSparse logical. Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for dichotomous traits. By default, TRUE 
#' @param IsOutputAFinCaseCtrl logical. Whether to output allele frequency in cases and controls. By default, FALSE
#' @param IsOutputNinCaseCtrl logical. Whether to output sample sizes in cases and controls. By default, FALSE
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out option. By default, FALSE
#' @param condition character. For conditional analysis. Genetic marker ids (chr:pos_ref/alt if sav/vcf dosage input , marker id if bgen input) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A, Note that currently conditional analysis is only for bgen,vcf,sav input.
#' @param sparseSigmaFile character. Path to the file containing the sparseSigma from step 1. The suffix of this file is ".mtx". 
#' @param groupFile character. Path to the file containing the group information for gene-based tests. Each line is for one gene/set of variants. The first element is for gene/set name. The rest of the line is for variant ids included in this gene/set. For vcf/sav, the genetic marker ids are in the format chr:pos_ref/alt. For bgen, the genetic marker ids should match the ids in the bgen file. Each element in the line is seperated by tab. 
#' @param kernel character. For gene-based test. By default, "linear.weighted". More options can be seen in the SKAT library 
#' @param method character. method for gene-based test p-values. By default, "optimal.adj". More options can be seen in the SKAT library
#' @param weights.beta.rare vector of numeric. parameters for the beta distribution to weight genetic markers with MAF <= weightMAFcutoff in gene-based tests.By default, "c(1,25)". More options can be seen in the SKAT library
#' @param weights.beta.common vector of numeric. parameters for the beta distribution to weight genetic markers with MAF > weightMAFcutoff in gene-based tests.By default, "c(0.5,0.5)". More options can be seen in the SKAT library 
#' @param weightMAFcutoff numeric. Between 0 and 0.5. See document above for weights.beta.rare and weights.beta.common. By default, 0.01
#' @param r.corr numeric. bewteen 0 and 1. parameters for gene-based tests.  By default, 0.  More options can be seen in the SKAT library
#' @param IsSingleVarinGroupTest logical. Whether to perform single-variant assoc tests for genetic markers included in the gene-based tests. By default, FALSE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used 
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(1.5,2.5,3.5,4.5,5.5,10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param dosageZerodCutoff numeric. In gene- or region-based tests, for each variants with MAC <= 10, dosages <= dosageZerodCutoff with be set to 0. By default, 0.2. 
#' @param IsOutputPvalueNAinGroupTestforBinary logical. In gene- or region-based tests for binary traits. if IsOutputPvalueNAinGroupTestforBinary is TRUE, p-values without accounting for case-control imbalance will be output. By default, FALSE 
#' @param IsAccountforCasecontrolImbalanceinGroupTest logical. In gene- or region-based tests for binary traits. If IsAccountforCasecontrolImbalanceinGroupTest is TRUE, p-values after accounting for case-control imbalance will be output. By default, TRUE
#' @return SAIGEOutputFile
#' @export
SPAGMMATtest = function(bgenFile = "",
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
		 IsDropMissingDosages = FALSE,	
		 minMAC = 0.5, 
                 minMAF = 0,
		 maxMAFforGroupTest = 0.5,
        	 minInfo = 0,
                 GMMATmodelFile = "", 
                 varianceRatioFile = "", 
                 SPAcutoff=2, 
                 SAIGEOutputFile = "",
		 numLinesOutput = 10000, 
		 IsSparse=TRUE,
		 IsOutputAFinCaseCtrl=FALSE,
		 IsOutputNinCaseCtrl=FALSE,
		 LOCO=FALSE,
		 condition="",
		 sparseSigmaFile="",
		 groupFile="",
		 kernel="linear.weighted",
		 method="optimal.adj",
		 weights.beta.rare = c(1,25), 
		 weights.beta.common = c(0.5,0.5), 
		 weightMAFcutoff = 0.01,
		 r.corr=0,
		 IsSingleVarinGroupTest = TRUE,
		 cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5), 
		 cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5),
		 dosageZerodCutoff = 0.2,	
		 IsOutputPvalueNAinGroupTestforBinary = FALSE,
		 IsAccountforCasecontrolImbalanceinGroupTest = TRUE){


  if(weightMAFcutoff < 0 | weightMAFcutoff > 0.5){
    stop("weightMAFcutoff needs to be between 0 and 0.5\n")
  }


  adjustCCratioinGroupTest=TRUE
  if(!IsAccountforCasecontrolImbalanceinGroupTest){
    IsOutputPvalueNAinGroupTestforBinary = TRUE
    adjustCCratioinGroupTest = FALSE
  }

  # if group file is specified, the region-based test will be performed, otherwise, the single-variant assoc test will be performed. 



  if(groupFile == ""){
    isGroupTest = FALSE
    cat("single-variant association test will be performed\n")
  }else{
    cat("group-based test will be performed\n")
    
  if(dosageZerodCutoff < 0){
    dosageZerodCutoff = 0
  }else if(dosageZerodCutoff >= 0 ){
    cat("Any dosages <= ", dosageZerodCutoff, " for genetic variants with MAC <= 10 are set to be 0 in group tests\n")
  }
    if(!file.exists(groupFile)){
      stop("ERROR! groupFile ", groupFile, " does not exsit\n")
    }else{
      isGroupTest = TRUE
    }
  }

  if(file.exists(SAIGEOutputFile)){
    file.remove(SAIGEOutputFile)
  }

  if(!file.exists(SAIGEOutputFile)){
    file.create(SAIGEOutputFile, showWarnings = TRUE)
  }

  #file for the glmm null model
  if(!file.exists(GMMATmodelFile)){
    stop("ERROR! GMMATmodelFile ", GMMATmodelFile, " does not exsit\n")
  }else{
    load(GMMATmodelFile)
    obj.glmm.null = modglmm
    rm(modglmm)

    sampleInModel = NULL
    sampleInModel$IID = obj.glmm.null$sampleID
    sampleInModel = data.frame(sampleInModel)
    sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
    cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
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
    ln = length(cateVarRatioMinMACVecExclude)
    hn = length(cateVarRatioMaxMACVecInclude)
    if(nrow(varRatioData) == 1){
      #ratioVec = rep(varRatioData[1,1],6)
      ratioVec = varRatioData[1,1]
      cat("Single variance ratio is provided, so categorical variance ratio won't be used!\n")

      if(isGroupTest){
	stop("ERROR! To perform gene-based tests, categorical variance ratios are required\n")
      }	
    }else{
      ratioVec = varRatioData[,1]
      nrv = length(ratioVec)
      if (nrv !=  ln){	
	stop("ERROR! The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude\n")
      }
      if (ln != (hn+1)){	
	stop("ERROR! The length of cateVarRatioMaxMACVecInclude does not match with the lenght of cateVarRatioMinMACVecExclude (-1)\n")
      }		
    }
    #cat("variance Ratio is ", varRatio, "\n")
    cat("variance Ratio is ", ratioVec, "\n")
  }

  #sample file
  if(!file.exists(sampleFile)){
    stop("ERROR! sampleFile ", sampleFile, " does not exsit\n")
  }else{
    sampleListinDosage = data.frame(data.table:::fread(sampleFile, header=F, stringsAsFactors=FALSE, colClasses=c("character")))
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


      #rm(sampleListinDosage)
      rm(dataMerge)
      rm(dataMerge_v2)
      rm(dataMerge_sort)
      rm(dataMerge_v2_sort)
      rm(sampleInModel)


    }
  }


  ####check and read files
  #sparseSigmaFile
  if(sparseSigmaFile == ""){
    sparseSigma = NULL
    cat("sparse kinship matrix is not used\n")  
  }else{
    cat("sparse kinship matrix is going to be used\n")
    if(!file.exists(sparseSigmaFile)){
      stop("ERROR! sparseSigmaFile ", sparseSigmaFile, " does not exsit\n")
    }else{
      sparseSigma = Matrix:::readMM(sparseSigmaFile)
      cat("sparseSigmaFile: ", sparseSigmaFile, "\n")
    }
  }

  ##Needs to check the number of columns and the number of samples in sample file
  if(bgenFile != ""){ 
    if(!file.exists(bgenFile)){
      stop("ERROR! bgenFile ", bgenFile, " does not exsit\n")
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


  if(IsDropMissingDosages){
     cat("Samples with missing dosages will be dropped from the analysis\n")
  }else{
    cat("Missing dosages will be mean imputed for the analysis\n")
  }


  ##############START TEST########################
  startTime = as.numeric(Sys.time())  # start time of the SPAGMMAT tests
  cat("Analysis started at ", startTime, "Seconds\n")

#  if(file.exists(SAIGEOutputFile)){file.remove(SAIGEOutputFile)}
#  gc(verbose=T, full=T)

  if(!isGroupTest){
    if(dosageFileType == "bgen"){
      dosageFilecolnamesSkip = c("CHR","POS","rsid","SNPID","Allele1","Allele2", "AC_Allele2", "AF_Allele2", "imputationInfo")
    }else if(dosageFileType == "vcf"){
      dosageFilecolnamesSkip = c("CHR","POS","SNPID","Allele1","Allele2", "AC_Allele2", "AF_Allele2", "imputationInfo")
    }
  }

  if(condition != ""){
    isCondition = TRUE
  }else{
    isCondition = FALSE
  }

  cat("isCondition is ", isCondition, "\n")

  if(isCondition){
    condition_original=unlist(strsplit(condition,","))

    if(length(condition_original) > 1){
    	condition_new=NULL
    	for(x in 1:length(condition_original)){
	  condition_new = rbind(condition_new, c(as.numeric(strsplit(strsplit(condition_original[x], ":")[[1]][2][1], "_")[[1]][1]), condition_original[x]))
    	}
    	condition_new2 = condition_new[order(as.numeric(condition_new[,1])),] 
    
    	conditionlist = paste(c("condMarkers",condition_new2[,2]),collapse="\t")
    }else{
    	conditionlist= paste(c("condMarkers",unlist(strsplit(condition,","))),collapse="\t")    
    }
#    conditionlist = paste(c("condMarkers",unlist(strsplit(condition,","))),collapse="\t")
    cat("conditionlist is ", conditionlist, "\n")

    if(dosageFileType == "vcf"){
      setMAFcutoffs(0, 0.5)
      isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
      SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N)
      Gx_cond = getGenoOfGene_vcf(conditionlist, 0)
      if(Gx_cond$cnt > 0){
        dosage_cond = Matrix:::sparseMatrix(i = as.vector(Gx_cond$iIndex), j = as.vector(Gx_cond$jIndex), x = as.vector(Gx_cond$dosages), symmetric = FALSE, dims = c(N, Gx_cond$cnt))
      }

    }else if(dosageFileType == "bgen"){
      SetSampleIdx(sampleIndex, N)
      Gx_cond = getGenoOfGene_bgen(bgenFile,bgenFileIndex, conditionlist)
      if(Gx_cond$cnt > 0){
        dosage_cond = matrix(Gx_cond$dosages, byrow=F, ncol = Gx_cond$cnt)	
        dosage_cond = as(dosage_cond, "sparseMatrix") 
      }	
    }else{
      stop("ERROR: conditional analysis can only work for dosageFileType vcf, sav or bgen\n")
    }

    cat("conditioning on ", unlist(Gx_cond$markerIDs), "\n")
    cntMarker = Gx_cond$cnt

    cat("isCondition is ", isCondition, "\n")

    if(cntMarker == 0){
        
      cat("WARNING: conditioning markers are not found in the provided dosage file \n")
      isCondition = FALSE
      dosage_cond = NULL
    }

  }else{#end of if(isCondition){
    dosage_cond = NULL
  }



  ########Binary traits####################
  if(traitType == "binary"){
    cat("It is a binary trait\n")
    if(!isGroupTest){
      if(!isCondition){
	resultHeader = c(dosageFilecolnamesSkip, "N", "BETA", "SE", "Tstat", "p.value", "p.value.NA", "Is.SPA.converge","varT","varTstar")
      }else{
	resultHeader = c(dosageFilecolnamesSkip, "N", "BETA", "SE", "Tstat", "p.value", "p.value.NA", "Is.SPA.converge","varT","varTstar", "Tstat_cond", "p.value_cond", "varT_cond", "BETA_cond", "SE_cond")
      }
 
      if(IsOutputAFinCaseCtrl){
        resultHeader = c(resultHeader, "AF.Cases", "AF.Controls")
      }
      if(IsOutputNinCaseCtrl){
	resultHeader = c(resultHeader, "N.Cases", "N.Controls")
      } 	
      write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
    } #if(!isGroupTest){ 

    if(SPAcutoff < 10^-2){
      Cutoff=10^-2
    }else{
      Cutoff = SPAcutoff
    }

    y = obj.glmm.null$obj.glm.null$y
    y1Index = which(y == 1)
    NCase = length(y1Index)
    y0Index = which(y == 0)
    NCtrl = length(y0Index)

    cat("Analyzing ", NCase, " cases and ",NCtrl, " controls \n")
    N = length(y)
    obj.glmm.null$obj.noK$XVX_inv_XV = obj.glmm.null$obj.noK$XXVX_inv * obj.glmm.null$obj.noK$V
    indChromCheck = FALSE
    if(!obj.glmm.null$LOCO){
      mu = obj.glmm.null$fitted.values
      mu.a<-as.vector(mu)
      mu2.a<-mu.a *(1-mu.a)
      obj.glmm.null$obj.noK$XVX = t(obj.glmm.null$obj.noK$X1) %*% (obj.glmm.null$obj.noK$X1 * mu2.a)
      obj.glmm.null$obj.noK$S_a = colSums(obj.glmm.null$obj.noK$X1 * (y - mu.a))

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
      obj.glmm.null$obj.noK$XVX = t(obj.glmm.null$obj.noK$X1) %*% (obj.glmm.null$obj.noK$X1 * mu2.a)
      obj.glmm.null$obj.noK$S_a = colSums(obj.glmm.null$obj.noK$X1 * (y - mu.a))
    }else{
      cat("WARNING: LOCO will be used, but chromosome for the dosage file is not specified. Will check each marker for its chromosome for LOCO!\n")
      indChromCheck = TRUE
    }
#####Quantitative traits##########

  }else if(traitType == "quantitative"){
    cat("It is a quantitative trait\n")
    adjustCCratioinGroupTest = FALSE	
    if(!isGroupTest){
      if(!isCondition){
        resultHeader = c(dosageFilecolnamesSkip,  "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar")
      }else{
        resultHeader = c(dosageFilecolnamesSkip,  "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar","Tstat_cond", "p.value_cond", "varT_cond", "BETA_cond", "SE_cond" )
      }
      write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
    }

    y = obj.glmm.null$obj.glm.null$y
    N = length(y)
    mu2.a = rep(1, N)
    tauVec = obj.glmm.null$theta
    obj.glmm.null$obj.noK$XVX = t(obj.glmm.null$obj.noK$X1) %*% (obj.glmm.null$obj.noK$X1)
    obj.glmm.null$obj.noK$XVX_inv_XV = obj.glmm.null$obj.noK$XXVX_inv * obj.glmm.null$obj.noK$V
    indChromCheck = FALSE


    #cat("obj.glmm.null$LOCO ", obj.glmm.null$LOCO, "\n")
    if(!obj.glmm.null$LOCO){
      mu = obj.glmm.null$fitted.values
      mu.a<-as.vector(mu)
      obj.glmm.null$obj.noK$S_a = colSums(obj.glmm.null$obj.noK$X1 * (y - mu.a))

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
      obj.glmm.null$obj.noK$S_a = colSums(obj.glmm.null$obj.noK$X1 * (y - mu.a))

    }else{
      cat("WARNING: LOCO will be used, but chromosome for the dosage file is not specified. Will check each marker for its chromosome for LOCO!\n")
      indChromCheck = TRUE
    }

  }else{
    stop("ERROR! The type of the trait has to be either binary or quantitative\n")
  }



  if(nrow(varRatioData) == 1){
    cateVarRatioMinMACVecExclude=c(0)
    cateVarRatioMaxMACVecInclude=c(2*N)
  }


  if(isCondition){
    condpre = getCovMandOUT_cond_pre(dosage_cond=dosage_cond, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec, obj.glmm.null = obj.glmm.null, sparseSigma = sparseSigma, IsSparse=IsSparse, mu = mu, mu.a = mu.a, mu2.a = mu2.a, Cutoff = Cutoff)
    OUT_cond = condpre$OUT_cond
    G2tilde_P_G2tilde_inv = condpre$G2tilde_P_G2tilde_inv	
  }else{# end of if(isCondition)
    OUT_cond = NULL
    G2tilde_P_G2tilde_inv = NULL
  }

  cat("isCondition is ", isCondition, "\n")

#  gc(verbose=T, full=T)
#determine minimum MAF for markers to be tested
  if(minMAC == 0){
    minMAC = 0.5
    cat("As minMAC is set to be 0, minMAC = 0.5 will be used\n")
  } ##01-19-2018
  cat("minMAC: ",minMAC,"\n")
  cat("minMAF: ",minMAF,"\n")
  minMAFBasedOnMAC = minMAC/(2*N) 
  testMinMAF = max(minMAFBasedOnMAC, minMAF) 
  cat("Minimum MAF of markers to be tested is ", testMinMAF, "\n")

  ##############START TEST########################
  startTime = as.numeric(Sys.time())  # start time of the SPAGMMAT tests
  cat("Analysis started at ", startTime, "Seconds\n")

  if(!isGroupTest){

    isVariant = TRUE
    if (dosageFileType == "bgen"){
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
        rangesExclude = data.table:::fread(rangestoExcludeFile, header=F, colClasses = c("character", "numeric", "numeric"))
        ranges_to_exclude = data.frame(rangesExclude)
        colnames(ranges_to_exclude) = c("chromosome","start","end")  
      }else{
        ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
      }

      if(rangestoIncludeFile != ""){
        rangesInclude = data.table:::fread(rangestoIncludeFile, header=F, colClasses = c("character", "numeric", "numeric"))
        ranges_to_include = data.frame(rangesInclude)
        colnames(ranges_to_include) = c("chromosome","start","end")
      }else{
        ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
      }

      Mtest = setgenoTest_bgenDosage(bgenFile,bgenFileIndex, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
      if(Mtest == 0){
        isVariant = FALSE
        stop("ERROR! Failed to open ", bgenFile, "\n")
      }
      isQuery = getQueryStatus()
      SetSampleIdx(sampleIndex, N)

	
      nsamplesinBgen = getSampleSizeinBgen()
      if(nrow(sampleListinDosage) != nsamplesinBgen){
	stop("ERROR! The number of samples specified in the sample file does not equal to the number of samples in the bgen file\n")
      }		


    }else if(dosageFileType == "vcf"){

      setgenoTest_vcfDosage(vcfFile,vcfFileIndex,vcfField,ids_to_exclude_vcf = idstoExcludeFile, ids_to_include_vcf = idstoIncludeFile, chrom, start, end)
    #setTestField(vcfField)
      isVariant = getGenoOfnthVar_vcfDosage_pre()
      SetSampleIdx_vcfDosage(sampleIndex, N)
      nsamplesinVCF = getSampleSizeinVCF()
     if(nrow(sampleListinDosage) != nsamplesinVCF){
        stop("ERROR! The number of samples specified in the sample file does not equal to the number of samples in the VCF file\nPlease check again. Please note that the sample file needs to have no header.")
    }


    }

    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
    OUT = NULL
    numPassMarker = 0
    mth = 0

    while(isVariant){
      mth = mth + 1
      if (dosageFileType == "bgen"){
        if(isQuery){
          Gx = getDosage_bgen_withquery()
        }else{
          Gx = getDosage_bgen_noquery()
        }
        markerInfo = getMarkerInfo()
        if(!(markerInfo >= 0 & markerInfo <= 1)){markerInfo=1}
        G0 = Gx$dosages
        AC = Gx$variants$AC
        AF = Gx$variants$AF
        Gx$variants$markerInfo = markerInfo
        rowHeader=as.vector(unlist(Gx$variants))
	#cat("rowHeader: ", rowHeader, "\n")
        if(indChromCheck){
	  CHR = Gx$variants$chromosome
	  cat("CHR ", CHR , "\n")	
        }	

        if(Mtest == mth){isVariant = FALSE}
        indexforMissing = Gx$indexforMissing

      }else if(dosageFileType == "vcf"){
        Gx = getGenoOfnthVar_vcfDosage(mth)
        G0 = Gx$dosages
        AC = Gx$variants$AC
        AF = Gx$variants$AF
        markerInfo = Gx$variants$markerInfo
        if(!(markerInfo >= 0 & markerInfo <= 1)){markerInfo=1; Gx$variants$markerInfo=1}
        rowHeader=as.vector(unlist(Gx$variants))
        if(indChromCheck){
          CHR = Gx$variants$chromosome
	  cat("CHR ", CHR , "\n")
        }
        isVariant = getGenoOfnthVar_vcfDosage_pre()
        indexforMissing = Gx$indexforMissing
      }

      MAC = AC
      MAF = AF
      if(AF > 0.5){
        MAC = 2*N - AC
        MAF = 1 - AF
      }
     #if(dosageZerodCutoff > 0){
     # 	if(MAC <= 10){
     #		G0[which(G0 <= dosageZerodCutoff)] = 0
     #		MAC = sum(G0)
     #		MAF = MAC/(2*length(G0))
     # 		cat("Any dosages <= ", dosageZerodCutoff, " are set to be 0")
     #	}
     #}  

      if(MAF >= testMinMAF & markerInfo >= minInfo){
         numPassMarker = numPassMarker + 1
         varRatio = getVarRatio(G0, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)

         if(indChromCheck){
           CHR = as.character(CHR)
           CHRv2 = as.numeric(gsub("[^0-9.]", "", CHR))
           cat("CHR is ", CHR, "\n") 
    
           if(obj.glmm.null$LOCOResult[[CHRv2]]$isLOCO){
             mu = obj.glmm.null$LOCOResult[[CHRv2]]$fitted.values
             mu.a<-as.vector(mu)
	     if(traitType == "binary"){
               mu2.a<-mu.a *(1-mu.a)
	     }else if(traitType == "quantitative"){
		mu2.a = rep(1,N)
	     }	
      	   }else{
             mu = obj.glmm.null$fitted.values
             mu.a<-as.vector(mu)
	     if(traitType == "binary"){
               mu2.a<-mu.a *(1-mu.a)
	     }else if(traitType == "quantitative"){
                mu2.a = rep(1,N)
             }
           }

      	   obj.glmm.null$obj.noK$XVX = t(obj.glmm.null$obj.noK$X1) %*% (obj.glmm.null$obj.noK$X1 * mu2.a)
           obj.glmm.null$obj.noK$XVX_inv_XV = obj.glmm.null$obj.noK$XXVX_inv * obj.glmm.null$obj.noK$V
           obj.glmm.null$obj.noK$S_a = colSums(obj.glmm.null$obj.noK$X1 * (y - mu.a))
         }


	if(IsDropMissingDosages & isCondition){
		indexforMissing = unique(c(indexforMissing, Gx_cond$indexforMissing))
	}


    if(IsDropMissingDosages & length(indexforMissing) > 0){
        missingind = seq(1, length(G0))[-(indexforMissing + 1)]
	cat("Removing ", length(indexforMissing), " samples with missing dosages/genotypes\n")

        G0 = G0[missingind]
	subsetModelResult = subsetModelFileforMissing(obj.glmm.null, missingind, mu, mu.a ,mu2.a)	
	obj.glmm.null.sub = subsetModelResult$obj.glmm.null.sub
	mu.a.sub = subsetModelResult$mu.a.sub
	mu.sub = subsetModelResult$mu.sub
	mu2.a.sub = subsetModelResult$mu2.a.sub
	rm(subsetModelResult)

        y.sub = obj.glmm.null.sub$obj.glm.null$y
	N.sub = length(G0)
	AC_Allele2.sub = sum(G0)
	AF_Allele2.sub = AC_Allele2.sub/(2*N.sub)
	rowHeader[6] = AC_Allele2.sub
	rowHeader[7] = AF_Allele2.sub

        if(traitType == "binary"){
          y1Index.sub = which(y.sub == 1)
          NCase.sub = length(y1Index.sub)
          y0Index.sub = which(y.sub == 0)
	  NCtrl.sub = length(y0Index.sub)
	}

	sparseSigma.sub = sparseSigma
	if(!is.null(sparseSigma)){sparseSigma.sub = sparseSigma[missingind, missingind]}
	####Update the conditional analysis after dropping missing genotypes
        if(isCondition){
        	cat("Removing ", length(indexforMissing), " samples from the conditional marker\n")
                dosage_cond.sub = dosage_cond[missingind, ]
                dosage_cond.sub = as(dosage_cond.sub, "sparseMatrix")
                ######re-test the conditional variants after removing samples with missing genotypes
		condpre.sub = getCovMandOUT_cond_pre(dosage_cond=dosage_cond.sub, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec, obj.glmm.null = obj.glmm.null.sub, sparseSigma = sparseSigma.sub, IsSparse=IsSparse,  mu = mu.sub, mu.a = mu.a.sub, mu2.a = mu2.a.sub, Cutoff = Cutoff)
    		OUT_cond.sub = condpre$OUT_cond
    		G2tilde_P_G2tilde_inv.sub = condpre.sub$G2tilde_P_G2tilde_inv
		condpre2.sub = getCovMandOUT_cond(G0 = G0, dosage_cond = dosage_cond.sub, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, obj.glmm.null = obj.glmm.null.sub, sparseSigma = sparseSigma.sub, covM = condpre.sub$covM, mu2.a = mu2.a.sub)
          	G1tilde_P_G2tilde.sub = condpre2.sub$G1tilde_P_G2tilde
          	GratioMatrixall.sub = condpre2.sub$GratioMatrixall
        }

	if(traitType == "binary"){
	  if (NCase.sub == 0 | NCtrl.sub == 0) {	
	   #out1 = c(rep(NA, 8), NCase.sub, NCtrl.sub)
	   out1 = c(rep(NA, 8))
	   OUTvec=c(rowHeader, N.sub, unlist(out1))
            # OUT = rbind(OUT, c(rowHeader, N.sub, unlist(out1)))	
	   if(IsOutputAFinCaseCtrl){
	     if(NCase.sub == 0){
		AFCase = NA
		AFCtrl = sum(G0[y0Index.sub])/(2*NCtrl.sub)		
	     }else if(NCtrl.sub == 0){
		AFCtrl = NA
		AFCase = sum(G0[y1Index.sub])/(2*NCase.sub)
	     }	
	     #OUT = rbind(OUT, c(rowHeader, N.sub, unlist(out1), AFCase, AFCtrl))	
	     OUTvec=c(OUTvec, AFCase, AFCtrl)
	   }

	   if(IsOutputNinCaseCtrl){
	     OUTvec=c(OUTvec, NCase.sub, NCtrl.sub)			
	   }
	   OUT = rbind(OUT, OUTvec)
	   OUTvec=NULL
	  }else{ #if (NCase.sub == 0 | NCtrl.sub == 0) {
           out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0, AC, AF, MAF, IsSparse, obj.glmm.null.sub$obj.noK, mu.a.sub, mu2.a.sub, y.sub, varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma.sub, isCondition=isCondition, OUT_cond=OUT_cond.sub, G1tilde_P_G2tilde = G1tilde_P_G2tilde.sub, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv.sub)
	  OUTvec=c(rowHeader, N.sub, unlist(out1))

	   #if(!IsOutputAFinCaseCtrl){
           #  OUT = rbind(OUT, c(rowHeader, N.sub, unlist(out1)))
           #}else{
	   if(IsOutputAFinCaseCtrl){
             AFCase = sum(G0[y1Index.sub])/(2*NCase.sub)
             AFCtrl = sum(G0[y0Index.sub])/(2*NCtrl.sub)
	     OUTvec=c(OUTvec, AFCase, AFCtrl)
             #OUT = rbind(OUT, c(rowHeader, N.sub, unlist(out1), AFCase, AFCtrl))
            }
	   if(IsOutputNinCaseCtrl){
             OUTvec=c(OUTvec, NCase.sub, NCtrl.sub)
           }
	   OUT = rbind(OUT, OUTvec)
		OUTvec=NULL
	  }
         }else if(traitType == "quantitative"){

           out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0,obj.glmm.null.sub$obj.noK, AC, AF, y.sub, mu.sub, varRatio, tauVec, sparseSigma=sparseSigma.sub, isCondition=isCondition, OUT_cond=OUT_cond.sub, G1tilde_P_G2tilde = G1tilde_P_G2tilde.sub, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv.sub)

           if(!isCondition){
             OUT = rbind(OUT, c(rowHeader, N.sub, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2))
           }else{
             OUT = rbind(OUT, c(rowHeader, N.sub, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2, out1$Tstat_c,  out1$p.value.c, out1$var1_c, out1$BETA_c, out1$SE_c))
           }
         }
	
     }else{ #if(IsDropMissingDosages & length(indexforMissing) > 0){
	          ##conditional analysis
         if(isCondition){
           condpre2 = getCovMandOUT_cond(G0 = G0, dosage_cond = dosage_cond, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, obj.glmm.null = obj.glmm.null, sparseSigma = sparseSigma, covM = condpre$covM, mu2.a = mu2.a)
           G1tilde_P_G2tilde = condpre2$G1tilde_P_G2tilde
           GratioMatrixall = condpre2$GratioMatrixall

         }else{ #end of if(isCondition)
           G1tilde_P_G2tilde = NULL
           GratioMatrixall = NULL
         }	





  	 if(traitType == "binary"){
           out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0, AC, AF, MAF, IsSparse, obj.glmm.null$obj.noK, mu.a, mu2.a, y, varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv)
	   OUTvec=c(rowHeader, N,unlist(out1))

    	   if(IsOutputAFinCaseCtrl){	     	
      	     AFCase = sum(G0[y1Index])/(2*NCase)
      	     AFCtrl = sum(G0[y0Index])/(2*NCtrl)
		OUTvec=c(OUTvec, AFCase, AFCtrl)

           }
	   if(IsOutputNinCaseCtrl){
	     OUTvec=c(OUTvec, NCase, NCtrl)
           }
	   OUT = rbind(OUT, OUTvec)
	   OUTvec=NULL

         }else if(traitType == "quantitative"){

           out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.glmm.null$obj.noK, AC, AF, y, mu, varRatio, tauVec, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv)

           if(!isCondition){
             OUT = rbind(OUT, c(rowHeader, N, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2))
           }else{
             OUT = rbind(OUT, c(rowHeader, N, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2, out1$Tstat_c,  out1$p.value.c, out1$var1_c, out1$BETA_c, out1$SE_c))
           }
         }
  
     } #end of else{ #if(IsDropMissingDosages & length(indexforMissing) > 0){


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
     } ####end of while(isVariant)

   }else{ #end if(!isGroupTest){
   #########Group Test

   
     OUT_single = NULL
     if(IsSingleVarinGroupTest){
       SAIGEOutputFile_single = paste0(SAIGEOutputFile, "_single")
     
       headerline = c("markerID", "AC", "AF", "N", "BETA", "SE", "Tstat", "p.value","varT","varTstar") 
	 if(traitType=="binary"){
	   headerline = c(headerline, "AF.Cases", "AF.Controls", "N.Cases", "N.Controls")	
	 }	
       write(headerline,file = SAIGEOutputFile_single, ncolumns = length(headerline))
     }
    	 
     if(dosageFileType == "bgen"){
       SetSampleIdx(sampleIndex, N)
     }else if(dosageFileType == "vcf"){
       setMAFcutoffs(testMinMAF, maxMAFforGroupTest)
       cat("genetic variants with ", testMinMAF, "<= MAF <= ", maxMAFforGroupTest, "are included for gene-based tests\n")
       isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
       SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N)
     }


     OUT = NULL
     if(traitType == "quantitative"){
       cat("It is a quantitative trait\n")
	IsOutputPvalueNAinGroupTestforBinary=TRUE
	adjustCCratioinGroupTest = FALSE
     }else if(traitType == "binary"){
       cat("It is a binary trait\n")
       #cat("WARNING!!!! Gene-based tests do not work for binary traits with unbalanced case-control ratios (disease prevalence < 20%)! \n")	
       #adjustCCratioinGroupTest = TRUE
       #if(isCondition){	
#	adjustCCratioinGroupTest = FALSE
#       	cat("WARNING!!!! Case-control imbalance is not adjusted for binary traits to perform conditional analysis. Do not specify condition= if needs to account for case-control imbalance\n")	
 #      }else 
	if(adjustCCratioinGroupTest){
          cat("Case-control imbalance is adjusted for binary traits.\n")		
	}
	if(IsOutputPvalueNAinGroupTestforBinary){
	  cat("P-values without case-control imbalance will be output.\n")
	}

	obj.glmm.null$obj_cc = SKAT::SKAT_Null_Model(obj.glmm.null$obj.glm.null$y ~ obj.glmm.null$obj.noK$X1-1, out_type="D", Adjustment = FALSE)
     }


       mth = 0
       MACcateNumHeader = paste0("Nmarker_MACCate_", seq(1,length(cateVarRatioMinMACVecExclude)))
       if(!isCondition){
	  if(adjustCCratioinGroupTest){	
           resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs")
	   if(method=="optimal.adj"){
	     resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_SKAT")	
	   }
	  }

	 if(IsOutputPvalueNAinGroupTestforBinary){
           if(!adjustCCratioinGroupTest){
             resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs")
             if(method=="optimal.adj"){
               resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_SKAT")
             }
           }else{
             resultHeader = c(resultHeader, "Pvalue.NA")
	     if(method=="optimal.adj"){
               resultHeader = c(resultHeader, "Pvalue_Burden.NA","Pvalue_SKAT.NA")
             }	
	   }
	 }
       }else{
	 if(adjustCCratioinGroupTest){
           resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs")
           if(method=="optimal.adj"){
             resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_Burden_cond","Pvalue_SKAT","Pvalue_SKAT_cond")
           }
          }

	if(IsOutputPvalueNAinGroupTestforBinary){
           if(!adjustCCratioinGroupTest){
	     resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs")	
		if(method=="optimal.adj"){

	     resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs", "Pvalue_Burden","Pvalue_Burden_cond","Pvalue_SKAT","Pvalue_SKAT_cond")	
             }	
	   }else{
			resultHeader = c(resultHeader,"Pvalue.NA", "Pvalue.NA_cond")
		if(method=="optimal.adj"){
			resultHeader = c(resultHeader,"Pvalue_Burden.NA","Pvalue_Burden.NA_cond","Pvalue_SKAT.NA","Pvalue_SKAT.NA_cond")
		}
	  }	
        } 

      }

	#if(adjustCCratioinGroupTest){
	#	if(IsOutputPvalueNAinGroupTestforBinary){
	#		if(method=="optimal.adj"){
	#			resultHeader = c(resultHeader, "Pvalue_skato_NA", "Pvalue_burden_NA", "Pvalue_skat_NA")
	#		}else{
	#			resultHeader = c(resultHeader, "Pvalue_NA")
	#		}
	#	}
		#resultHeader = c(resultHeader, "Pvalue_skato_old", "Pvalue_burden_old", "Pvalue_skat_old", "Pvalue_skato_new", "Pvalue_burden_new", "Pvalue_skat_new", "Pvalue_skato_new2", "Pvalue_burden_new2", "Pvalue_skat_new2")	
	#}


       write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
#	gc(verbose=T, full=T)	

       cat("isCondition is ", isCondition, "\n")


       gf = file(groupFile, "r")
       while ( TRUE ) {
         marker_group_line = readLines(gf, n = 1)
         if(length(marker_group_line) == 0 ){
	   break	
         }else{
	   marker_group_line_list = strsplit(marker_group_line, split="\t")[[1]] 
	   geneID = marker_group_line_list[1]
           #geneID = strsplit(marker_group_line, split="\t")[[1]][1]
	   cat("geneID: ", geneID, "\n")	
	   if(length(marker_group_line_list) <= 1){
	    stop("no marker IDs are found for gene ", geneID, ". Please make sure the group file is tab delimited.", "\n")
	   }

           if(dosageFileType == "vcf"){
             Gx = getGenoOfGene_vcf(marker_group_line, minInfo)

           }else if(dosageFileType == "bgen"){
	     print(marker_group_line)
	     cat("genetic variants with ", testMinMAF, "<= MAF <= ", maxMAFforGroupTest, "are included for gene-based tests\n") 
             Gx = getGenoOfGene_bgen(bgenFile,bgenFileIndex, marker_group_line, testMinMAF, maxMAFforGroupTest, minInfo)
           }
           cntMarker = Gx$cnt
           cat("cntMarker: ", cntMarker, "\n")
           if(cntMarker > 0){
             #Gmat = matrix(Gx$dosages, byrow=F, ncol = cntMarker)
	     #Gx$dosages = NULL	

		if(dosageFileType == "vcf"){
			Gmat = Matrix:::sparseMatrix(i = as.vector(Gx$iIndex), j = as.vector(Gx$jIndex), x = as.vector(Gx$dosages), symmetric = FALSE, dims = c(N, cntMarker))
		}else{
			Gmat = matrix(Gx$dosages, byrow=F, ncol = cntMarker)
			Gmat = as(Gmat, "sparseMatrix")	
		}



		if(isCondition){
			indexforMissing = unique(c(Gx$indexforMissing, Gx_cond$indexforMissing))
		}else{
			indexforMissing = unique(Gx$indexforMissing)
		}
		



	     if(IsDropMissingDosages & length(indexforMissing) > 0){	

	        missingind = seq(1, nrow(Gmat))[-(indexforMissing + 1)]
		subsetModelResult = subsetModelFileforMissing(obj.glmm.null, missingind, mu, mu.a ,mu2.a)
        	obj.glmm.null.sub = subsetModelResult$obj.glmm.null.sub
        	mu.a.sub = subsetModelResult$mu.a.sub
        	mu.sub = subsetModelResult$mu.sub
       		mu2.a.sub = subsetModelResult$mu2.a.sub
        	rm(subsetModelResult)

		y.sub = obj.glmm.null.sub$obj.glm.null$y

		obj.glmm.null.sub$obj_cc = SKAT::SKAT_Null_Model(y.sub ~ obj.glmm.null.sub$obj.noK$X1-1, out_type="D", Adjustment = FALSE)

        	if(traitType == "binary"){
          		y1Index.sub = which(y.sub == 1)
          		NCase.sub = length(y1Index.sub)
          		y0Index.sub = which(y.sub == 0)
          		NCtrl.sub = length(y0Index.sub)
        	}

        	sparseSigma.sub = sparseSigma
        	if(!is.null(sparseSigma)){sparseSigma.sub = sparseSigma[missingind, missingind]}
	       	Gmat = Gmat[missingind,]

		if(isCondition){
                	cat("Removing ", length(indexforMissing), " samples from the conditional marker\n")
                	dosage_cond.sub = dosage_cond[missingind, ]
                	dosage_cond.sub = as(dosage_cond.sub, "sparseMatrix")


			condpre.sub = getCovMandOUT_cond_pre(dosage_cond=dosage_cond.sub, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec, obj.glmm.null = obj.glmm.null.sub, sparseSigma = sparseSigma.sub, IsSparse=IsSparse,  mu = mu.sub, mu.a = mu.a.sub, mu2.a = mu2.a.sub, Cutoff = Cutoff)
	                OUT_cond.sub = condpre.sub$OUT_cond
        	        G2tilde_P_G2tilde_inv.sub = condpre.sub$G2tilde_P_G2tilde_inv

		}else{ #if(isCondition){
			dosage_cond.sub = NULL
			OUT_cond.sub = NULL
		}

              } #if(IsDropMissingDosages & length(indexforMissing) > 0){ 


	      rmMarkerIndex = NULL
	      if(dosageZerodCutoff > 0 & sum(Gx$MACs <= 10) > 0){
		zerodIndex = which(Gx$MACs <= 10)
		for(z in zerodIndex){
			cat("z is ", z, "\n")
			cat("dim(Gmat) is ", dim(Gmat), "\n")
			if(dim(Gmat)[2] > 1){
				replaceindex = which(Gmat[,z] <= dosageZerodCutoff)
				if(length(replaceindex) > 0){	
					Gmat[replaceindex,z] = 0 
					cat("dim(Gmat) is ", dim(Gmat), "\n")
					if(sum(Gmat[,z])/(2*nrow(Gmat)) < testMinMAF){rmMarkerIndex = c(rmMarkerIndex, z)}
				}
			}else{
				Gmat[which(Gmat <= dosageZerodCutoff)] = 0
				if(sum(Gmat)/(2*length(Gmat)) < testMinMAF){rmMarkerIndex = c(rmMarkerIndex, z)} 
				Gmat = matrix(Gmat, ncol=1)
				Gmat = as(Gmat, "sparseMatrix")
			}
		}
		cat("length(rmMarkerIndex): ", length(rmMarkerIndex), "\n")
		if(length(rmMarkerIndex) > 0){
			cntMarker = cntMarker - length(rmMarkerIndex)
			if(cntMarker > 0){
				Gmat = Gmat[,-rmMarkerIndex]
				cntMarker = cntMarker - length(rmMarkerIndex)
				Gx$markerIDs = Gx$markerIDs[-rmMarkerIndex]
				Gx$markerAFs = Gx$markerAFs[-rmMarkerIndex]
			}
		}
		Gmat = as(Gmat, "sparseMatrix")		
	     }
         
          }#if(cntMarker > 0){


	    if(cntMarker > 0){
		if(IsDropMissingDosages & length(indexforMissing) > 0){
		  cat("isCondition is ", isCondition, "\n")
		groupTestResult = groupTest(Gmat = Gmat, obj.glmm.null = obj.glmm.null.sub, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, G2_cond = dosage_cond.sub, G2_cond_es = OUT_cond.sub[,1], kernel = kernel, method = method, weights.beta.rare = weights.beta.rare, weights.beta.common = weights.beta.common, weightMAFcutoff = weightMAFcutoff, r.corr = r.corr, max_maf = maxMAFforGroupTest, sparseSigma = sparseSigma.sub, mu.a = mu.a.sub, mu2.a = mu2.a.sub, IsSingleVarinGroupTest = IsSingleVarinGroupTest, markerIDs = Gx$markerIDs, markerAFs = Gx$markerAFs, IsSparse= IsSparse, geneID = geneID, Cutoff = Cutoff, adjustCCratioinGroupTest = adjustCCratioinGroupTest, IsOutputPvalueNAinGroupTestforBinary = IsOutputPvalueNAinGroupTestforBinary)
	      }else{#if(IsDropMissingDosages & length(indexforMissing) > 0){	
		cat("isCondition is ", isCondition, "\n")
		#Gmat0 = as.matrix(Gmat)
		#write.table(Gmat0, "/net/hunt/disk2/zhowei/project/SAIGE_SKAT/typeIError_simuUsingRealData/quantitative/SAIGE/step2/C1orf122_seed256Phneo.geno.txt", col.names=F, row.names=F, quote=F)	
		#Gmat = round(Gmat)
		groupTestResult = groupTest(Gmat = Gmat, obj.glmm.null = obj.glmm.null, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, G2_cond = dosage_cond, G2_cond_es = OUT_cond[,1], kernel = kernel, method = method, weights.beta.rare = weights.beta.rare, weights.beta.common = weights.beta.common, weightMAFcutoff = weightMAFcutoff, r.corr = r.corr, max_maf = maxMAFforGroupTest, sparseSigma = sparseSigma, mu.a = mu.a, mu2.a = mu2.a, IsSingleVarinGroupTest = IsSingleVarinGroupTest, markerIDs = Gx$markerIDs, markerAFs = Gx$markerAFs, IsSparse= IsSparse, geneID = geneID, Cutoff = Cutoff, adjustCCratioinGroupTest = adjustCCratioinGroupTest, IsOutputPvalueNAinGroupTestforBinary = IsOutputPvalueNAinGroupTestforBinary)	
	     }
	    outVec = groupTestResult$outVec
	    OUT = rbind(OUT, outVec)
	    if(IsSingleVarinGroupTest){
            	outsingle = as.data.frame(groupTestResult$OUT_single)
            	OUT_single = rbind(OUT_single, outsingle)
            }


            mth = mth + 1
            if(mth %% numLinesOutput == 0){
              ptm <- proc.time()
              print(ptm)
              print(mth)
              OUT = as.data.frame(OUT)
              write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
              OUT = NULL
	      if(IsSingleVarinGroupTest){	
	        write.table(OUT_single, SAIGEOutputFile_single, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)	
	        OUT_single = NULL
              }	
            }
          }#if(cntMarker > 0){
      }#end of else for if(length(line) == 0 )
    } # end of while ( TRUE ) {

    if(!is.null(OUT)){
      OUT = as.data.frame(OUT)
      write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
      OUT = NULL
      if(IsSingleVarinGroupTest){
        #OUT_single = as.data.frame(OUT_single)
        write.table(OUT_single, SAIGEOutputFile_single, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
        OUT_single = NULL
      }	

    }
  #}else{
  #  stop("ERROR! The type of the trait has to be quantitative\n")
  #}  
 
 
}#if(groupTest)

  if (dosageFileType == "bgen"){
    closetestGenoFile_bgenDosage()
  }else if(dosageFileType == "vcf"){
    closetestGenoFile_vcfDosage()
  }
  summary(warnings())
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
  g1<-G[idx_no0]
  noCov = FALSE
  if(dim(obj.null$X1)[2] == 1){
    noCov = TRUE 
  }


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
    out1 = SPAtest:::Saddle_Prob_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="P")
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
  ##As g was not devided by sqrt(AC) 
  logOR = (Tstat/var1)/sqrt(AC)
  #logOR = Tstat/var1

  SE = abs(logOR/qnorm(out1$p.value/2))
  out1 = c(out1, BETA = logOR, SE = SE, Tstat = Tstat)


  ##condition
  Tstat1c = Tstat*sqrt(AC) - innerProduct(as.vector(OUT_cond[,1]), covariateVec)

  if(ncol(covM) > 2){
  G2_tildeG2_tilde1_v1  = t(covM[2:ncol(covM),2:ncol(covM)])
  diag(G2_tildeG2_tilde1_v1) = 0
  G2_tildeG2_tilde = G2_tildeG2_tilde1_v1 + covM[2:ncol(covM),2:ncol(covM)]


  covMsub = matrix(as.vector(covM[1,2:ncol(covM)]), nrow=1)
  var1c = varRatio*covM[1,1] - varRatio*covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub)
  }else{
    G2_tildeG2_tilde = covM[2:ncol(covM),2:ncol(covM)]
    covMsub = matrix(as.vector(covM[1,2:ncol(covM)]), nrow=1)
    var1c = varRatio*covM[1,1] - varRatio*(covMsub %*% solve(G2_tildeG2_tilde) %*% t(covMsub))
  }
  #var1c = var1 - varRatio*(covM[1,2:ncol(covM)])%*% solve(G2_tildeG2_tilde) %*% t(covM[1,2:ncol(covM)])
  cat("var1c: ", var1c, "\n")

if(var1c > (.Machine$double.xmin)^2){
  qtilde1c = ((Tstat1c)/sqrt(AC))/sqrt(var1c/AC) * sqrt(var2/AC) + m1
  
  if(length(NAset)/length(g) < 0.5){
    #print("OK")
    out1c = SPAtest:::Saddle_Prob(q=qtilde1c, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
  }else{
    #print("OK1")
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
    noCov = FALSE
    #if(is.null(dim(obj.null$X1))){
    #  noCov = TRUE
    #  obj.null$X1 = matrix(obj.null$X1)
    #  obj.noK$XVX_inv_XV = matrix(obj.noK$XVX_inv_XV)
    #}else{
      if(dim(obj.noK$X1)[2] == 1){
        noCov = TRUE
      }
    #}



    g1<-G0[idx_no0]
    A1<-obj.noK$XVX_inv_XV[idx_no0,]
    X1<-obj.noK$X1[idx_no0,]
    mu1<-mu[idx_no0]
    y1<-obj.noK$y[idx_no0]

    #noCov = FALSE
    #if(dim(obj.noK$X1)[2] == 1){
    # noCov = TRUE
    #}

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
  pcginvSigma<-solve(sparseSigma, g, sparse=T)
  var2 = as.matrix(t(g) %*% pcginvSigma) 
  var1 = var2 * varRatio 

}
#cat("Tstat is ", Tstat, "\n")

if(isCondition){

  T2stat = OUT_cond[,2]
  #m_all = nrow(GratioMatrixall)
 
#  cat("Tstat: ", Tstat, "\n")
  G1tilde_P_G2tilde = matrix(G1tilde_P_G2tilde,nrow=1)
  #Tstat_c = Tstat - covM[1,c(2:m_all)] %*% (solve(covM[c(2:m_all),c(2:m_all)])) %*% T2stat
  Tstat_c = Tstat - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% T2stat
  var1_c = var1 - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% t(G1tilde_P_G2tilde)

}


if(AF > 0.5){
    Tstat = (-1)*Tstat
    if(isCondition){
      Tstat_c = (-1)*Tstat_c
    }	
}

if(var1 < (.Machine$double.xmin)){
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
    if(var1_c <= (.Machine$double.xmin)){
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


scoreTest_SAIGE_binaryTrait_cond_sparseSigma=function(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y,varRatio, Cutoff, rowHeader, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

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

if(!isCondition){
  if(IsSparse==TRUE){
    if(MAF < 0.05){
       out.score<-Score_Test_Sparse(obj.noK, G0, mu.a, mu2.a, varRatio );
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

       outVec = list(BETA = out.score$BETA, SE = out.score$SE, Tstat = out.score$Tstat, p.value = out.score$pval.noadj, p.value.NA = out.score$pval.noadj, Is.converge = 1, var1 = out.score$var1, var2 = out.score$var2)
       Run1=FALSE

     }
  }
}


  if(Run1){
    G0 = matrix(G0, ncol = 1)
    XVG0 = eigenMapMatMult(obj.noK$XV, G0)
    G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
    #g = G/sqrt(AC2)
    g = G
    NAset = which(G0==0)
#    out1 = scoreTest_SPAGMMAT_binaryTrait(g, AC2, NAset, y, mu.a, varRatio, Cutoff = Cutoff)
#    if(AF > 0.5){
#      out1$BETA = (-1)*out1$BETA
#      out1$Tstat = (-1)*out1$Tstat
#    }
    out1 = scoreTest_SPAGMMAT_binaryTrait_cond_sparseSigma(g, AC2, AC,NAset, y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv)

    if(isCondition){
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"], Tstat_c = out1["Tstat_c"], p.value.c = out1["p.value.c"], var1_c = out1["var1_c"], BETA_c = out1["BETA_c"], SE_c = out1["SE_c"]) 

    }else{
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"])	
     #outVec = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
   }


  }

  return(outVec)
}


scoreTest_SPAGMMAT_binaryTrait_cond_sparseSigma=function(g, AC, AC_true, NAset, y, mu, varRatio, Cutoff, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL){

  #g = G/sqrt(AC)
  q = innerProduct(g, y)
  m1 = innerProduct(g, mu)
  Tstat = q-m1
  var2 = innerProduct(mu*(1-mu), g*g)
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
    out1 = SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
  }else{
    out1 = SPAtest:::Saddle_Prob_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="p")
  }


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
        out1_c = SPAtest:::Saddle_Prob(q=qtilde_c, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8)
      }else{
        out1_c = SPAtest:::Saddle_Prob_fast(q=qtilde_c,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="p")
      }
    #01-27-2019
    #logOR_c = (Tstat_c/var1_c)/sqrt(AC)
    logOR_c = Tstat_c/var1_c
    SE_c = abs(logOR_c/qnorm(out1_c$p.value/2))	
    out1 = c(out1, var1_c = var1_c,BETA_c = logOR_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = out1_c$p.value, p.value.NA.c = out1_c$p.value.NA) 
    }

  }


  return(out1)
}


subsetModelFileforMissing=function(obj.glmm.null, missingind, mu, mu.a, mu2.a){
        obj.glmm.null.sub = obj.glmm.null
        obj.glmm.null.sub$residuals = obj.glmm.null.sub$residuals[missingind]

        if(!is.null(obj.glmm.null.sub$P)){
                obj.glmm.null.sub$P = obj.glmm.null.sub$P[missingind, missingind]
        }
	noCov = FALSE
  if(is.null(dim(obj.glmm.null.sub$obj.noK$X1))){
	noCov = TRUE
  }else{
    if(dim(obj.glmm.null.sub$obj.noK$X1)[2] == 1){
      noCov = TRUE
    }
  }

	#if(noCov){
    	#	obj.glmm.null.sub$obj.noK$X1 = as.matrix(obj.glmm.null.sub$obj.noK$X1)
	#	obj.glmm.null.sub$obj.noK$XXVX_inv = as.matrix(obj.glmm.null.sub$obj.noK$XXVX_inv)
	#	obj.glmm.null.sub$obj.noK$XVX_inv_XV = as.matrix(obj.glmm.null.sub$obj.noK$XVX_inv_XV)
	#}

        obj.glmm.null.sub$obj.noK$X1 = obj.glmm.null.sub$obj.noK$X1[missingind,]
        obj.glmm.null.sub$obj.noK$XXVX_inv = obj.glmm.null.sub$obj.noK$XXVX_inv[missingind,]
        obj.glmm.null.sub$obj.noK$V = obj.glmm.null.sub$obj.noK$V[missingind]
        obj.glmm.null.sub$obj.noK$XV = obj.glmm.null.sub$obj.noK$XV[,missingind]
        obj.glmm.null.sub$obj.noK$y = obj.glmm.null.sub$obj.noK$y[missingind]
        obj.glmm.null.sub$obj.noK$XVX_inv_XV = obj.glmm.null.sub$obj.noK$XVX_inv_XV[missingind,]
        ##fitted values
        obj.glmm.null.sub$fitted.values = obj.glmm.null.sub$fitted.values[missingind]
        ##

	mu.sub = mu[missingind]
	mu.a.sub = mu.a[missingind]
	mu2.a.sub = mu2.a[missingind]

#        if(obj.glmm.null.sub$traitType == "binary"){
       obj.glmm.null.sub$obj.noK$XVX = t(obj.glmm.null.sub$obj.noK$X1) %*% (obj.glmm.null.sub$obj.noK$X1 *mu2.a.sub)
#        }

	obj.glmm.null.sub$obj.glm.null$y = obj.glmm.null.sub$obj.glm.null$y[missingind]
	if(!is.null(dim(obj.glmm.null.sub$obj.noK$X1))){
	  obj.glmm.null.sub$obj.noK$S_a = colSums(obj.glmm.null.sub$obj.noK$X1 * (obj.glmm.null.sub$obj.glm.null$y -  mu.a.sub))
	}else{
          obj.glmm.null.sub$obj.noK$S_a = sum(obj.glmm.null.sub$obj.noK$X1 * (obj.glmm.null.sub$obj.glm.null$y -  mu.a.sub))
        }
#	 obj.glmm.null.sub$residuals = obj.glmm.null.sub$residuals[missingind]
	if(noCov){
                obj.glmm.null.sub$obj.noK$X1 = as.matrix(obj.glmm.null.sub$obj.noK$X1)
                obj.glmm.null.sub$obj.noK$XXVX_inv = as.matrix(obj.glmm.null.sub$obj.noK$XXVX_inv)
                obj.glmm.null.sub$obj.noK$XVX_inv_XV = as.matrix(obj.glmm.null.sub$obj.noK$XVX_inv_XV)
        }

        return(subsertforMissingResult = list(obj.glmm.null.sub = obj.glmm.null.sub, mu.sub = mu.sub, mu.a.sub = mu.a.sub, mu2.a.sub = mu2.a.sub))
}

getCovMandOUT_cond_pre = function(dosage_cond, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, obj.glmm.null, sparseSigma, IsSparse=TRUE, mu, mu.a, mu2.a, Cutoff){
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
                varRatio = getVarRatio(G0, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)

                rowHeader = paste0("condMarker_",i)

                if(obj.glmm.null$traitType == "binary"){
                        out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0, AC, AF, MAF, IsSparse, obj.glmm.null$obj.noK, mu.a, mu2.a, obj.glmm.null$obj.glm.null$y, varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma)
                        OUT_cond = rbind(OUT_cond, c(as.numeric(out1$BETA), as.numeric(out1$Tstat), as.numeric(out1$var1)))

                }else if(obj.glmm.null$traitType == "quantitative"){
                        out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.glmm.null$obj.noK, AC, AF, obj.glmm.null$obj.glm.null$y, mu, varRatio, tauVec = obj.glmm.null$theta, sparseSigma=sparseSigma)
                        OUT_cond = rbind(OUT_cond, c(as.numeric(out1$BETA), as.numeric(out1$Tstat), as.numeric(out1$var1)))
                }

                OUT_cond = as.matrix(OUT_cond)
        } #end of for(i in 1:ncol(dosage_cond)){

        Mcond = ncol(dosage_cond)
        covM = matrix(0,nrow=Mcond+1, ncol = Mcond+1)

        covMsub = getCovM_nopcg(G1 = dosage_cond, G2 = dosage_cond, obj.glmm.null$obj.noK$XV, obj.glmm.null$obj.noK$XXVX_inv, sparseSigma=sparseSigma, mu2 = mu2.a)

        covM[2:(Mcond+1), 2:(Mcond+1)] = covMsub
        GratioMatrix_cond = getVarRatio(dosage_cond, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)
        G2tilde_P_G2tilde_inv = solve(covMsub * GratioMatrix_cond)

        return(condpre = list(covM = covM, OUT_cond = OUT_cond, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv))
}



getCovMandOUT_cond = function(G0, dosage_cond, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, obj.glmm.null, sparseSigma, covM, mu2.a){
        Gall = cbind(G0, dosage_cond)

        GratioMatrixall = getVarRatio(Gall, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)
        N = nrow(Gall)
        AF = sum(G0)/(2*N)
        if(AF > 0.5){
             G0_v2 = 2 - G0
        }else{
             G0_v2 = G0
        }
        G0_v2 = matrix(G0_v2, ncol=1)

        covM[1,2:ncol(covM)] = getCovM_nopcg(G1 = G0_v2, G2 = dosage_cond, obj.glmm.null$obj.noK$XV, obj.glmm.null$obj.noK$XXVX_inv, sparseSigma=sparseSigma, mu2 = mu2.a)
        G1tilde_P_G2tilde = covM[1,c(2:ncol(covM))]*(GratioMatrixall[1,c(2:ncol(covM))])
        return(condpre2 = list(covM = covM, GratioMatrixall = GratioMatrixall, G1tilde_P_G2tilde = G1tilde_P_G2tilde))
}




groupTest = function(Gmat, obj.glmm.null, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, G2_cond, G2_cond_es, kernel, method, weights.beta.rare, weights.beta.common, weightMAFcutoff, r.corr, max_maf, sparseSigma, mu.a, mu2.a, IsSingleVarinGroupTest, markerIDs, markerAFs, IsSparse, geneID, Cutoff, adjustCCratioinGroupTest, IsOutputPvalueNAinGroupTestforBinary){

        testtime <- system.time({saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.glmm.null,  cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude,ratioVec, G2_cond=G2_cond, G2_cond_es=G2_cond_es, kernel=kernel, method = method, weights.beta.rare=weights.beta.rare, weights.beta.common=weights.beta.common, weightMAFcutoff = weightMAFcutoff,  r.corr = r.corr, max_maf = max_maf, sparseSigma = sparseSigma, mu2 = mu2.a, adjustCCratioinGroupTest = adjustCCratioinGroupTest, mu = mu.a, IsOutputPvalueNAinGroupTestforBinary = IsOutputPvalueNAinGroupTestforBinary)})
	
        if(is.null(G2_cond)){
                isCondition = FALSE

        }else{
                isCondition = TRUE

        }

        OUT_single = NULL

        cat("time for SAIGE_SKAT_withRatioVec\n")
        print(testtime)
        if(length(saigeskatTest$indexNeg) > 0){
                Gmat = Gmat[,-saigeskatTest$indexNeg]
                Gmat = as.matrix(Gmat)
                markerIDs = markerIDs[-saigeskatTest$indexNeg]
                markerAFs = markerAFs[-saigeskatTest$indexNeg]
        }
        cat("saigeskatTest$p.value: ", saigeskatTest$p.value, "\n")

        if(ncol(Gmat) > 0){
	     N = nrow(Gmat)
             if(IsSingleVarinGroupTest){
		 if(obj.glmm.null$traitType == "binary"){
		   caseIndex = which(obj.glmm.null$obj.glm.null$y == 1)
		   numofCase = length(caseIndex)	
		   ctrlIndex = which(obj.glmm.null$obj.glm.null$y == 0)	
		   numofCtrl = length(ctrlIndex)
		 }
               for(nc in 1:ncol(Gmat)){
                 G0_single = Gmat[,nc]

                 AC = sum(G0_single)
                 AF = AC/(2*length(G0_single))
                 MAC = min(AC, 2*length(G0_single)-AC)
                 MAF = MAC/(2*N)
                 varRatio = getVarRatio(G0_single, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)

                 if(obj.glmm.null$traitType == "quantitative"){
                   out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0_single, obj.glmm.null$obj.noK, AC, AF, y = obj.glmm.null$obj.glm.null$y, mu = mu.a, varRatio, tauVec = obj.glmm.null$theta, sparseSigma=sparseSigma)

                  }else if(obj.glmm.null$traitType == "binary"){
		    freqinCase = sum(G0_single[caseIndex])/(2*numofCase)
		    freqinCtrl = sum(G0_single[ctrlIndex])/(2*numofCtrl)
                    out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0_single, AC, AF, MAF, IsSparse, obj.glmm.null$obj.noK, mu.a = mu.a, mu2.a = mu2.a, obj.glmm.null$obj.glm.null$y,varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma)

                  }

		outsingle = c(as.character((markerIDs)[nc]), as.numeric(AC), as.numeric((markerAFs)[nc]), as.numeric(N), as.numeric(out1$BETA), as.numeric(out1$SE), as.numeric(out1$Tstat), as.numeric(out1$p.value), as.numeric(out1$var1), as.numeric(out1$var2))

		  if(obj.glmm.null$traitType == "binary"){
			outsingle = c(outsingle, freqinCase, freqinCtrl, numofCase, numofCtrl)
		  }

		 OUT_single = rbind(OUT_single, outsingle)

                }
              }
            }# if(length(Gx$markerIDs > 0)){


            if(isCondition){
               if(saigeskatTest$m > 0){

		if(adjustCCratioinGroupTest){

			outVec = c(geneID, saigeskatTest$Out_ccadj$p.value, saigeskatTest$condOut_ccadj$p.value, saigeskatTest$markerNumbyMAC, paste(markerIDs, collapse=";"), paste(markerAFs, collapse=";"))
			if(method=="optimal.adj" & saigeskatTest$m > 0){
                        	if(saigeskatTest$m > 1){
                                	if(!is.null(saigeskatTest$Out_ccadj$p.val.each)){

                                        	p.val.vec.ccadj = saigeskatTest$Out_ccadj$param$p.val.each
                                        	rho.val.vec.ccadj = saigeskatTest$Out_ccadj$param$rho
                                        	outVec = c(outVec, p.val.vec.ccadj[which(rho.val.vec.ccadj == 1)], p.val.vec.ccadj[which(rho.val.vec.ccadj == 0)])

                                        	if(!is.null(saigeskatTest$condOut_ccadj$param$p.val.each)){
                                                	p.val.cond.vec.ccadj = saigeskatTest$condOut_ccadj$param$p.val.each
							print("p.val.cond.vec.ccadj")
							print(p.val.cond.vec.ccadj)
                                                	rho.val.cond.vec.ccadj = saigeskatTest$condOut_ccadj$param$rho
                                                	outVec = c(outVec, p.val.cond.vec.ccadj[which(rho.val.cond.vec.ccadj == 1)], p.val.cond.vec.ccadj[which(rho.val.cond.vec.ccadj == 0)])
                                        	}else{
                                                	outVec = c(outVec, 1, 1)
                                        	}

                                	}else{
                                        	outVec = c(outVec, NA, NA, NA, NA)

                                	}
                        	}else{
                                	outVec = c(outVec, NA, NA, NA, NA)
                        	}
                	}
		}
		


		if(IsOutputPvalueNAinGroupTestforBinary){
			if(!adjustCCratioinGroupTest){
				#saigeskatTest$Out_ccadj$p.value, saigeskatTest$condOut_ccadj$p.value
               			outVec = c(geneID, saigeskatTest$p.value, saigeskatTest$condOut$p.value, saigeskatTest$markerNumbyMAC, paste(markerIDs, collapse=";"), paste(markerAFs, collapse=";"))
			}else{
				outVec = c(outVec, saigeskatTest$p.value, saigeskatTest$condOut$p.value)
			}
                	if(method=="optimal.adj" & saigeskatTest$m > 0){
                        	if(saigeskatTest$m > 1){
                                	if(!is.null(saigeskatTest$param$p.val.each)){

                                        	p.val.vec = saigeskatTest$param$p.val.each
                                        	rho.val.vec = saigeskatTest$param$rho
                                        	outVec = c(outVec, p.val.vec[which(rho.val.vec == 1)], p.val.vec[which(rho.val.vec == 0)])

                                        	if(!is.null(saigeskatTest$condOut$param$p.val.each)){

                                                	p.val.cond.vec = saigeskatTest$condOut$param$p.val.each
                                                	rho.val.cond.vec = saigeskatTest$condOut$param$rho
                                                	outVec = c(outVec, p.val.cond.vec[which(rho.val.cond.vec == 1)], p.val.cond.vec[which(rho.val.cond.vec == 0)])
                                        	}else{
                                                	outVec = c(outVec, 1, 1)
                                        	}

                                	}else{
                                        	outVec = c(outVec, NA, NA, NA, NA)

                                	}
                        	}else{
                                	outVec = c(outVec, NA, NA, NA, NA)
                        	}
                	}
		}


             }else{ #end if(saigeskatTest$m > 0){
		if(adjustCCratioinGroupTest){		
                	outVec = c(geneID, NA, NA,  saigeskatTest$markerNumbyMAC, NA, NA)
			if(method=="optimal.adj"){
                        	outVec = c(outVec, NA, NA, NA, NA)
                	}

		}

		if(IsOutputPvalueNAinGroupTestforBinary){
                        if(!adjustCCratioinGroupTest){
				outVec = c(geneID, NA, NA,  saigeskatTest$markerNumbyMAC, NA, NA)
				if(method=="optimal.adj"){
                                	outVec = c(outVec, NA, NA, NA , NA)
                        	}
			}else{
				outVec = c(outVec, NA, NA)
				if(method=="optimal.adj"){
                                        outVec = c(outVec, NA, NA, NA, NA)
                                }
			}
		}	


            }

            }else{ #end of if(isCondition){


                if(saigeskatTest$m > 0){

			if(adjustCCratioinGroupTest){
	                        outVec = c(geneID, saigeskatTest$Out_ccadj$p.value, saigeskatTest$markerNumbyMAC, paste(markerIDs, collapse=";"), paste(markerAFs, collapse=";"))
        	                if(method=="optimal.adj" & saigeskatTest$m > 0){
                	                if(saigeskatTest$m > 1){
                                        	if(!is.null(saigeskatTest$Out_ccadj$param$p.val.each)){

                                                p.val.vec.ccadj = saigeskatTest$Out_ccadj$param$p.val.each
                                                rho.val.vec.ccadj = saigeskatTest$Out_ccadj$param$rho
                                                outVec = c(outVec, p.val.vec.ccadj[which(rho.val.vec.ccadj == 1)], p.val.vec.ccadj[which(rho.val.vec.ccadj == 0)])

                                        }else{
                                                outVec = c(outVec, NA, NA)

                                        }
                                }else{
                                        outVec = c(outVec, NA, NA)
                                }
                        	}
                	}

			if(IsOutputPvalueNAinGroupTestforBinary){
                        if(!adjustCCratioinGroupTest){
                                outVec = c(geneID, saigeskatTest$p.value, saigeskatTest$markerNumbyMAC, paste(markerIDs, collapse=";"), paste(markerAFs, collapse=";"))
                        }else{
                                outVec = c(outVec, saigeskatTest$p.value)
                        }
                        if(method=="optimal.adj" & saigeskatTest$m > 0){
                                if(saigeskatTest$m > 1){
                                        if(!is.null(saigeskatTest$param$p.val.each)){

                                                p.val.vec = saigeskatTest$param$p.val.each
                                                rho.val.vec = saigeskatTest$param$rho
                                                outVec = c(outVec, p.val.vec[which(rho.val.vec == 1)], p.val.vec[which(rho.val.vec == 0)])

                                        }else{
                                                outVec = c(outVec, NA, NA)

                                        }
                                }else{
                                        outVec = c(outVec, NA, NA)
                                }
                        }
                	}


	


                }else{#end of if(saigeskatTest$m > 0){
        
                if(adjustCCratioinGroupTest){
                        outVec = c(geneID, NA, saigeskatTest$markerNumbyMAC, NA, NA)
                        if(method=="optimal.adj"){
                                outVec = c(outVec, NA, NA)
                        }

                }

                if(IsOutputPvalueNAinGroupTestforBinary){
                        if(!adjustCCratioinGroupTest){
                                outVec = c(geneID, NA,  saigeskatTest$markerNumbyMAC, NA, NA)
                                if(method=="optimal.adj"){
                                        outVec = c(outVec, NA, NA)
                                }
                        }else{
                                outVec = c(outVec, NA)
                                if(method=="optimal.adj"){
                                        outVec = c(outVec, NA, NA)
                                }
                        }
                }


            }


           } # end of }else{ #end of if(isCondition){


      return(groupTestResult = list(OUT_single = OUT_single, outVec = outVec))

}





