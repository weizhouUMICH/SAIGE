#' Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
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
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the bgen file with NO header
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param SPAcutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param SAIGEOutputFile character. Path to the output file containing assoc test results
#' @param numLinesOutput numeric. Number of  markers to be output each time. By default, 10000   
#' @param IsSparse logical. Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for dichotomous traits. By default, TRUE 
#' @param IsOutputAFinCaseCtrl logical. Whether to output allele frequency in cases and controls. By default, FALSE
#' @param IsOutputNinCaseCtrl logical. Whether to output sample sizes in cases and controls. By default, FALSE
#' @param IsOutputHetHomCountsinCaseCtrl logical. Whether to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns "homN_Allele2_cases", "hetN_Allele2_cases", "homN_Allele2_ctrls", "hetN_Allele2_ctrls" will be output.
#' @param IsOutputlogPforSingle logical. Whether to output log(Pvalue) for single-variant assoc tests. By default, FALSE. If TRUE, the log(Pvalue) instead of original P values will be output 
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out option. By default, FALSE
#' @param condition character. For conditional analysis. Genetic marker ids (chr:pos_ref/alt if sav/vcf dosage input , marker id if bgen input) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A, Note that currently conditional analysis is only for bgen,vcf,sav input.
#' @param sparseSigmaFile character. Path to the file containing the sparseSigma from step 1. The suffix of this file is ".mtx". 
#' @param groupFile character. Path to the file containing the group information for gene-based tests. Each line is for one gene/set of variants. The first element is for gene/set name. The rest of the line is for variant ids included in this gene/set. For vcf/sav, the genetic marker ids are in the format chr:pos_ref/alt. For bgen, the genetic marker ids should match the ids in the bgen file. Each element in the line is seperated by tab. 
#' @param kernel character. For gene-based test. By default, "linear.weighted". More options can be seen in the SKAT library 
#' @param method character. method for gene-based test p-values. By default, "optimal.adj". More options can be seen in the SKAT library
#' @param weights.beta.rare vector of numeric. parameters for the beta distribution to weight genetic markers with MAF <= weightMAFcutoff in gene-based tests.By default, "c(1,25)". More options can be seen in the SKAT library
#' @param weights.beta.common vector of numeric. parameters for the beta distribution to weight genetic markers with MAF > weightMAFcutoff in gene-based tests.By default, "c(1,25)". More options can be seen in the SKAT library. NOTE: this argument is not fully developed. currently, weights.beta.common is euqal to weights.beta.rare 
#' @param weightMAFcutoff numeric. Between 0 and 0.5. See document above for weights.beta.rare and weights.beta.common. By default, 0.01
#' @param weightsIncludeinGroupFile logical. Whether to specify customized weight for makers in gene- or region-based tests. If TRUE, weights are included in the group file. For vcf/sav, the genetic marker ids and weights are in the format chr:pos_ref/alt;weight. For bgen, the genetic marker ids should match the ids in the bgen filE, e.g. SNPID;weight. Each element in the line is seperated by tab. By default, FALSE
#' @param weights_for_G2_cond vector of float. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, delimited by comma. By default, "c(1,2)" 
#' @param r.corr numeric. bewteen 0 and 1. parameters for gene-based tests.  By default, 0.  More options can be seen in the SKAT library
#' @param IsSingleVarinGroupTest logical. Whether to perform single-variant assoc tests for genetic markers included in the gene-based tests. By default, FALSE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used 
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(1.5,2.5,3.5,4.5,5.5,10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param dosageZerodCutoff numeric. In gene- or region-based tests, for each variants with MAC <= 10, dosages <= dosageZerodCutoff with be set to 0. By default, 0.2. 
#' @param IsOutputPvalueNAinGroupTestforBinary logical. In gene- or region-based tests for binary traits. if IsOutputPvalueNAinGroupTestforBinary is TRUE, p-values without accounting for case-control imbalance will be output. By default, FALSE 
#' @param IsAccountforCasecontrolImbalanceinGroupTest logical. In gene- or region-based tests for binary traits. If IsAccountforCasecontrolImbalanceinGroupTest is TRUE, p-values after accounting for case-control imbalance will be output. By default, TRUE
#' @param IsOutputBETASEinBurdenTest logical. Output effect size (BETA and SE) for burden tests. By default, FALSE
#' @param X_PARregion character. ranges of (pseudoautosomal) PAR region on chromosome X, which are seperated by comma and in the format start:end. By default: '60001-2699520,154931044-155260560' in the UCSC build hg19. For males, there are two X alleles in the PAR region, so PAR regions are treated the same as autosomes. In the NON-PAR regions (outside the specified PAR regions on chromosome X), for males, there is only one X allele. If is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in the NON-PAR regions on chromosome X will be multiplied by 2. 
#' @param is_rewrite_XnonPAR_forMales logical. Whether to rewrite gentoypes or dosages of variants in the NON-PAR regions on chromosome X for males (multiply by 2). By default, FALSE. Note, only use is_rewrite_XnonPAR_forMales=TRUE when the specified VCF or Bgen file only has variants on chromosome X. When is_rewrite_XnonPAR_forMales=TRUE, the program does not check the chromosome value by assuming all variants are on chromosome X 
#' @param sampleFile_male character. Path to the file containing one column for IDs of MALE samples in the bgen or vcf file with NO header. Order does not matter  
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
		 IsOutputHetHomCountsinCaseCtrl=FALSE,
		 IsOutputNinCaseCtrl=FALSE,
		 IsOutputlogPforSingle=FALSE,
		 LOCO=FALSE,
		 condition="",
		 sparseSigmaFile="",
		 groupFile="",
		 kernel="linear.weighted",
		 method="optimal.adj",
		 weights.beta.rare = c(1,25), 
		 weights.beta.common = c(1,25), 
		 weightMAFcutoff = 0.01,
		 weightsIncludeinGroupFile=FALSE,
		 weights_for_G2_cond = NULL, 
		 r.corr=0,
		 IsSingleVarinGroupTest = TRUE,
		 cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5), 
		 cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5),
		 dosageZerodCutoff = 0.2,	
		 IsOutputPvalueNAinGroupTestforBinary = FALSE,
		 IsAccountforCasecontrolImbalanceinGroupTest = TRUE,
		 IsOutputBETASEinBurdenTest = FALis_rewrite_XnonPAR_forMalesSE,
		 X_PARregion="60001-2699520,154931044-155270560",
		 is_rewrite_XnonPAR_forMales=FALSE,
		 sampleFile_male=""){


  if(weightMAFcutoff < 0 | weightMAFcutoff > 0.5){
    stop("weightMAFcutoff needs to be between 0 and 0.5\n")
  }

  adjustCCratioinGroupTest=TRUE
  if(!IsAccountforCasecontrolImbalanceinGroupTest){
    IsOutputPvalueNAinGroupTestforBinary = TRUE
    adjustCCratioinGroupTest = FALSE
  }

  if(sum(weights.beta.rare!=weights.beta.common) > 0){
    cat("WARNING:The option for weights.beta.common is not fully developed\n")
    cat("weights.beta.common is set to be equal to weights.beta.rare\n")
    weights.beta.common = weights.beta.rare		
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

  splitfun_weight = function(x){return(strsplit(x, split=";")[[1]][2])}
  splitfun_markerID = function(x){return(strsplit(x, split=";")[[1]][1])}

  #file for the glmm null model
  if(!file.exists(GMMATmodelFile)){
    stop("ERROR! GMMATmodelFile ", GMMATmodelFile, " does not exsit\n")
  }else{
    load(GMMATmodelFile)
    #ytemp=modglmm$y
    #modglmm$obj.glm.null = NULL
    #reduce model size
    modglmm$Y = NULL
    #modglmm$obj.glm.null = list(y=ytemp)
    modglmm$linear.predictors = NULL
    modglmm$coefficients = NULL
    modglmm$cov = NULL
    obj.glmm.null = modglmm
    rm(modglmm)
    gc(T)

    sampleInModel = NULL
    sampleInModel$IID = obj.glmm.null$sampleID
    sampleInModel = data.frame(sampleInModel)
    sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
    cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
    traitType = obj.glmm.null$traitType
    if(traitType == "quantitative"){
      IsOutputHetHomCountsinCaseCtrl = FALSE
    }

    y = obj.glmm.null$y
    X = obj.glmm.null$X
    N = length(y)
    tauVec = obj.glmm.null$theta


    indChromCheck = FALSE
    if(!LOCO){
      print("Leave-one-chromosome-out is not applied")
    }else{	    
        if(!obj.glmm.null$LOCO){
          stop("LOCO is TRUE but the null model file .rda does not contain LOCO results. In order to apply Leave-one-chromosome-out, please run Step 1 using LOCO. Otherwise, please set LOCO=FALSE in this step (Step 2).\n")
	}else{
           if(isGroupTest){ 		 
             if(chrom == ""){
               stop("chrom needs to be specified in order to apply Leave-one-chromosome-out on gene- or region-based tests")
	     }else{
	       chrom_v2 = as.character(chrom)
	       chrom_v2 = gsub("CHR", "", chrom_v2, ignore.case=T)
               chrom_v3 = as.numeric(gsub("[^0-9.]", "", chrom_v2))
               if(chrom_v3 > length(obj.glmm.null$LOCOResult) | chrom_v3 < 1){
	         stop("chromosome ", chrom, " is out of the range of null model LOCO results\n")
	       }else{	       
	         cat("Leave chromosome ", chrom_v3, " out will be applied\n")
	       } 
  	     } 
	   }else{
            if(chrom == ""){
	      if(condition != ""){
	        cat("Conditional test will be conducted and LOCO is TRUE\n")
                stop("chromosome is needed by specifying chrom for LOCO in conditioning analysis. We assume conditioning markers and testing markers are on the same chromosome")
	      }else{ 	      
                cat("WARNING: LOCO will be used, but chromosome for the dosage file is not specified. Will check each marker for its chromosome for LOCO!\n")
                indChromCheck = TRUE  
               } 
	    }else{
               chrom_v2 = as.character(chrom)
	       chrom_v2 = gsub("CHR", "", chrom_v2, ignore.case=T)
	       chrom_v3 = as.numeric(gsub("[^0-9.]", "", chrom_v2))
               if(chrom_v3 > length(obj.glmm.null$LOCOResult) | chrom_v3 < 1){
                 stop("chromosome ", chrom, " is out of the range of null model LOCO results\n")
               }else{
                 cat("Leave chromosome ", chrom_v3, " out will be applied\n")
               }
	    }	    
           
	  }		    
       }
   }#if(LOCO){

 }#if(file.exists(GMMATmodelFile)){


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


  #sample file
  if (dosageFileType == "bgen"){
    if(!file.exists(sampleFile)){
      stop("ERROR! The dosage file type is bgen but sampleFile ", sampleFile, " does not exsit\n")
    }else{
      sampleListinDosage = data.frame(data.table:::fread(sampleFile, header=F, stringsAsFactors=FALSE, colClasses=c("character")))
      sampleListinDosage$IndexDose = seq(1,nrow(sampleListinDosage), by=1)
      cat(nrow(sampleListinDosage), " sample IDs are found in sample file\n")
      colnames(sampleListinDosage)[1] = "IIDDose"
    }
  }


  if(condition != ""){
    isCondition = TRUE
  }else{
    isCondition = FALSE
  }

  cat("isCondition is ", isCondition, "\n")

   CHRv2 = NULL
   obj.model = NULL
   if(LOCO){
      if(!indChromCheck){
        if(obj.glmm.null$LOCOResult[[chrom_v3]]$isLOCO){
          obj.model = list(obj.noK = obj.glmm.null$LOCOResult[[chrom_v3]]$obj.noK, mu = as.vector(obj.glmm.null$LOCOResult[[chrom_v3]]$fitted.values))
	  #CHRv2 = chrom_v3
        }else{
	  obj.model = list(obj.noK = obj.glmm.null$obj.noK, mu  = as.vector(obj.glmm.null$fitted.values))
        }
      }
   }else{
      obj.model = list(obj.noK = obj.glmm.null$obj.noK, mu  = as.vector(obj.glmm.null$fitted.values))
   }

  if(!is.null(obj.model)){
    if(traitType == "binary"){
       obj.model$mu2 = (obj.model$mu)* (1-obj.model$mu)
    }else if(traitType == "quantitative"){	     
       obj.model$mu2 = (1/tauVec[1])*rep(1, N)
    }	     
  } 

  if(IsOutputlogPforSingle){
    cat("IsOutputlogPforSingle = TRUE. NOTE: log(Pvalue) will be output ONLY for single-variant assoc tests\n")
  }

  if (dosageFileType == "vcf"){
    vcffileopen=FALSE
    if(isCondition){ 
      isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
      sampleListinDosage_vec = getSampleIDlist_vcfMatrix()
    }else{
      if(!isGroupTest){
        setgenoTest_vcfDosage(vcfFile,vcfFileIndex,vcfField,ids_to_exclude_vcf = idstoExcludeFile, ids_to_include_vcf = idstoIncludeFile, chrom, start, end)
        isVariant = getGenoOfnthVar_vcfDosage_pre()
        sampleListinDosage_vec = getSampleIDlist() 
	vcffileopen=TRUE
      }else{
        isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField) 
        sampleListinDosage_vec = getSampleIDlist_vcfMatrix()	
      }
    }
    sampleListinDosage = data.frame(IIDDose = sampleListinDosage_vec)
    sampleListinDosage$IndexDose = seq(1,nrow(sampleListinDosage), by=1) 
    cat(nrow(sampleListinDosage), " sample IDs are found in the vcf file\n") 
  }


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
    #rm(sampleInModel)
  }

  #read in male sample IDs for assoc tests for X chromosome
  if(is_rewrite_XnonPAR_forMales){
    cat("is_rewrite_XnonPAR_forMales is TRUE, so genotypes/dosages in the non-PAR regions of X chromosome for males will be multiplied by 2\n")
    if(!file.exists(sampleFile_male)){
      stop("ERROR! The sample file for male IDs ", sampleFile_male, " does not exist\n") 
    }else{
      sampleList_male = data.frame(data.table:::fread(sampleFile_male, header=F, stringsAsFactors=FALSE, colClasses=c("character"), data.table=F))
      colnames(sampleList_male) = c("sampleID_male")
      cat(nrow(sampleList_male), " sample IDs are found in ", sampleFile_male, "\n")
      indexInModel_male = sampleInModel[sampleInModel$IID %in% (sampleList_male$sampleID_male), c("IndexInModel")]
      cat(length(indexInModel_male), " males are found in the test\n")	
      if(length(indexInModel_male) == 0){
	is_rewrite_XnonPAR_forMales=FALSE
        if(nrow(sampleList_male) > 0){
		cat("WARNING: no male IDs specified in the ", sampleFile_male, " are found sample IDs used to fit in the null model in Step 1\n")
	}	
      }else{
        cat("is_rewrite_XnonPAR_forMales=TRUE and minInfo and minMAF won't be applied to all X chromosome variants\n")
        minInfo = 0
        minMAF = 0
      }	      
    }

    X_PARregion_list = unlist(strsplit(X_PARregion, split=","))
    X_PARregion_mat = NULL
    if(length(X_PARregion_list) > 0){
      for(lxp in 1:length(X_PARregion_list)){
	X_PARregion_list_sub = as.numeric(unlist(strsplit(X_PARregion_list[lxp], split="-"))) 
        X_PARregion_mat = rbind(X_PARregion_mat, X_PARregion_list_sub)
      }
    }else{
      cat("PAR region on X chromosome is not specified\n")
    }	    
  }

  rm(sampleInModel) 
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

  if(IsDropMissingDosages){
    cat("Samples with missing dosages will be dropped from the analysis\n")
  }else{
    cat("Missing dosages will be mean imputed for the analysis\n")
  }

  setIsDropMissingDosages_bgen(IsDropMissingDosages)
  setIsDropMissingDosages_vcf(IsDropMissingDosages)
  

  ##############START TEST########################
  startTime = as.numeric(Sys.time())  # start time of the SPAGMMAT tests
  cat("Analysis started at ", startTime, "Seconds\n")

  if(minMAC == 0){
    minMAC = 0.5
    cat("As minMAC is set to be 0, minMAC = 0.5 will be used\n")
  } ##01-19-2018

  cat("minMAC: ",minMAC,"\n")
  cat("minMAF: ",minMAF,"\n")
  minMAFBasedOnMAC = minMAC/(2*N)
  testMinMAF = max(minMAFBasedOnMAC, minMAF)
  cat("Minimum MAF of markers to be tested is ", testMinMAF, "\n")

#  if(file.exists(SAIGEOutputFile)){file.remove(SAIGEOutputFile)}
#  gc(verbose=T, full=T)

  if(!isGroupTest){
    if(dosageFileType == "bgen"){
      dosageFilecolnamesSkip = c("CHR","POS","rsid","SNPID","Allele1","Allele2", "AC_Allele2", "AF_Allele2", "imputationInfo")

    }else if(dosageFileType == "vcf"){
      dosageFilecolnamesSkip = c("CHR","POS","SNPID","Allele1","Allele2", "AC_Allele2", "AF_Allele2", "imputationInfo")
    }
  }


  if(isCondition){
    condition_original=unlist(strsplit(condition,","))

    if(weightsIncludeinGroupFile){
      if(!is.null(weights_for_G2_cond)){
	#weights_for_G2_cond = unlist(strsplit(weights_for_G2_cond,","))
	if(length(weights_for_G2_cond) != length(condition_original)){
	  stop("Number of weights specified for conditioning marker(s) is different from the number of conditioning marker(s)\n")
	}	
        weights_for_G2_cond_specified = tryCatch(expr = as.numeric(weights_for_G2_cond), warning = function(w) { message("The vector is not numeric."); return(NULL)})
        if(is.null(weights_for_G2_cond_specified)){
          stop("Weights specified for conditioning marker(s) are not numeric\n")
        }		
      }else{
        stop("Weights is not specified for the conditioning marker(s)\n")
      }
    }

    if(length(condition_original) > 1){
      condition_new=NULL
      for(x in 1:length(condition_original)){
        condition_new = rbind(condition_new, c(as.numeric(strsplit(strsplit(condition_original[x], ":")[[1]][2][1], "_")[[1]][1]), condition_original[x]))
      }
      condition_new2 = condition_new[order(as.numeric(condition_new[,1])),]  

      if(weightsIncludeinGroupFile){   
        weights_for_G2_cond_specified = weights_for_G2_cond_specified[order(as.numeric(condition_new[,1]))]
	condition_specified = condition_new2[,2]
      }
      conditionlist = paste(c("condMarkers",condition_new2[,2]),collapse="\t")

    }else{
      conditionlist= paste(c("condMarkers",unlist(strsplit(condition,","))),collapse="\t") 
      if(weightsIncludeinGroupFile){
	condition_specified = unlist(strsplit(condition,","))
      }   
    }
#    conditionlist = paste(c("condMarkers",unlist(strsplit(condition,","))),collapse="\t")
    cat("conditionlist is ", conditionlist, "\n")

    if(dosageFileType == "vcf"){
      #setMAFcutoffs(0, 0.5)
      setMAFcutoffs(testMinMAF, 0.5)
      #isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
      SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N)
      Gx_cond = getGenoOfGene_vcf(conditionlist, minInfo)
      if(Gx_cond$cnt > 0){
        dosage_cond = Matrix:::sparseMatrix(i = as.vector(Gx_cond$iIndex), j = as.vector(Gx_cond$jIndex), x = as.vector(Gx_cond$dosages), symmetric = FALSE, dims = c(N, Gx_cond$cnt))
      }
    }else if(dosageFileType == "bgen"){
      SetSampleIdx(sampleIndex, N)
      Gx_cond = getGenoOfGene_bgen(bgenFile,bgenFileIndex, conditionlist, testMinMAF, 0.5, minInfo)
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
      stop("Conditioning markers are not found in the provided dosage file \n")
      isCondition = FALSE
      dosage_cond = NULL
    }else{
      if(is_rewrite_XnonPAR_forMales){
        dosage_cond = processMale_XnonPAR(indexInModel_male, dosage_cond, Gx_cond$positions, X_PARregion_mat) 
      }
    }    
  }else{#end of if(isCondition){
    dosage_cond = NULL
  }

  if(isCondition){
    if(weightsIncludeinGroupFile){
      re_index_cond = match(Gx_cond$markerIDs, condition_specified)
      weights_for_G2_cond_specified = weights_for_G2_cond_specified[re_index_cond]
      cat("Weights specified for conditioning marker(s) ", Gx_cond$markerIDs, " is ", weights_for_G2_cond_specified, "\n")
    }
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

      if(IsOutputHetHomCountsinCaseCtrl){
	resultHeader = c(resultHeader, "homN_Allele2_cases", "hetN_Allele2_cases", "homN_Allele2_ctrls", "hetN_Allele2_ctrls")
      }
	
      write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
    } #if(!isGroupTest){ 

    if(SPAcutoff < 10^-2){
      Cutoff=10^-2
    }else{
      Cutoff = SPAcutoff
    }

    #y = obj.glmm.null$y
    y1Index = which(y == 1)
    NCase = length(y1Index)
    y0Index = which(y == 0)
    NCtrl = length(y0Index)

    cat("Analyzing ", NCase, " cases and ",NCtrl, " controls \n")
    #N = length(y)
    #if(!LOCO | (LOCO & !indChromCheck)){
    #  mu2.a<-mu.a *(1-mu.a)
    #}
	    
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

    #y = obj.glmm.null$y
    #N = length(y)
    #tauVec = obj.glmm.null$theta
    #mu2.a = (1/(tauVec[1]))*rep(1, N)
    #obj.glmm.null$obj.noK$XVX = t(obj.glmm.null$obj.noK$X1) %*% (obj.glmm.null$obj.noK$X1)
    #obj.glmm.null$obj.noK$XVX_inv_XV = obj.glmm.null$obj.noK$XXVX_inv * obj.glmm.null$obj.noK$V

  }else{
    stop("ERROR! The type of the trait has to be either binary or quantitative\n")
  }



  if(nrow(varRatioData) == 1){
    cateVarRatioMinMACVecExclude=c(0)
    cateVarRatioMaxMACVecInclude=c(2*N)
  }


  if(isCondition){

    condpre = getCovMandOUT_cond_pre(dosage_cond=dosage_cond, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec, obj.model = obj.model, y = y, X = X, sparseSigma = sparseSigma, IsSparse=IsSparse, Cutoff = Cutoff, traitType = traitType, tauVec = tauVec)
    OUT_cond = condpre$OUT_cond
    G2tilde_P_G2tilde_inv = condpre$G2tilde_P_G2tilde_inv	
  }else{# end of if(isCondition)
    OUT_cond = NULL
    G2tilde_P_G2tilde_inv = NULL
  }

  cat("isCondition is ", isCondition, "\n")

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
      if(!vcffileopen){
        setgenoTest_vcfDosage(vcfFile,vcfFileIndex,vcfField,ids_to_exclude_vcf = idstoExcludeFile, ids_to_include_vcf = idstoIncludeFile, chrom, start, end)
        isVariant = getGenoOfnthVar_vcfDosage_pre()	
      } 
      SetSampleIdx_vcfDosage(sampleIndex, N)
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
        if(markerInfo >= 0 & markerInfo <= 1){
		markerInfo0 = markerInfo
	}else{
		markerInfo0 = 1
		if(markerInfo == ""){
			markerInfo = NA
		}
	}	
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
        if(markerInfo >= 0 & markerInfo <= 1){
		markerInfo0=markerInfo
	}else{	
		markerInfo0=1
		if(markerInfo ==""){
			markerInfo=NA
			Gx$variants$markerInfo=markerInfo
		}	
	}	
	#Gx$variants$markerInfo=1
        rowHeader=as.vector(unlist(Gx$variants))
        if(indChromCheck){
          CHR = Gx$variants$chromosome
	  cat("CHR ", CHR , "\n")
        }
        isVariant = getGenoOfnthVar_vcfDosage_pre()
        indexforMissing = Gx$indexforMissing
      }
  

      if(is_rewrite_XnonPAR_forMales){
	G0 = processMale_XnonPAR(indexInModel_male, G0, Gx$variants$position, X_PARregion_mat)	
      }


      MAC = min(AC, 2*N - AC)
      MAF = min(AF, 1-AF)

      if(MAF >= testMinMAF & markerInfo0 >= minInfo){
         numPassMarker = numPassMarker + 1
         varRatio = getVarRatio(G0, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)

         if(indChromCheck){
           CHR = as.character(CHR)
	   CHRv2 = gsub("CHR", "", CHR, ignore.case=T)
	   CHRv2 = as.numeric(gsub("[^0-9.]", "", CHRv2))
           if(CHRv2 > length(obj.glmm.null$LOCOResult) | CHRv2 < 1){
             stop("chromosome ", CHRv2, " is out of the range of null model LOCO results\n")
           }else{
             cat("Leave chromosome ", CHRv2, " out will be applied\n")
           }
    
           if(obj.glmm.null$LOCOResult[[CHRv2]]$isLOCO){
             obj.model = list(obj.noK = obj.glmm.null$LOCOResult[[CHRv2]]$obj.noK, mu = as.vector(obj.glmm.null$LOCOResult[[CHRv2]]$fitted.values))
      	   }else{
             obj.model = list(obj.noK = obj.glmm.null$obj.noK, mu = as.vector(obj.glmm.null$fitted.values))
           }

	   if(traitType == "binary"){
             obj.model$mu2 = (obj.model$mu) *(1-obj.model$mu)
	   }else if(traitType == "quantitative"){
	     obj.model$mu2 = (1/tau[1])*rep(1,N)
	   }	
         }


	if(IsDropMissingDosages & isCondition){
		indexforMissing = unique(c(indexforMissing, Gx_cond$indexforMissing))
	}

	if(IsOutputHetHomCountsinCaseCtrl){
		G0round = round(G0)
	}	

    if(IsDropMissingDosages & length(indexforMissing) > 0){
        missingind = seq(1, length(G0))[-(indexforMissing + 1)]
	cat("Removing ", length(indexforMissing), " samples with missing dosages/genotypes\n")

        G0 = G0[missingind]
        if(IsOutputHetHomCountsinCaseCtrl){
          G0round = G0round[missingind]
        }  
	subsetModelResult = subsetModelFileforMissing(obj.model, missingind, y, X)	
	obj.model.sub = subsetModelResult$obj.model
	#mu.a.sub = subsetModelResult$mu
	#mu.sub = mu.a.sub
        y.sub = subsetModelResult$y
	X.sub = subsetModelResult$X
	N.sub = length(G0)
	#if(traitType == "binary"){
        #	mu2.a.sub<-mu.a.sub *(1-mu.a.sub)
        #}else if(traitType == "quantitative"){
        #	mu2.a.sub = (1/tau[1])*rep(1,N.sub)
        #}

	#mu2.a.sub = subsetModelResult$mu2.a.sub
	rm(subsetModelResult)

	AC_Allele2.sub = sum(G0)
	AF_Allele2.sub = AC_Allele2.sub/(2*N.sub)
	MAF.sub = min(AF_Allele2.sub, 1-AF_Allele2.sub)

	 if(dosageFileType == "bgen"){
		rowHeader[7] = AC_Allele2.sub
		rowHeader[8] = AF_Allele2.sub

    	  }else if(dosageFileType == "vcf"){
		rowHeader[6] = AC_Allele2.sub
		rowHeader[7] = AF_Allele2.sub
    	  }




        if(traitType == "binary"){
          y1Index.sub = which(y.sub == 1)
          NCase.sub = length(y1Index.sub)
          y0Index.sub = which(y.sub == 0)
	  NCtrl.sub = length(y0Index.sub)
	  
	  if(IsOutputHetHomCountsinCaseCtrl){
	    homN_Allele2_cases = sum(G0round[y1Index.sub] == 2)
	    #print(which(G0round[y1Index.sub] == 2))
            hetN_Allele2_cases = sum(G0round[y1Index.sub] == 1) 
	    #print(which(G0round[y1Index.sub] == 1))
	    homN_Allele2_ctrls = sum(G0round[y0Index.sub] == 2)
	    #print(which(G0round[y0Index.sub] == 2))
            hetN_Allele2_ctrls = sum(G0round[y0Index.sub] == 1) 
	    #print(which(G0round[y0Index.sub] == 1))
          }	

	}

	sparseSigma.sub = sparseSigma
	if(!is.null(sparseSigma)){sparseSigma.sub = sparseSigma[missingind, missingind]}
	####Update the conditional analysis after dropping missing genotypes
        if(isCondition){
        	cat("Removing ", length(indexforMissing), " samples from the conditional marker\n")
                dosage_cond.sub = dosage_cond[missingind, ]
                dosage_cond.sub = as(dosage_cond.sub, "sparseMatrix")
                ######re-test the conditional variants after removing samples with missing genotypes
		condpre.sub = getCovMandOUT_cond_pre(dosage_cond=dosage_cond.sub, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec, obj.model = obj.model.sub, y = y.sub, X = X.sub, sparseSigma = sparseSigma.sub, IsSparse=IsSparse, Cutoff = Cutoff, traitType = traitType,tauVec=tauVec)
    		OUT_cond.sub = condpre$OUT_cond
    		G2tilde_P_G2tilde_inv.sub = condpre.sub$G2tilde_P_G2tilde_inv
		condpre2.sub = getCovMandOUT_cond(G0 = G0, dosage_cond = dosage_cond.sub, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, obj.model = obj.model.sub, sparseSigma = sparseSigma.sub, covM = condpre.sub$covM)
          	G1tilde_P_G2tilde.sub = condpre2.sub$G1tilde_P_G2tilde
          	GratioMatrixall.sub = condpre2.sub$GratioMatrixall
        }

	if(traitType == "binary"){
	  if (NCase.sub == 0 | NCtrl.sub == 0) {	
	   out1 = c(rep(NA, 8))
	   OUTvec=c(rowHeader, N.sub, unlist(out1))
	   if(IsOutputAFinCaseCtrl){
	     if(NCase.sub == 0){
		AFCase = NA
		AFCtrl = sum(G0[y0Index.sub])/(2*NCtrl.sub)
	     }else if(NCtrl.sub == 0){
		AFCtrl = NA
		AFCase = sum(G0[y1Index.sub])/(2*NCase.sub)
	     }	
	     OUTvec=c(OUTvec, AFCase, AFCtrl)
	   }

	   if(IsOutputNinCaseCtrl){
	     OUTvec=c(OUTvec, NCase.sub, NCtrl.sub)			
	   }

	   if(IsOutputHetHomCountsinCaseCtrl){
             OUTvec=c(OUTvec, homN_Allele2_cases, hetN_Allele2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls)
           }

	   OUT = rbind(OUT, OUTvec)
	   OUTvec=NULL
	  }else{ #if (NCase.sub == 0 | NCtrl.sub == 0) {
           out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0, AC_Allele2.sub, AF_Allele2.sub, MAF.sub, IsSparse, obj.model.sub$obj.noK, obj.model.sub$mu, obj.model.sub$mu2, y.sub, X.sub, varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma.sub, isCondition=isCondition, OUT_cond=OUT_cond.sub, G1tilde_P_G2tilde = G1tilde_P_G2tilde.sub, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv.sub, IsOutputlogPforSingle = IsOutputlogPforSingle)
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

	   if(IsOutputHetHomCountsinCaseCtrl){
            OUTvec=c(OUTvec, homN_Allele2_cases, hetN_Allele2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls)
           } 

	   OUT = rbind(OUT, OUTvec)
	   OUTvec=NULL
	  }
         }else if(traitType == "quantitative"){

           out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.model.sub$obj.noK, AC_Allele2.sub, AF_Allele2.sub, y.sub, X.sub, obj.model.sub$mu, varRatio, tauVec, sparseSigma=sparseSigma.sub, isCondition=isCondition, OUT_cond=OUT_cond.sub, G1tilde_P_G2tilde = G1tilde_P_G2tilde.sub, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv.sub)

           if(!isCondition){
             OUT = rbind(OUT, c(rowHeader, N.sub, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2))
           }else{
             OUT = rbind(OUT, c(rowHeader, N.sub, out1$BETA, out1$SE, out1$Tstat, out1$p.value, out1$var1, out1$var2, out1$Tstat_c,  out1$p.value.c, out1$var1_c, out1$BETA_c, out1$SE_c))
           }
         }
	
     }else{ #if(IsDropMissingDosages & length(indexforMissing) > 0){

  	  if(is_rewrite_XnonPAR_forMales){
		AC = sum(G0)
	  	AF = sum(G0)/(2*length(G0))
	        MAF = min(AF, 1-AF)	
	 	if(dosageFileType == "bgen"){
			rowHeader[7] = AC
			rowHeader[8] = AF

    	  	}else if(dosageFileType == "vcf"){
			rowHeader[6] = AC
			rowHeader[7] = AF
    	  	}	
	   }
	          ##conditional analysis
         if(isCondition){
           condpre2 = getCovMandOUT_cond(G0 = G0, dosage_cond = dosage_cond, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, obj.model = obj.model, sparseSigma = sparseSigma, covM = condpre$covM)
           G1tilde_P_G2tilde = condpre2$G1tilde_P_G2tilde
           GratioMatrixall = condpre2$GratioMatrixall

         }else{ #end of if(isCondition)
           G1tilde_P_G2tilde = NULL
           GratioMatrixall = NULL
         }	





  	 if(traitType == "binary"){
	  

           out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0, AC, AF, MAF, IsSparse, obj.model$obj.noK, obj.model$mu, obj.model$mu2, y, X, varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv, IsOutputlogPforSingle=IsOutputlogPforSingle)
	   OUTvec=c(rowHeader, N,unlist(out1))


    	   if(IsOutputAFinCaseCtrl){	     	
      	     AFCase = sum(G0[y1Index])/(2*NCase)
      	     AFCtrl = sum(G0[y0Index])/(2*NCtrl)
		OUTvec=c(OUTvec, AFCase, AFCtrl)
           }

	   if(IsOutputNinCaseCtrl){
	     OUTvec=c(OUTvec, NCase, NCtrl)
           }


	   if(IsOutputHetHomCountsinCaseCtrl){
            homN_Allele2_cases = sum(G0round[y1Index] == 2)
            hetN_Allele2_cases = sum(G0round[y1Index] == 1)
            homN_Allele2_ctrls = sum(G0round[y0Index] == 2)
            hetN_Allele2_ctrls = sum(G0round[y0Index] == 1)
	    OUTvec = c(OUTvec, homN_Allele2_cases, hetN_Allele2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls)	
          }	
	  
	   OUT = rbind(OUT, OUTvec)
	   OUTvec=NULL

         }else if(traitType == "quantitative"){

           out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.model$obj.noK, AC, AF, y, X, obj.model$mu, varRatio, tauVec, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv)

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
       #isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
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

	#obj.model$obj_cc = SKAT::SKAT_Null_Model(y ~ X-1, out_type="D", Adjustment = FALSE)
     }


       mth = 0
       MACcateNumHeader = paste0("Nmarker_MACCate_", seq(1,length(cateVarRatioMinMACVecExclude)))
       if(!isCondition){
	  if(adjustCCratioinGroupTest){	
           resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs")
	   if(method=="optimal.adj"){
	     if(IsOutputBETASEinBurdenTest){
	       resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_SKAT", "BETA_Burden", "SE_Burden")
	     }else{
	       resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_SKAT")	
	     }	
	   }
	  }

	 if(IsOutputPvalueNAinGroupTestforBinary){
           if(!adjustCCratioinGroupTest){
             resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs")
             if(method=="optimal.adj"){
	       if(IsOutputBETASEinBurdenTest){
                 resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_SKAT", "BETA_Burden", "SE_Burden")
	       }else{
		resultHeader = c("Gene", "Pvalue", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_SKAT")
	       }


             }
           }else{
             resultHeader = c(resultHeader, "Pvalue.NA")
	     if(method=="optimal.adj"){
		if(IsOutputBETASEinBurdenTest){
                  resultHeader = c(resultHeader, "Pvalue_Burden.NA","Pvalue_SKAT.NA", "BETA_Burden.NA", "SE_Burden.NA")
		}else{
		  resultHeader = c(resultHeader, "Pvalue_Burden.NA","Pvalue_SKAT.NA")	
		}

             }	
	   }
	 }
       }else{
	 if(adjustCCratioinGroupTest){
           resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs")
           if(method=="optimal.adj"){
	     if(IsOutputBETASEinBurdenTest){
             	resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_Burden_cond","Pvalue_SKAT","Pvalue_SKAT_cond", "BETA_Burden", "SE_Burden", "BETA_Burden_cond", "SE_Burden_cond")
	     }else{	
		resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs" , "Pvalue_Burden","Pvalue_Burden_cond","Pvalue_SKAT","Pvalue_SKAT_cond")
	     }	
           }
          }

	if(IsOutputPvalueNAinGroupTestforBinary){
           if(!adjustCCratioinGroupTest){
	     resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs")	
		if(method=="optimal.adj"){
		  if(IsOutputBETASEinBurdenTest){
	     	    resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs", "Pvalue_Burden","Pvalue_Burden_cond","Pvalue_SKAT","Pvalue_SKAT_cond", "BETA_Burden", "SE_Burden", "BETA_Burden_cond", "SE_Burden_cond")
		  }else{
		    resultHeader = c("Gene", "Pvalue", "Pvalue_cond", MACcateNumHeader ,"markerIDs","markerAFs", "Pvalue_Burden","Pvalue_Burden_cond","Pvalue_SKAT","Pvalue_SKAT_cond")	
		  }	
             }	
	   }else{
			resultHeader = c(resultHeader,"Pvalue.NA", "Pvalue.NA_cond")
		if(method=="optimal.adj"){
		  if(IsOutputBETASEinBurdenTest){
		    resultHeader = c(resultHeader,"Pvalue_Burden.NA","Pvalue_Burden.NA_cond","Pvalue_SKAT.NA","Pvalue_SKAT.NA_cond", "BETA_Burden.NA", "SE_Burden.NA", "BETA_Burden.NA_cond", "SE_Burden.NA_cond")
		  }else{
		    resultHeader = c(resultHeader,"Pvalue_Burden.NA","Pvalue_Burden.NA_cond","Pvalue_SKAT.NA","Pvalue_SKAT.NA_cond")
		  }
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
       cat("isCondition is ", isCondition, "\n")
       gf = file(groupFile, "r")
       while ( TRUE ) {
         marker_group_line = readLines(gf, n = 1)

         if(length(marker_group_line) == 0 ){
	   break	
         }else{
	   marker_group_line_list = strsplit(marker_group_line, split="\t")[[1]] 
	   geneID = marker_group_line_list[1]

	   cat("geneID: ", geneID, "\n")	
	   if(length(marker_group_line_list) <= 1){
	    stop("no marker IDs are found for gene ", geneID, ". Please make sure the group file is tab delimited.", "\n")
	   }
	  if(weightsIncludeinGroupFile){
		marker_group_line_list_v2 = marker_group_line_list[-1]	
		weights_specified_tmp = unlist(lapply(marker_group_line_list_v2, splitfun_weight))
		markerID_specified_tmp = unlist(lapply(marker_group_line_list_v2, splitfun_markerID))
		if(length(weights_specified_tmp) != length(markerID_specified_tmp)){stop("The length of weights is not equal to the length of markers in the group file\n")}	
		weights_specified = tryCatch(expr = as.numeric(weights_specified_tmp), warning = function(w) { message("The vector is not numeric."); return(NULL)})	
		if(is.null(weights_specified)){
			stop("Weights specified for gene ", geneID, " are not numeric\n")
		}
		marker_group_line = paste(c(geneID,markerID_specified_tmp),collapse="\t")
	  }

           if(dosageFileType == "vcf"){
             Gx = getGenoOfGene_vcf(marker_group_line, minInfo)

           }else if(dosageFileType == "bgen"){
	     print(marker_group_line)
	     cat("genetic variants with ", testMinMAF, "<= MAF <= ", maxMAFforGroupTest, "are included for gene-based tests\n") 
	     Gx = getGenoOfGene_bgen_Sparse(bgenFile,bgenFileIndex, marker_group_line, testMinMAF, maxMAFforGroupTest, minInfo)
           }
           cntMarker = Gx$cnt
           cat("cntMarker: ", cntMarker, "\n")
           if(cntMarker > 0){
		if(dosageFileType == "vcf"){
			Gmat = Matrix:::sparseMatrix(i = as.vector(Gx$iIndex), j = as.vector(Gx$jIndex), x = as.vector(Gx$dosages), symmetric = FALSE, dims = c(N, cntMarker))
		}else{
			Gmat = Matrix:::sparseMatrix(i = as.vector(Gx$iIndex), j = as.vector(Gx$jIndex), x = as.vector(Gx$dosages), symmetric = FALSE, dims = c(N, cntMarker))	
		}
		Gx$iIndex=NULL
		Gx$jIndex=NULL
		Gx$dosages=NULL

		if(is_rewrite_XnonPAR_forMales){
			Gmat = as.matrix(Gmat)
			Gmat = processMale_XnonPAR(indexInModel_male, Gmat, Gx$positions, X_PARregion_mat)
			Gmat = as(Gmat, "sparseMatrix")
		}



		if(isCondition){
			indexforMissing = unique(c(Gx$indexforMissing, Gx_cond$indexforMissing))
			#cat("indexforMissing: ", indexforMissing, "\n")
		}else{
			indexforMissing = unique(Gx$indexforMissing)

			#cat("indexforMissing: ", indexforMissing, "\n")
		}
		

		if(is_rewrite_XnonPAR_forMales | (IsDropMissingDosages & length(indexforMissing) > 0)){
	        	Gx$ACs = colSums(Gmat)
	        	Gx$markerAFs = Gx$ACs/(2*nrow(Gmat))
	        	ACtemp = 2*nrow(Gmat) - Gx$ACs
	        	Gx$MACs = pmin(Gx$ACs, ACtemp)
		}


	     if(IsDropMissingDosages & length(indexforMissing) > 0){	
		cat("Removing ", length(indexforMissing), " samples with missing dosages/genotypes in the gene\n")
		N_sub = N - length(indexforMissing) 
		cat(N_sub, " samples are left\n")
	      }else{
		N_sub = N
	      }		

	     if(N_sub > 0){

		if(IsDropMissingDosages & length(indexforMissing) > 0){

	        missingind = seq(1, nrow(Gmat))[-(indexforMissing + 1)]
		subsetModelResult = subsetModelFileforMissing(obj.model, missingind, y, X)

		y.sub = subsetModelResult$y
		X.sub = subsetModelResult$X
		obj.model.sub = subsetModelResult$obj.model
        	if(traitType == "binary"){
          		y1Index.sub = which(y.sub == 1)
          		NCase.sub = length(y1Index.sub)
          		y0Index.sub = which(y.sub == 0)
          		NCtrl.sub = length(y0Index.sub)
        	}

        	sparseSigma.sub = sparseSigma
        	if(!is.null(sparseSigma)){sparseSigma.sub = sparseSigma[missingind, missingind]}
		Gmat = Gmat[missingind,,drop = FALSE]	

		if(isCondition){
                	cat("Removing ", length(indexforMissing), " samples from the conditional marker\n")
			dosage_cond.sub = dosage_cond[missingind, , drop = FALSE]
			condpre.sub = getCovMandOUT_cond_pre(dosage_cond=dosage_cond.sub, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, ratioVec=ratioVec, obj.model = obj.model.sub, y = y.sub, X = X.sub, sparseSigma = sparseSigma.sub, IsSparse=IsSparse, Cutoff = Cutoff, traitType = traitType, tauVec=tauVec)

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
			replaceindex = which(Gmat[,z] <= dosageZerodCutoff & Gmat[,z] >0)
			if(length(replaceindex) > 0){	
				Gmat[replaceindex,z] = 0 
			}
		}
	      }	
	     cm = colMeans(Gmat)/2
	     cm[which(cm > 0.5)] = 1 - cm[which(cm > 0.5)]	
	     rmMarkerIndex = which(cm < testMinMAF | cm > maxMAFforGroupTest)	
	     if(length(rmMarkerIndex) > 0){
		cat(length(rmMarkerIndex), " marker(s) is(are) further removed\n")
		cntMarker = cntMarker - length(rmMarkerIndex)
		if(cntMarker > 0){
			Gmat = Gmat[,-rmMarkerIndex,drop = FALSE]
			print(dim(Gmat))
			Gx$markerIDs = Gx$markerIDs[-rmMarkerIndex]
			Gx$markerAFs = Gx$markerAFs[-rmMarkerIndex]
			Gmat = as(Gmat, "sparseMatrix")		
		}
	     }else{
			
		Gmat = as(Gmat, "sparseMatrix")		

	     }
	  if(weightsIncludeinGroupFile){	
	  	re_index = match(Gx$markerIDs, markerID_specified_tmp)
		weights_specified = weights_specified[re_index]
	  } 	

	    if(cntMarker > 0){
	      if(IsDropMissingDosages & length(indexforMissing) > 0){
		  cat("isCondition is ", isCondition, "\n")
	  	groupTestResult = groupTest(Gmat = Gmat, obj.model = obj.model.sub, y = y.sub, X = X.sub, tauVec, traitType, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, G2_cond = dosage_cond.sub, G2_cond_es = OUT_cond.sub[,1], kernel = kernel, method = method, weights.beta.rare = weights.beta.rare, weights.beta.common = weights.beta.common, weightMAFcutoff = weightMAFcutoff, r.corr = r.corr, max_maf = maxMAFforGroupTest, sparseSigma = sparseSigma.sub, IsSingleVarinGroupTest = IsSingleVarinGroupTest, markerIDs = Gx$markerIDs, markerAFs = Gx$markerAFs, IsSparse= IsSparse, geneID = geneID, Cutoff = Cutoff, adjustCCratioinGroupTest = adjustCCratioinGroupTest, IsOutputPvalueNAinGroupTestforBinary = IsOutputPvalueNAinGroupTestforBinary, weights_specified = weights_specified, weights_for_G2_cond = weights_for_G2_cond_specified, weightsIncludeinGroupFile = weightsIncludeinGroupFile, IsOutputBETASEinBurdenTest = IsOutputBETASEinBurdenTest, IsOutputlogPforSingle=IsOutputlogPforSingle )
	      }else{#if(IsDropMissingDosages & length(indexforMissing) > 0){	
		cat("isCondition is ", isCondition, "\n")
		groupTestResult = groupTest(Gmat = Gmat, obj.model = obj.model, y = y, X = X, tauVec, traitType,cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude, ratioVec = ratioVec, G2_cond = dosage_cond, G2_cond_es = OUT_cond[,1], kernel = kernel, method = method, weights.beta.rare = weights.beta.rare, weights.beta.common = weights.beta.common, weightMAFcutoff = weightMAFcutoff, r.corr = r.corr, max_maf = maxMAFforGroupTest, sparseSigma = sparseSigma, IsSingleVarinGroupTest = IsSingleVarinGroupTest, markerIDs = Gx$markerIDs, markerAFs = Gx$markerAFs, IsSparse= IsSparse, geneID = geneID, Cutoff = Cutoff, adjustCCratioinGroupTest = adjustCCratioinGroupTest, IsOutputPvalueNAinGroupTestforBinary = IsOutputPvalueNAinGroupTestforBinary, weights_specified = weights_specified, weights_for_G2_cond = weights_for_G2_cond_specified, weightsIncludeinGroupFile = weightsIncludeinGroupFile, IsOutputBETASEinBurdenTest = IsOutputBETASEinBurdenTest, IsOutputlogPforSingle = IsOutputlogPforSingle)

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
          }else{ ##if(cntMarker > 0){
		print("No markers are left!")
	  }
	}else{ ##if(N_sub > 0){
		print("No samples are left after removing samples with missing dosages/genotypes of variants in the gene")
	}	     
       }else{
		 print("No markers are left!")
	   }       
    } # end of while ( TRUE ) {
       }
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
scoreTest_SAIGE_quantitativeTrait=function(G0, obj.noK, AC, AF, y, X, mu, varRatio, tauVec){
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
    X<-X[idx_no0,,drop=F]
    mu1<-mu[idx_no0]
    y1<-obj.noK$y[idx_no0]
 
    noCov = FALSE
    if(dim(obj.noK$X1)[2] == 1){
     noCov = TRUE
    }

## V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv
    if(length(idx_no0) > 1){
      Z = t(A1) %*% g1
      B<-X %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      if(!noCov){
        S_a2 = obj.noK$S_a - colSums(X * (y1 - mu1))
      }else{
        S_a2 = obj.noK$S_a - crossprod(X, y1 - mu1)
      }
      #S_a2 = obj.noK$S_a - colSums(X1 * (y1 - mu1))
      S2 = -S_a2 %*% Z
    }else{
      Z = A1 * g1    
      B<-X %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      S_a2 = obj.noK$S_a - X * (y1 - mu1)
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


Score_Test_Sparse<-function(obj.null, y, X1, G, mu, mu2, varRatio, IsOutputlogPforSingle){
  # mu=mu.a; mu2= mu2.a; G=G0; obj.null=obj.noK
  idx_no0<-which(G>0)
  g1<-G[idx_no0]
  #print(length(g1))
  #noCov = FALSE
  #if(dim(obj.null$X1)[2] == 1){
  #  noCov = TRUE 
  #}
  #print("OK")
  #print(dim(X1))
  X1 = X1[idx_no0,,drop=F] 
  #print("OK")
  #print(dim(X1))
  #V = obj.null$V[idx_no0]
  #XV = obj.null$XV[,idx_no0,drop=F]
  #XV = t(X1 * V)
  #XVX = t(X1) %*% t(XV)
  #print(XVX)
  #XVX_inv = solve(XVX)
  #if(class(XVX_inv) == "try-error"){
  #	XVX_inv = ginv(XVX)
  #}
	  #else{
#	XVX_inv = solve(XVX)
 # }	  
  #XXVX_inv = X1 %*% XVX_inv
  #A1 = XXVX_inv * V 
  A1 =  obj.null$XVX_inv_XV[idx_no0,]
  #A1<-obj.null$XVX_inv_XV[idx_no0,]
  #X1<-obj.null$X1[idx_no0,]
  mu21<-mu2[idx_no0]
  mu1<-mu[idx_no0]
  y1<-y[idx_no0]

  if(length(idx_no0) > 1){
#    cat("idx_no0 ", idx_no0, "\n")
    Z = t(A1) %*% g1
    #print(dim(Z))
    B<-X1 %*% Z
    #cat("dim(Z) ", Z, "\n")
    g_tilde1 = g1 - B
    var2 = t(Z) %*% obj.null$XVX %*% Z - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
    var1 = var2 * varRatio
    S1 = crossprod(y1-mu1, g_tilde1)

    #if(!noCov){
    S_a2 = obj.null$S_a - colSums(X1 * (y1 - mu1))
    #}else{
    #  S_a2 = obj.null$S_a - crossprod(X1, y1 - mu1)
    #}

    S2 = -S_a2 %*% Z
  }else{
    Z = A1 * g1
    B<-X1 %*% Z
    g_tilde1 = g1 - B
    var2 = t(Z) %*% obj.null$XVX %*% Z - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
    var1 = var2 * varRatio
    S1 = crossprod(y1-mu1, g_tilde1)
    S_a2 = obj.null$S_a - X1 * (y1 - mu1)
    S2 = -S_a2 %*% Z
  }

  S<- S1+S2
	
  pval.noadj<-pchisq((S)^2/(var1), lower.tail = FALSE, df=1, log.p=IsOutputlogPforSingle)
  ##add on 10-25-2017
  BETA = S/var1
  if(!IsOutputlogPforSingle){
    SE = abs(BETA/qnorm(pval.noadj/2))
  }else{
    SE = abs(BETA/qnorm(exp(pval.noadj)/2))
  }
  Tstat = S
  #return(c(BETA, SE, Tstat, pval.noadj, pval.noadj, 1, var1, var2))
  return(list(BETA=BETA, SE=SE, Tstat=Tstat, pval.noadj=pval.noadj, pval.noadj=pval.noadj, is.converge=TRUE, var1=var1, var2=var2))	
}




Score_Test<-function(obj.null, G, y, mu, mu2, varRatio, IsOutputlogPforSingle){
  #print("NO SPARSE")
  g<-G  -  obj.null$XXVX_inv %*%  (obj.null$XV %*% G)
  q<-crossprod(g, y) 
  m1<-crossprod(mu, g)
  var2<-crossprod(mu2, g^2)
  var1 = var2 * varRatio
  S = q-m1
  #cat("S is ", S, "\n")
  pval.noadj<-pchisq((S)^2/var1, lower.tail = FALSE, df=1, log.p=IsOutputlogPforSingle)

  ##add on 10-25-2017
  BETA = S/var1
  if(!IsOutputlogPforSingle){
    SE = abs(BETA/qnorm(pval.noadj/2))
  }else{
    SE = abs(BETA/qnorm(exp(pval.noadj)/2))
  }
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
scoreTest_SAIGE_binaryTrait=function(G0, y, X1, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, varRatio, Cutoff, rowHeader){
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
       out.score<-Score_Test_Sparse(obj.noK, y, X1, G0,mu.a, mu2.a, varRatio );
	  #  if(is.na(as.numeric(unlist(out.score["var1"])[1]))){	
     #    out.score<-Score_Test(obj.noK, G0, y, mu.a, mu2.a, varRatio)
     #  }
     }else{
       out.score<-Score_Test(obj.noK, G0, y, mu.a, mu2.a, varRatio );
     }
     #if(out.score["pval.noadj"] > 0.05){
   if(as.numeric(unlist(out.score["var1"])[1]) <= 0){
     Run1=TRUE
   }else{	   
    if(abs(as.numeric(unlist(out.score["Tstat"])[1])/sqrt(as.numeric(unlist(out.score["var1"])[1]))) < Cutoff){
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



scoreTest_SAIGE_quantitativeTrait_sparseSigma=function(G0, obj.noK, AC, AF, y, X, mu, varRatio, tauVec, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL, IsOutputlogPforSingle=FALSE){

  N = length(G0)
  if(AF > 0.5){
    G0 = 2-G0
    AC2 = 2*N - AC
  }else{
    AC2 = AC
  }
  maf = min(AF, 1-AF)
#  cat("HERE2\n")

isSparse=FALSE
if(maf < 0.05){isSparse=TRUE}


if(isSparse){
#  cat("HERE2a\n")
    idx_no0<-which(G0>0)
    g1<-G0[idx_no0]
    X = X[idx_no0,,drop=F]
    V = obj.noK$V[idx_no0]
    XV = obj.noK$XV[,idx_no0,drop=F]
    XVX = t(X) %*% t(XV)

    
    XVX_inv = try(solve(XVX),silent=T)
    if(class(XVX_inv) == "try-error"){
    	isSparse=FALSE
    }
    if(isSparse){    
    XXVX_inv = X %*% XVX_inv
    A1 = XXVX_inv * V
    mu1<-mu[idx_no0]
    y1<-y[idx_no0]
## V = V, X1 = X1, XV = XV, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv
    if(length(idx_no0) > 1){
      Z = t(A1) %*% g1
      B<-X %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      #if(!noCov){
      S_a2 = obj.noK$S_a - colSums(X * (y1 - mu1))
      #}else{
      #  S_a2 = obj.noK$S_a - crossprod(X, y1 - mu1)
      #}
      #S_a2 = obj.noK$S_a - colSums(X1 * (y1 - mu1))
      S2 = -S_a2 %*% Z
    }else{
      Z = A1 * g1
      B<-X %*% Z
      g_tilde1 = g1 - B
      var2 = t(Z) %*% obj.noK$XVX %*% Z - sum(B^2) + sum(g_tilde1^2)
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      S_a2 = obj.noK$S_a - X * (y1 - mu1)
      S2 = -S_a2 %*% Z
    }
    S<- S1+S2
    Tstat = S/tauVec[1]
    }
}




if(!isSparse){
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
  if(!IsOutputlogPforSingle){
    p.value = 1
  }else{
    p.value=0
  }
  BETA = NA
  SE = NA
}else{
  if(!IsOutputlogPforSingle){
    p.value = pchisq(Tstat^2/var1, lower.tail = FALSE, df=1)
    BETA = (Tstat/var1)
    SE = abs(BETA/qnorm(p.value/2))
  }else{
    p.value = pchisq(Tstat^2/var1, lower.tail = FALSE, df=1, log.p=IsOutputlogPforSingle)
    BETA = (Tstat/var1)
    SE = abs(BETA/qnorm(exp(p.value)/2))
  }
}


  if(isCondition){
    if(var1_c <= (.Machine$double.xmin)){
      if(!IsOutputlogPforSingle){
        p.value.c = 1
      }else{
        p.value.c = 0
      }
      BETA_c = NA
      SE_c = NA
    }else{
      if(!IsOutputlogPforSingle){
        p.value.c = pchisq(Tstat_c^2/var1_c, lower.tail = FALSE, df=1)
        BETA_c = (Tstat_c/var1_c)
        SE_c = abs(BETA_c/qnorm(p.value.c/2))
      }else{ 
        p.value.c = pchisq(Tstat_c^2/var1_c, lower.tail = FALSE, df=1, log.p=IsOutputlogPforSingle)
        BETA_c = (Tstat_c/var1_c)
        SE_c = abs(BETA_c/qnorm(exp(p.value.c)/2))
      }
    }
  }
  
  if(isCondition){
    out1 = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2, BETA_c = BETA_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = p.value.c, var1_c = var1_c)
  }else{
    out1 = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
  }
  return(out1)
}


scoreTest_SAIGE_binaryTrait_cond_sparseSigma=function(G0, AC, AF, MAF, IsSparse, obj.noK, mu.a, mu2.a, y, X, varRatio, Cutoff, rowHeader, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL, IsOutputlogPforSingle=FALSE){

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

 if(!is.null(sparseSigma)){IsSparse=FALSE}



if(!isCondition){
  if(IsSparse==TRUE){
    if(MAF < 0.05){
       out.score<-try(Score_Test_Sparse(obj.noK, y, X, G0, mu.a, mu2.a, varRatio, IsOutputlogPforSingle=IsOutputlogPforSingle), silent=TRUE)
       if(class(out.score) == "try-error"){
           out.score<-Score_Test(obj.noK, G0, y, mu.a, mu2.a, varRatio, IsOutputlogPforSingle=IsOutputlogPforSingle)
	   #print("no sparse score here")	
       }	       
    }else{
           out.score<-Score_Test(obj.noK, G0, y, mu.a, mu2.a, varRatio, IsOutputlogPforSingle=IsOutputlogPforSingle)
    }

    if(out.score["pval.noadj"] > 0.05){

        if(as.numeric(unlist(out.score["var1"])[1]) <= 0){
          Run1=TRUE
        }else{	    
          if(abs(as.numeric(unlist(out.score["Tstat"])[1])/sqrt(as.numeric(unlist(out.score["var1"])[1]))) < Cutoff){
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
    out1 = scoreTest_SPAGMMAT_binaryTrait_cond_sparseSigma(g, AC2, AC,NAset, y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv, IsOutputlogPforSingle=IsOutputlogPforSingle)

    if(isCondition){
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"], Tstat_c = out1["Tstat_c"], p.value.c = out1["p.value.c"], var1_c = out1["var1_c"], BETA_c = out1["BETA_c"], SE_c = out1["SE_c"]) 

    }else{
     outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], var2 = out1["var2"])	
     #outVec = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
   }


  }

  return(outVec)
}


scoreTest_SPAGMMAT_binaryTrait_cond_sparseSigma=function(g, AC, AC_true, NAset, y, mu, varRatio, Cutoff, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL, IsOutputlogPforSingle=FALSE){

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
    out1 = SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8, log.p=IsOutputlogPforSingle)
  }else{
    out1 = SPAtest:::Saddle_Prob_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="p", log.p=IsOutputlogPforSingle)
  }


  out1$var1 = var1
  out1$var2 = var2

  #01-27-2019
  #as g is not divided by sqrt(AC), the sqrt(AC) is removed from the denominator
  #logOR = (Tstat/var1)/sqrt(AC)
  logOR = Tstat/var1
  if(!IsOutputlogPforSingle){
    SE = abs(logOR/qnorm(out1$p.value/2))
  }else{
    SE = abs(logOR/qnorm(exp(out1$p.value)/2))
  }
#  out1 = c(out1, BETA = logOR, SE = SE, Tstat = Tstat)
  out1$BETA=logOR
  out1$SE=SE
  out1$Tstat = Tstat

  if(isCondition){
    if(var1_c <= (.Machine$double.xmin)^2){
      if(!IsOutputlogPforSingle){
        out1 = c(out1, var1_c = var1_c,BETA_c = NA, SE_c = NA, Tstat_c = Tstat_c, p.value.c = 1, p.value.NA.c = 1)	
      }else{
        out1 = c(out1, var1_c = var1_c,BETA_c = NA, SE_c = NA, Tstat_c = Tstat_c, p.value.c = 0, p.value.NA.c = 0)
      }

    }else{

      qtilde_c = Tstat_c/sqrt(var1_c) * sqrt(var2) + m1
      if(length(NAset)/length(g) < 0.5){
        out1_c = SPAtest:::Saddle_Prob(q=qtilde_c, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8, log.p=IsOutputlogPforSingle)
      }else{
        out1_c = SPAtest:::Saddle_Prob_fast(q=qtilde_c,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="p", log.p=IsOutputlogPforSingle)
      }
    #01-27-2019
    #logOR_c = (Tstat_c/var1_c)/sqrt(AC)
    logOR_c = Tstat_c/var1_c
    if(!IsOutputlogPforSingle){
      SE_c = abs(logOR_c/qnorm(out1_c$p.value/2))	
    }else{
      SE_c = abs(logOR_c/qnorm(exp(out1_c$p.value)/2))	
    }

    out1 = c(out1, var1_c = var1_c,BETA_c = logOR_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = out1_c$p.value, p.value.NA.c = out1_c$p.value.NA) 
    }

  }


  return(out1)
}



subsetModelFileforMissing=function(obj.model, missingind, y, X){

        y = y[missingind]
        X = X[missingind,,drop=FALSE]

        mu = obj.model$mu[missingind]
	mu2 = obj.model$mu2[missingind]
	obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)

        return(list(obj.model = list(obj.noK = obj.noK, mu = mu, mu2 = mu2), y = y, X = X))
        #return(subsertforMissingResult = list(obj.glmm.null.sub = obj.glmm.null.sub, mu.sub = mu.sub, mu.a.sub = mu.a.sub, mu2.a.sub = mu2.a.sub))
}


        #obj.glmm.null.sub$obj.noK$V = obj.glmm.null.sub$obj.noK$V[missingind]
        #obj.glmm.null.sub$obj.noK$XV = obj.glmm.null.sub$obj.noK$XV[,missingind, drop=FALSE]
	#obj.glmm.null.sub$obj.noK$XV = t(obj.glmm.null.sub$obj.noK$X1 * obj.glmm.null.sub$obj.noK$V)
        #obj.glmm.null.sub$obj.noK$XVX_inv_XV_old = obj.glmm.null.sub$obj.noK$XVX_inv_XV[missingind,]
        ##fitted values
	#obj.glmm.null.sub$obj.noK$XVX = t(obj.glmm.null.sub$obj.noK$X1)  %*% t(obj.glmm.null.sub$obj.noK$XV)
	#obj.glmm.null.sub$obj.noK$XVX_inv = solve(obj.glmm.null.sub$obj.noK$XVX)
	#obj.glmm.null.sub$obj.noK$XXVX_inv = obj.glmm.null.sub$obj.noK$X1 %*% obj.glmm.null.sub$obj.noK$XVX_inv 
	#obj.glmm.null.sub$obj.noK$XVX_inv_XV = obj.glmm.null.sub$obj.noK$XXVX_inv * obj.glmm.null.sub$obj.noK$V
	#print("XVX_inv_XV_old")
	#print(obj.glmm.null.sub$obj.noK$XVX_inv_XV_old)
	#print("XVX_inv_XV")
	#print(obj.glmm.null.sub$obj.noK$XVX_inv_XV)
	#print("XXVX_inv_old")
	#print(obj.glmm.null.sub$obj.noK$XXVX_inv_old)
	#print("XXVX_inv")
	#print(obj.glmm.null.sub$obj.noK$XXVX_inv)
        ##

	#mu.sub = mu[missingind]
	#mu.a.sub = mu.a[missingind]
	#mu2.a.sub = mu2.a[missingind]
	#mu2.a.sub = (1-mu.a.sub)*mu.a.sub
        #if(obj.glmm.null.sub$traitType == "binary"){
       #obj.glmm.null.sub$obj.noK$XVX = t(obj.glmm.null.sub$obj.noK$X1) %*% (obj.glmm.null.sub$obj.noK$X1 *mu2.a.sub)
        #}

	#obj.glmm.null.sub$obj.glm.null$y = obj.glmm.null.sub$obj.glm.null$y[missingind]
	#if(!is.null(dim(obj.glmm.null.sub$obj.noK$X1))){
	#  obj.glmm.null.sub$obj.noK$S_a = colSums(obj.glmm.null.sub$obj.noK$X1 * (obj.glmm.null.sub$obj.glm.null$y -  mu.a.sub))
	#}else{
        #  obj.glmm.null.sub$obj.noK$S_a = sum(obj.glmm.null.sub$obj.noK$X1 * (obj.glmm.null.sub$obj.glm.null$y -  mu.a.sub))
        #}
#	 obj.glmm.null.sub$residuals = obj.glmm.null.sub$residuals[missingind]
	#if(noCov){
        #        obj.glmm.null.sub$obj.noK$X1 = as.matrix(obj.glmm.null.sub$obj.noK$X1)
        #        obj.glmm.null.sub$obj.noK$XXVX_inv = as.matrix(obj.glmm.null.sub$obj.noK$XXVX_inv)
        #        obj.glmm.null.sub$obj.noK$XVX_inv_XV = as.matrix(obj.glmm.null.sub$obj.noK$XVX_inv_XV)
        #}


getCovMandOUT_cond_pre = function(dosage_cond, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, obj.model, y, X, sparseSigma, IsSparse=TRUE, Cutoff, traitType, tauVec=tauVec){
        OUT_cond = NULL
        for(i in 1:ncol(dosage_cond)){
                G0  = dosage_cond[,i]
                AC = sum(G0)
                N  = length(G0)
                AF = AC/(2*N)
                MAF = min(AF, 1-AF)
                MAC = min(AC, 2*N - AC)

                varRatio = getVarRatio(G0, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)

                rowHeader = paste0("condMarker_",i)

                if(traitType == "binary"){
                        out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0, AC, AF, MAF, IsSparse, obj.model$obj.noK, obj.model$mu, obj.model$mu2, y, X, varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma)
                        OUT_cond = rbind(OUT_cond, c(as.numeric(out1$BETA), as.numeric(out1$Tstat), as.numeric(out1$var1)))

                }else if(traitType == "quantitative"){
                        out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0, obj.model$obj.noK, AC, AF, y, X, obj.model$mu,  varRatio, tauVec = tauVec, sparseSigma=sparseSigma)
                        OUT_cond = rbind(OUT_cond, c(as.numeric(out1$BETA), as.numeric(out1$Tstat), as.numeric(out1$var1)))
                }

                OUT_cond = as.matrix(OUT_cond)
        } #end of for(i in 1:ncol(dosage_cond)){

        Mcond = ncol(dosage_cond)
        covM = matrix(0,nrow=Mcond+1, ncol = Mcond+1)

        covMsub = getCovM_nopcg(G1 = dosage_cond, G2 = dosage_cond, obj.model$obj.noK$XV, obj.model$obj.noK$XXVX_inv, sparseSigma=sparseSigma, mu2 = obj.model$mu2)

        covM[2:(Mcond+1), 2:(Mcond+1)] = covMsub
        GratioMatrix_cond = getVarRatio(dosage_cond, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)
        G2tilde_P_G2tilde_inv = solve(covMsub * GratioMatrix_cond)

        return(condpre = list(covM = covM, OUT_cond = OUT_cond, G2tilde_P_G2tilde_inv = G2tilde_P_G2tilde_inv))
}



getCovMandOUT_cond = function(G0, dosage_cond, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, obj.model, sparseSigma, covM){
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

        covM[1,2:ncol(covM)] = getCovM_nopcg(G1 = G0_v2, G2 = dosage_cond, obj.model$obj.noK$XV, obj.model$obj.noK$XXVX_inv, sparseSigma=sparseSigma, mu2 = obj.model$mu2)
        G1tilde_P_G2tilde = covM[1,c(2:ncol(covM))]*(GratioMatrixall[1,c(2:ncol(covM))])
        return(condpre2 = list(covM = covM, GratioMatrixall = GratioMatrixall, G1tilde_P_G2tilde = G1tilde_P_G2tilde))
}




groupTest = function(Gmat, obj.model, y, X, tauVec, traitType, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, G2_cond, G2_cond_es, kernel, method, weights.beta.rare, weights.beta.common, weightMAFcutoff, r.corr, max_maf, sparseSigma, IsSingleVarinGroupTest, markerIDs, markerAFs, IsSparse, geneID, Cutoff, adjustCCratioinGroupTest, IsOutputPvalueNAinGroupTestforBinary, weights_specified, weights_for_G2_cond, weightsIncludeinGroupFile, IsOutputBETASEinBurdenTest, IsOutputlogPforSingle=FALSE){
	obj.model$theta = tauVec
	obj.model$residuals = as.vector(y-obj.model$mu)
        testtime <- system.time({saigeskatTest = SAIGE_SKAT_withRatioVec(Gmat, obj.model, y, X, tauVec, cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude,ratioVec, G2_cond=G2_cond, G2_cond_es=G2_cond_es, kernel=kernel, method = method, weights.beta.rare=weights.beta.rare, weights.beta.common=weights.beta.common, weightMAFcutoff = weightMAFcutoff,  r.corr = r.corr, max_maf = max_maf, sparseSigma = sparseSigma, mu2 = obj.model$mu2, adjustCCratioinGroupTest = adjustCCratioinGroupTest, mu = obj.model$mu, IsOutputPvalueNAinGroupTestforBinary = IsOutputPvalueNAinGroupTestforBinary, weights_specified = weights_specified, weights_for_G2_cond = weights_for_G2_cond, weightsIncludeinGroupFile = weightsIncludeinGroupFile, IsOutputBETASEinBurdenTest=IsOutputBETASEinBurdenTest)})
	
        if(is.null(G2_cond)){
                isCondition = FALSE

        }else{
                isCondition = TRUE

        }

        OUT_single = NULL

        cat("time for SAIGE_SKAT_withRatioVec\n")
        print(testtime)
        if(length(saigeskatTest$indexNeg) > 0){
                #Gmat = Gmat[,-saigeskatTest$indexNeg]
		Gmat = array(a, dim = c(nrow(Gmat), ncol(Gmat)))[,-saigeskatTest$indexNeg, drop=F]
                #Gmat = Gmat[,-saigeskatTest$indexNeg]
                Gmat = as.matrix(Gmat)
                markerIDs = markerIDs[-saigeskatTest$indexNeg]
                markerAFs = markerAFs[-saigeskatTest$indexNeg]
        }
        cat("saigeskatTest$p.value: ", saigeskatTest$p.value, "\n")
	print("OK1")

        if(ncol(Gmat) > 0){
	     N = nrow(Gmat)
             if(IsSingleVarinGroupTest){
		 if(traitType == "binary"){
		   caseIndex = which(y == 1)
		   numofCase = length(caseIndex)	
		   ctrlIndex = which(y == 0)	
		   numofCtrl = length(ctrlIndex)
		 }
               for(nc in 1:ncol(Gmat)){
                 G0_single = Gmat[,nc]

                 AC = sum(G0_single)
                 AF = AC/(2*length(G0_single))
                 MAC = min(AC, 2*length(G0_single)-AC)
                 MAF = MAC/(2*N)
                 varRatio = getVarRatio(G0_single, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec)

                 if(traitType == "quantitative"){
                   out1 = scoreTest_SAIGE_quantitativeTrait_sparseSigma(G0_single, obj.model$obj.noK, AC, AF, y = y, X=X,  mu = obj.model$mu, varRatio, tauVec = tauVec, sparseSigma=sparseSigma, IsOutputlogPforSingle)

                  }else if(traitType == "binary"){
		    freqinCase = sum(G0_single[caseIndex])/(2*numofCase)
		    freqinCtrl = sum(G0_single[ctrlIndex])/(2*numofCtrl)
                    out1 = scoreTest_SAIGE_binaryTrait_cond_sparseSigma(G0_single, AC, AF, MAF, IsSparse, obj.model$obj.noK, mu.a = obj.model$mu, mu2.a = obj.model$mu2, y, X, varRatio, Cutoff, rowHeader, sparseSigma=sparseSigma, IsOutputlogPforSingle=IsOutputlogPforSingle)

                  }

		outsingle = c(as.character((markerIDs)[nc]), as.numeric(AC), as.numeric((markerAFs)[nc]), as.numeric(N), as.numeric(out1$BETA), as.numeric(out1$SE), as.numeric(out1$Tstat), as.numeric(out1$p.value), as.numeric(out1$var1), as.numeric(out1$var2))

		  if(traitType == "binary"){
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
                                	if(!is.null(saigeskatTest$Out_ccadj$param$p.val.each)){

                                        	p.val.vec.ccadj = saigeskatTest$Out_ccadj$param$p.val.each
                                        	rho.val.vec.ccadj = saigeskatTest$Out_ccadj$param$rho
                                        	outVec = c(outVec, p.val.vec.ccadj[which(rho.val.vec.ccadj == 1)], p.val.vec.ccadj[which(rho.val.vec.ccadj == 0)])
						if(IsOutputBETASEinBurdenTest){
							BETA_Burden = saigeskatTest$Score_sum/(saigeskatTest$Phi_ccadj_sum)
							SE_Burden = abs(BETA_Burden/qnorm(p.val.vec.ccadj[which(rho.val.vec.ccadj == 1)]/2)) 
						}


                                        	if(!is.null(saigeskatTest$condOut_ccadj$param$p.val.each)){
                                                	p.val.cond.vec.ccadj = saigeskatTest$condOut_ccadj$param$p.val.each
                                                	rho.val.cond.vec.ccadj = saigeskatTest$condOut_ccadj$param$rho
                                                	outVec = c(outVec, p.val.cond.vec.ccadj[which(rho.val.cond.vec.ccadj == 1)], p.val.cond.vec.ccadj[which(rho.val.cond.vec.ccadj == 0)])
							if(IsOutputBETASEinBurdenTest){	
								BETA_Burden_cond = saigeskatTest$Score_cond_ccadj_sum/(saigeskatTest$Phi_cond_ccadj_sum)
								SE_Burden_cond = abs(BETA_Burden_cond/qnorm(p.val.cond.vec.ccadj[which(rho.val.cond.vec.ccadj == 1)]/2))
								outVec = c(outVec, BETA_Burden, SE_Burden, BETA_Burden_cond, SE_Burden_cond)

							}
                                        	}else{
                                                	outVec = c(outVec, 1, 1)
							if(IsOutputBETASEinBurdenTest){
								outVec = c(outVec, BETA_Burden, SE_Burden, NA, NA)
							}

                                        	}

                                	}else{
                                        	outVec = c(outVec, NA, NA, NA, NA)
						if(IsOutputBETASEinBurdenTest){
							outVec = c(outVec, NA, NA, NA, NA)
						}

                                	}
                        	}else{
                                	outVec = c(outVec, NA, NA, NA, NA)
					if(IsOutputBETASEinBurdenTest){
						outVec = c(outVec, NA, NA, NA, NA)
					}
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

						if(IsOutputBETASEinBurdenTest){
                                                        BETA_Burden = saigeskatTest$Score_sum/(saigeskatTest$Phi_sum)
                                                        SE_Burden = abs(BETA_Burden/qnorm(p.val.vec[which(rho.val.vec == 1)]/2))
								
                                                }


                                        	if(!is.null(saigeskatTest$condOut$param$p.val.each)){

                                                	p.val.cond.vec = saigeskatTest$condOut$param$p.val.each
                                                	rho.val.cond.vec = saigeskatTest$condOut$param$rho
                                                	outVec = c(outVec, p.val.cond.vec[which(rho.val.cond.vec == 1)], p.val.cond.vec[which(rho.val.cond.vec == 0)])
							if(IsOutputBETASEinBurdenTest){	
								BETA_Burden_cond =  saigeskatTest$Score_cond_sum/(saigeskatTest$Phi_cond_sum)
                                                        	SE_Burden_cond = abs(BETA_Burden_cond/qnorm(p.val.cond.vec[which(rho.val.cond.vec == 1)]/2))
								outVec = c(outVec, BETA_Burden, SE_Burden, BETA_Burden_cond, SE_Burden_cond)	
							}	
                                        	}else{
                                                	outVec = c(outVec, 1, 1)
							if(IsOutputBETASEinBurdenTest){
                                                                outVec = c(outVec, BETA_Burden, SE_Burden, NA, NA)
                                                        }
                                        	}

                                	}else{
                                        	outVec = c(outVec, NA, NA, NA, NA)
						if(IsOutputBETASEinBurdenTest){
                                                	outVec = c(outVec, NA, NA, NA, NA)
                                                }
                                	}
                        	}else{
                                	outVec = c(outVec, NA, NA, NA, NA)
					if(IsOutputBETASEinBurdenTest){
                                                outVec = c(outVec, NA, NA, NA, NA)
                                        }
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
						if(IsOutputBETASEinBurdenTest){
							BETA_Burden = saigeskatTest$Score_sum/(saigeskatTest$Phi_ccadj_sum)
							SE_Burden = abs(BETA_Burden/qnorm(p.val.vec.ccadj[which(rho.val.vec.ccadj == 1)]/2)) 
							outVec = c(outVec, BETA_Burden, SE_Burden)
						}
                                        }else{
                                                outVec = c(outVec, NA, NA)
						if(IsOutputBETASEinBurdenTest){
							outVec = c(outVec, NA, NA)
						}

                                        }
                                }else{
                                        outVec = c(outVec, NA, NA)
					if(IsOutputBETASEinBurdenTest){
						outVec = c(outVec, NA, NA)
					}
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
						if(IsOutputBETASEinBurdenTest){
							BETA_Burden.NA = saigeskatTest$Score_sum/(saigeskatTest$Phi_sum)
							SE_Burden.NA = abs(BETA_Burden.NA/qnorm(p.val.vec[which(rho.val.vec == 1)]/2)) 
							outVec = c(outVec, BETA_Burden.NA, SE_Burden.NA)
						}

                                        }else{
                                                outVec = c(outVec, NA, NA)
						if(IsOutputBETASEinBurdenTest){
							outVec = c(outVec, NA, NA)
						}

                                        }
                                }else{
                                        outVec = c(outVec, NA, NA)
					if(IsOutputBETASEinBurdenTest){
						outVec = c(outVec, NA, NA)
					}
                                }
                        }
                	}


	


                }else{#end of if(saigeskatTest$m > 0){
        
                if(adjustCCratioinGroupTest){
                        outVec = c(geneID, NA, saigeskatTest$markerNumbyMAC, NA, NA)
                        if(method=="optimal.adj"){
                                outVec = c(outVec, NA, NA)
				if(IsOutputBETASEinBurdenTest){
					outVec = c(outVec, NA, NA)
				}
                        }

                }

                if(IsOutputPvalueNAinGroupTestforBinary){
                        if(!adjustCCratioinGroupTest){
                                outVec = c(geneID, NA,  saigeskatTest$markerNumbyMAC, NA, NA)
                                if(method=="optimal.adj"){
                                        outVec = c(outVec, NA, NA)
					if(IsOutputBETASEinBurdenTest){
						outVec = c(outVec, NA, NA)
					}
                                }
                        }else{
                                outVec = c(outVec, NA)
                                if(method=="optimal.adj"){
                                        outVec = c(outVec, NA, NA)
					if(IsOutputBETASEinBurdenTest){
						outVec = c(outVec, NA, NA)
					}
                                }
                        }
                }


            }


           } # end of }else{ #end of if(isCondition){


      return(groupTestResult = list(OUT_single = OUT_single, outVec = outVec))

}


processMale_XnonPAR = function(maleIDindex, Gx, positionL, XPARregion){
	print(positionL)
	for(i in 1:length(positionL)){
		inPAR = FALSE
		if(!is.null(XPARregion)){
		  for (j in 1:nrow(XPARregion)){
		    if(inPAR == FALSE){	
			if(positionL[i] <= XPARregion[j,2] & positionL[i] >= XPARregion[j,1]){
				inPAR = TRUE
			}
		    }
		  }
		}
		if(!inPAR){
			if(length(positionL) > 1){
			  Gx[maleIDindex,i] = 2*Gx[maleIDindex,i] 
			}else{
			  Gx[maleIDindex] = 2*Gx[maleIDindex]	
			}	
		}	
	}
	return(Gx)	
}	


