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
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out option. By default, TRUE
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
#' @param method_to_CollapseUltraRare character. Method to collpase the ultra rare variants in the set-based association tests for BINARY traits only. This argument can be "absence_or_presence", "sum_geno", or "". absence_or_presence:  For the resulted collpased marker, any individual having dosage >= DosageCutoff_for_UltraRarePresence for any ultra rare variant has 1 in the genotype vector, otherwise 0. sum_geno: Ultra rare variants with MAC <=  MACCutoff_to_CollapseUltraRare will be collpased for set-based tests in the 'sum_geno' way and the resulted collpased marker's genotype equals weighted sum of the genotypes of all ultra rare variants. NOTE: this option sum_geno currently is NOTE active. By default, "".
#' @param MACCutoff_to_CollapseUltraRare numeric. MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10.
#' @param DosageCutoff_for_UltraRarePresence numeric. Dosage cutoff to determine whether the ultra rare variants are absent or present in the samples. Dosage >= DosageCutoff_for_UltraRarePresence indicates the varaint in present in the sample. 0< DosageCutoff_for_UltraRarePresence <= 2. By default, 0.5.
#' @return SAIGEOutputFile
#' @export
SPAGMMATtest = function(bgenFile = "",
                 bgenFileIndex = "",
		 sampleFile = "",
		 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",

		 savFile = "",
                 savFileIndex = "",

		 bedFile="",
                 bimFile="",
                 famFile="",
                 AlleleOrder = "alt-first", #new

		 idstoExcludeFile = NULL,
                 idstoIncludeFile = NULL,
                 rangestoExcludeFile = NULL,
                 rangestoIncludeFile = NULL,

		 chrom = "", #for vcf file
                 start = 1, #for vcf
                 end = 250000000, #for vcf

		 max_missing = 0.15,  #new
		 impute_method = "mean",  #"drop", "mean", "minor"     #new

		 min_MAC = 0.5,
                 min_MAF = 0,
                 min_Info = 0,
		 is_imputed_data = FALSE, #new

		 GMMATmodelFile = "",
                 LOCO=TRUE,
                 varianceRatioFile = "",
                 cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5),
                 cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5),

		 SPAcutoff=2,
                 SAIGEOutputFile = "",
                 numLinesOutput = 10,
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
                 r.corr=NULL,
		 dosage_zerod_cutoff = 0.2,
		 dosage_zerod_MAC_cutoff = 10,
		 is_output_moreDetails = FALSE, #new
                 X_PARregion="60001-2699520,154931044-155270560",   ##not activate
                 is_rewrite_XnonPAR_forMales=FALSE, #not activate
                 sampleFile_male="", #not activate

		 method_to_CollapseUltraRare="absence_or_presence",  #saige-gene+
		 MACCutoff_to_CollapseUltraRare = 10,
                 DosageCutoff_for_UltraRarePresence = 0.5,


                 function_group_test =c("lof", "missense", "synonymous"),  #new
                 maxMAFforGroupTest = c(0.1),

                 max_markers_region = 100   #new
		 ){

   if(!(impute_method %in% c("mean","minor","drop"))){
     stop("impute_method should be 'mean', 'minor', or 'drop'.")
   }


   checkArgsListBool(is_imputed_data = is_imputed_data,
                     LOCO = LOCO,
		     is_output_moreDetails = is_output_moreDetails,
		     is_rewrite_XnonPAR_forMales = is_rewrite_XnonPAR_forMales)



   checkArgsListNumeric(start = start,
                     end = end,
		     max_missing = max_missing,
                     min_MAC = min_MAC,
                     min_MAF = min_MAF,
                     min_Info = min_Info,
                     SPAcutoff = SPAcutoff,
                     numLinesOutput = numLinesOutput,
                     dosage_zerod_cutoff = dosage_zerod_cutoff,
		     dosage_zerod_MAC_cutoff = dosage_zerod_MAC_cutoff)
    if(file.exists(SAIGEOutputFile)) {print("ok -2 file exist")} 



    ##check and create the output file
    #Check_OutputFile_Create(SAIGEOutputFile)
    OutputFile = SAIGEOutputFile
    OutputFileIndex=NULL
    if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index") 

    ##check the variance ratio file and extract the variance ratio vector

    if(groupFile == ""){
      isGroupTest = FALSE
      cat("single-variant association test will be performed\n")
	  setMarker_GlobalVarsInCPP(impute_method,
				    max_missing,
                            min_MAF,
                            min_MAC,
                            min_Info,
                            1,
			    is_output_moreDetails,
			    numLinesOutput,
			    dosage_zerod_cutoff,
			    dosage_zerod_MAC_cutoff
                            )

    }else{
      isGroupTest = TRUE
      Check_File_Exist(groupFile, "groupFile")
      cat("group-based test will be performed\n")
      checkArgsList_for_Region(method_to_CollapseUltraRare,
                                    MACCutoff_to_CollapseUltraRare,
                                    DosageCutoff_for_UltraRarePresence,
                                    maxMAFforGroupTest = maxMAFforGroupTest,
				    max_markers_region = max_markers_region)




    if(file.exists(SAIGEOutputFile)) {print("ok -1 file exist")} 

        IsOutputlogPforSingle = FALSE   #to check
        OUT_Filename_Single<-sprintf("%s.single",SAIGEOutputFile)
        Check_OutputFile_Create(OUT_Filename_Single)
      if (sum(weights.beta.rare != weights.beta.common) > 0) {
        cat("WARNING:The option for weights.beta.common is not fully developed\n")
        cat("weights.beta.common is set to be equal to weights.beta.rare\n")
        weights.beta.common = weights.beta.rare
      }

      setRegion_GlobalVarsInCPP(impute_method,
                                max_missing,
				maxMAFforGroupTest,
				max_markers_region,
				1,
				method_to_CollapseUltraRare,
				MACCutoff_to_CollapseUltraRare,
				DosageCutoff_for_UltraRarePresence,
				dosage_zerod_cutoff,
                            	dosage_zerod_MAC_cutoff
                            )
     cat("dosage_zerod_cutoff is ", dosage_zerod_cutoff, "\n")
     cat("dosage_zerod_MAC_cutoff is ", dosage_zerod_MAC_cutoff, "\n")


    }
    
    ratioVec = Get_Variance_Ratio(varianceRatioFile, sparseSigmaFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest) #readInGLMM.R
  
    sparseSigmaRList = list()
    isSparseGRM = TRUE
    if(sparseSigmaFile != ""){ 
      sparseSigmaRList = setSparseSigma(sparseSigmaFile)
      isSparseGRM = TRUE
    }else{
      sparseSigmaRList = list(nSubj = 0, locations = matrix(0,nrow=2,ncol=2), values = rep(0,2))  
      isSparseGRM = FALSE 
    }	    

    obj.model = ReadModel(GMMATmodelFile, chrom, LOCO) #readInGLMM.R


    nsample = length(obj.model$y)
    cateVarRatioMaxMACVecInclude = c(cateVarRatioMaxMACVecInclude, nsample)	
    #print(names(obj.model$obj.noK))


     #in Geno.R
    objGeno = setGenoInput(bgenFile = bgenFile,
                 bgenFileIndex = bgenFileIndex,
                 vcfFile = vcfFile,   #not activate yet
                 vcfFileIndex = vcfFileIndex,
                 vcfField = vcfField,
                 savFile = savFile,
                 savFileIndex = savFileIndex,
                 sampleFile = sampleFile,
                 bedFile=bedFile,
                 bimFile=bimFile,
                 famFile=famFile,
                 idstoExcludeFile = idstoExcludeFile,
                 idstoIncludeFile = idstoIncludeFile,
                 rangestoExcludeFile = rangestoExcludeFile,
                 rangestoIncludeFile = rangestoIncludeFile,
                 chrom = chrom,
                 start = start,
                 end = end,
                 AlleleOrder = AlleleOrder,
                 sampleInModel = obj.model$sampleID)

    markerInfo = objGeno$markerInfo
    genoIndex = markerInfo$genoIndex
    genoType = objGeno$dosageFileType


    cat("isSparseGRM ", isSparseGRM, "\n")

   if (condition != "") {
        isCondition = TRUE
        #n = length(obj.model$y) #sample size
        #print(n)
        #assign_conditionMarkers_factors(genoType, condition_genoIndex,  n)
    }
    else {
        isCondition = FALSE
    }
    
    condition_genoIndex = c(-1)
    cat("isCondition ", isCondition, "\n") 
    #set up the SAIGE object based on the null model results
    setSAIGEobjInCPP(t_XVX=obj.model$obj.noK$XVX,
		     t_XXVX_inv=obj.model$obj.noK$XXVX_inv,
		     t_XV=obj.model$obj.noK$XV,
		     t_XVX_inv_XV=obj.model$obj.noK$XVX_inv_XV,
		     t_X=obj.model$X,
		     t_S_a=obj.model$obj.noK$S_a,
		     t_res=obj.model$residuals,
		     t_mu2=obj.model$mu2,
		     t_mu=obj.model$mu,
		     t_varRatio = as.vector(ratioVec),
		     t_cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
		     t_cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
		     t_SPA_Cutoff = SPAcutoff,
		     t_tauvec = obj.model$theta,
		     t_traitType = obj.model$traitType,
		     t_y = obj.model$y,
		     t_impute_method = impute_method, 
		     t_flagSparseGRM = isSparseGRM,
        	     t_locationMat = as.matrix(sparseSigmaRList$locations),
        	     t_valueVec = sparseSigmaRList$values,
        	     t_dimNum = sparseSigmaRList$nSubj, 
		     t_isCondition = isCondition,
		     t_condition_genoIndex = condition_genoIndex)


   #process condition
    if (isCondition) {
	 #print("OK1")
        n = length(obj.model$y) #sample size

        condition_genoIndex = extract_genoIndex_condition(condition, markerInfo, genoType)
	assign_conditionMarkers_factors(genoType, condition_genoIndex,  n)
	 #print("OK2")
	if(obj.model$traitType == "binary" & isGroupTest){
		outG2cond = RegionSetUpConditional_binary_InCPP()
	 #print("OK3")


	G2condList = get_newPhi_scaleFactor(q.sum = outG2cond$qsum_G2_cond, mu.a = obj.model$mu, g.sum = outG2cond$gsum_G2_cond, p.new = outG2cond$pval_G2_cond, Score = outG2cond$Score_G2_cond, Phi = outG2cond$VarMat_G2_cond)
	#print(G2condList)
	scaleFactorVec = as.vector(G2condList$scaleFactor)
	#print(scaleFactorVec)
	assign_conditionMarkers_factors_binary_region(scaleFactorVec)
	 #print("OK5")
	}	
    }


    if(file.exists(SAIGEOutputFile)) {print("ok 0 file exist")} 


    #cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
    #cat("Number of markers in each chunk:\t", numLinesOutput, "\n")
    #cat("Number of chunks for all markers:\t", nChunks, "\n")
    #}

    if(!isGroupTest){
    OutputFile = SAIGEOutputFile

    nMarkersEachChunk = numLinesOutput
    if(file.exists(SAIGEOutputFile)) {print("ok 2 file exist")}
        SAIGE.Marker(obj.model,
                   objGeno,
                   OutputFile,
                   OutputFileIndex,
                   nMarkersEachChunk,
                   is_output_moreDetails,
		   is_imputed_data,
                   LOCO,
                   chrom,
		   isCondition)
    }else{
	SAIGE.Region(obj.model,
		     objGeno,
		     sparseSigma,
		     OutputFile,
		     method_to_CollapseUltraRare,
		     MACCutoff_to_CollapseUltraRare,
                     DosageCutoff_for_UltraRarePresence,
                     groupFile,
                     function_group_test,
                     maxMAFforGroupTest,
                     max_markers_region,
                     genoType,
                     markerInfo,
		     obj.model$traitType,
		     is_imputed_data,
		     isCondition, 
		     r.corr)
    }	    
}
