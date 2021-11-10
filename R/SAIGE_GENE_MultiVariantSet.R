#' Run gene- or region-based score tests with SPA based on the linear/logistic mixed model.
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
#' @param IsOutputMAFinCaseCtrlinGroupTest logical. Whether to output minor allele frequency in cases and controls in set-based tests By default, FALSE
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
#' @param method_to_CollapseUltraRare character. Method to collpase the ultra rare variants in the set-based association tests. This argument can be 'absence_or_presence', 'sum_geno', or ''. absence_or_presence:  For the resulted collpased marker, any individual having DosageCutoff_for_UltraRarePresence <= dosage < 1+DosageCutoff_for_UltraRarePresence for any ultra rare variant has 1 in the genotype vector, having dosage >= 1+DosageCutoff_for_UltraRarePresence for any ultra rare variant has 2 in the genotype vector, otherwise 0. sum_geno: Ultra rare variants with MAC <=  MACCutoff_to_CollapseUltraRare will be collpased for set-based tests in the 'sum_geno' way and the resulted collpased marker's genotype equals weighted sum of the genotypes of all ultra rare variants. NOTE: this option sum_geno currently is NOT active. By default, "absence_or_presence". 
#' @param MACCutoff_to_CollapseUltraRare numeric. MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10. 
#' @param DosageCutoff_for_UltraRarePresence numeric. Dosage cutoff to determine whether the ultra rare variants are absent or present in the samples. Dosage >= DosageCutoff_for_UltraRarePresence indicates the varaint in present in the sample. 0< DosageCutoff_for_UltraRarePresence <= 2. By default, 0.5.  
#' @return SAIGEOutputFile
#' @export
SAIGE_GENE_MultiVariantSets = function(bgenFile = "",
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
		 maxMAFforGroupTest = 0.5,
     minInfo = 0,
     GMMATmodelFile = "",
     varianceRatioFile = "",
     SPAcutoff=2,
     SAIGEOutputFile = "",
		 IsSparse=TRUE,
		 LOCO=TRUE,
		 sparseSigmaFile="",
		 groupFile="",
		 kernel="linear.weighted",
		 method="optimal.adj",
		 weights.beta.rare = c(1,25),
		 weights.beta.common=c(1,25),
		 weightMAFcutoff = 0.01,
		 r.corr=0,
		 IsSingleVarinGroupTest = TRUE,
		 cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5),
		 cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5),
		 dosageZerodCutoff = 0.2,
		 method_to_CollapseUltraRare="absence_or_presence",
		 MACCutoff_to_CollapseUltraRare = 10,
		 DosageCutoff_for_UltraRarePresence = 0.5,
		 function_group_test =c("lof", "missense", "synonymous"),
		 MAF_cutoff=c(0.0001,0.001,0.01)){

		
		# weights_specified is currently not used 
		weights_specified=NULL;
		 
	  # currently prespecified 
	  IsDropMissingDosages =FALSE
	
	
    if (weightMAFcutoff < 0 | weightMAFcutoff > 0.5) {
        stop("weightMAFcutoff needs to be between 0 and 0.5\n")
    }
    adjustCCratioinGroupTest = TRUE
    if (dosageZerodCutoff < 0) {
        dosageZerodCutoff = 0
    } else if (dosageZerodCutoff >= 0) {
        cat("Any dosages <= ", dosageZerodCutoff, " for genetic variants with MAC <= 10 are set to be 0 in group tests\n")
    }
    
    if (file.exists(SAIGEOutputFile)) {
        file.remove(SAIGEOutputFile)
    }
    if (!file.exists(SAIGEOutputFile)) {
        file.create(SAIGEOutputFile, showWarnings = TRUE)
    }
    
    # Check file existence
    Check_File_Exist(groupFile, "groupFile")
    Check_File_Exist(GMMATmodelFile, "GMMATmodelFile")
    
    # load GMMATmodelFile
    load(GMMATmodelFile)
    modglmm$Y = NULL
    modglmm$offset = modglmm$linear.predictors - modglmm$coefficients[1]
    modglmm$linear.predictors = NULL
    modglmm$coefficients = NULL
    modglmm$cov = NULL
    obj.glmm.null = modglmm
    rm(modglmm)
    gc(T)
        

    traitType = obj.glmm.null$traitType
    y = obj.glmm.null$y
    X = obj.glmm.null$X
    N = length(y)
    tauVec = obj.glmm.null$theta
    
    indChromCheck = FALSE
    
    # check genotypes and samples 
 	  obj_geno_sample = Check_Genotypes_and_Samples(bgenFile, vcfFile, vcfField, vcfFileIndex, savFile, sampleFile, obj.glmm.null$sampleID, chrom)
 	  dosageFileType = obj_geno_sample$dosageFileType
 	  vcfFile = obj_geno_sample$vcfFile
 	  sampleIndex = obj_geno_sample$sampleIndex
 	


 	
  	# get variance ratio 
   	ratioVec = Get_Variance_Ratio(varianceRatioFile, sparseSigmaFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)


    # Check loco
    if (!LOCO) {
        print("Leave-one-chromosome-out is not applied")
        if (obj.glmm.null$LOCO) {
            for (chr in 1:22) {
                obj.glmm.null$LOCOResult[chr] = list(NULL)
                cat("chromosome ", chr, " model results are removed to save memory\n")
                gc()
            }
        }
    }else {
        if (!obj.glmm.null$LOCO) {
            stop("LOCO is TRUE but the null model file .rda does not contain LOCO results. In order to apply Leave-one-chromosome-out, please run Step 1 using LOCO. Otherwise, please set LOCO=FALSE in this step (Step 2).\n")
        }else {
            if (chrom == "") {
                stop("chrom needs to be specified in order to apply Leave-one-chromosome-out on gene- or region-based tests")
            }else {
                chrom_v2 = as.character(chrom)
                chrom_v2 = gsub("CHR", "", chrom_v2, ignore.case = T)
                chrom_v3 = as.numeric(gsub("[^0-9.]", "", chrom_v2))
                if (chrom_v3 > length(obj.glmm.null$LOCOResult) | chrom_v3 < 1) {
                    stop("chromosome ", chrom, " is out of the range of null model LOCO results\n")
                }
                else {
                    cat("Leave chromosome ", chrom_v3, " out will be applied\n")
                }
            }
        }
    }


    
    CHRv2 = NULL
    obj.model = NULL
    if (LOCO) {
        if (obj.glmm.null$LOCOResult[[chrom_v3]]$isLOCO) {
            obj.model = list(obj.noK = obj.glmm.null$LOCOResult[[chrom_v3]]$obj.noK, mu = as.vector(obj.glmm.null$LOCOResult[[chrom_v3]]$fitted.values))
        }
    } else {
        obj.model = list(obj.noK = obj.glmm.null$obj.noK, mu = as.vector(obj.glmm.null$fitted.values))
    }
    obj.model$offset = obj.glmm.null$offset
   
    if (traitType == "binary") {
        obj.model$mu2 = (obj.model$mu) * (1 - obj.model$mu)
    } else if (traitType == "quantitative") {
        obj.model$mu2 = (1/tauVec[1]) * rep(1, N)
    }
    IsOutputlogPforSingle= FALSE

    
    # check Sparse sigma
    if (sparseSigmaFile == "") {
        sparseSigma = NULL
        cat("sparse kinship matrix is not used\n")
    } else {
        cat("sparse kinship matrix is going to be used\n")
        Check_File_Exist(sparseSigmaFile, "sparseSigmaFile")
        
        sparseSigma = Matrix:::readMM(sparseSigmaFile)
        cat("sparseSigmaFile: ", sparseSigmaFile, "\n")
        
    }
    
     # check dropmissing doage 
    if (IsDropMissingDosages) {
        cat("Samples with missing dosages will be dropped from the analysis\n")
    } else {
        cat("Missing dosages will be mean imputed for the analysis\n")
    }
    
    SAIGE:::setIsDropMissingDosages_bgen(IsDropMissingDosages)
    SAIGE:::setIsDropMissingDosages_vcf(IsDropMissingDosages)
    
    startTime = as.numeric(Sys.time())
    cat("Analysis started at ", startTime, "Seconds\n")
    if (minMAC == 0) {
        minMAC = 0.5
        cat("As minMAC is set to be 0, minMAC = 0.5 will be used\n")
    }
    cat("minMAC: ", minMAC, "\n")
    cat("minMAF: ", minMAF, "\n")
    minMAFBasedOnMAC = minMAC/(2 * N)
    testMinMAF = max(minMAFBasedOnMAC, minMAF)
    cat("Minimum MAF of markers to be tested is ", testMinMAF,"\n")

	
	  adjustCCratioinGroupTest = TRUE
	  obj_cc =NULL
    if (traitType == "binary") {
        cat("It is a binary trait\n")
        
        if (SPAcutoff < 10^-2) {
            Cutoff = 10^-2
        }
        else {
            Cutoff = SPAcutoff
        }
        y1Index = which(y == 1)
        NCase = length(y1Index)
        y0Index = which(y == 0)
        NCtrl = length(y0Index)
        cat("Analyzing ", NCase, " cases and ", NCtrl, " controls \n")
        
        #obj_cc
        obj_cc<-SKAT:::SKAT_Null_Model(y~X-1, out_type="D", Adjustment=FALSE)
        obj_cc$mu = obj.model$mu
        obj_cc$res= y-obj_cc$mu
        obj_cc$pi_1 = obj_cc$mu*(1-obj_cc$mu)
       
    } else if (traitType == "quantitative") {
        cat("It is a quantitative trait\n")
        adjustCCratioinGroupTest = FALSE
    } else {
        stop("ERROR! The type of the trait has to be either binary or quantitative\n")
    }
    
    
    if (length(ratioVec) == 1) {
        cateVarRatioMinMACVecExclude = c(0, 2 * N)
        cateVarRatioMaxMACVecInclude = c(2 * N)
    }
    
    OUT_cond = NULL
    G2tilde_P_G2tilde_inv = NULL
    
    #
    # Group test
    #
    
    if (dosageFileType == "bgen") {
        SAIGE:::SetSampleIdx(sampleIndex, N)
    } else if (dosageFileType == "vcf") {
        SAIGE:::setMAFcutoffs(testMinMAF, maxMAFforGroupTest)
        cat("genetic variants with ", testMinMAF, "<= MAF <= ",  maxMAFforGroupTest, "are included for gene-based tests\n")
        SAIGE:::SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N)
    }
        
    OUT = NULL
    
    if (method_to_CollapseUltraRare != "") {
        if (MACCutoff_to_CollapseUltraRare <= 0) {
            stop("MACCutoff_to_CollapseUltraRare needs to be larger than 0\n")
        }
        if (DosageCutoff_for_UltraRarePresence <= 0 | DosageCutoff_for_UltraRarePresence > 2) {
            stop("DosageCutoff_for_UltraRarePresenc needs be to larger than 0 and less or equal to 2\n")
        }
        
        if (method_to_CollapseUltraRare == "absence_or_presence") {
                  cat("Ultra rare variants with MAC <= ", MACCutoff_to_CollapseUltraRare, 
                    " will be collpased for set-based tests in the 'absence or presence' way. ", 
                    "For the resulted collpased marker, any individual having ", 
                    DosageCutoff_for_UltraRarePresence, "<= dosage < ", 
                    (1 + DosageCutoff_for_UltraRarePresence), 
                    " for any ultra rare variant has 1 in the genotype vector, having dosage >= ", 
                    (1 + DosageCutoff_for_UltraRarePresence), 
                    " for any ultra rare variant has 2 in the genotype vector, otherwise 0. \n")
        }else if (method_to_CollapseUltraRare == "sum_geno") {
                  cat("Ultra rare variants with MAC <= ", MACCutoff_to_CollapseUltraRare, 
                    " will be collpased for set-based tests in the 'sum_geno' way. ", 
                    "The resulted collpased marker equals weighted sum of the genotypes of all ultra rare variantsi. NOTE: this option currently is not active\n")
        }
    }else {
        cat("Ultra rare variants won't be collpased for set-based association tests\n")
    }
    
    # read group file
    group_info_list =  ReadGroupFile(groupFile)
    
    ###########################
    # Main Test
    gene_base_test_df_A<-NULL
    single_test_df_A<-NULL
    Append1=FALSE
    ngroup<-length(group_info_list)
    for(i in 1:ngroup){
      geneID = group_info_list[[i]]$geneID
      markerIDs = group_info_list[[i]]$var
      marker_group_line<-paste(c(geneID, markerIDs), collapse = "\t")
      if (dosageFileType == "vcf") {
          Gx = getGenoOfGene_vcf(marker_group_line, minInfo)
      } else if (dosageFileType == "bgen") {
          Gx = getGenoOfGene_bgen_Sparse(bgenFile, bgenFileIndex, marker_group_line, testMinMAF, maxMAFforGroupTest, minInfo)
      }
      
      cntMarker = Gx$cnt
      cat("cntMarker: ", cntMarker, "\n")
      cat("marker_group_line", marker_group_line, "\n")
      cat("dosageFileType", dosageFileType, "\n")
      if (cntMarker == 0){
        next
      } 
      
      markerIDs = Gx$markerIDs
      Gmat = Matrix:::sparseMatrix(i = as.vector(Gx$iIndex), j = as.vector(Gx$jIndex), x = as.vector(Gx$dosages), symmetric = FALSE, dims = c(N, cntMarker))
      Gx$iIndex = NULL
      Gx$jIndex = NULL
      Gx$dosages = NULL
                        
      indexforMissing = unique(Gx$indexforMissing)
      N_sub = N  
      if (dosageZerodCutoff > 0 & sum(Gx$MACs <= 10) > 0) {
          zerodIndex = which(Gx$MACs <= 10)
          for (z in zerodIndex) {
            if (Gx$markerAFs[z] <= 0.5) {
              replaceindex = which(Gmat[, z] <= dosageZerodCutoff & Gmat[, z] > 0)
              if (length(replaceindex) > 0) {
                Gmat[replaceindex, z] = 0
              } 
            } else {
              replaceindex = which(Gmat[, z] >= (2 - dosageZerodCutoff) & Gmat[, z] < 2)
                if (length(replaceindex) > 0) {
                  Gmat[replaceindex, z] = 2
                }
              }
            }
      }
  
        
    	# Note: we don't remove marker here...
    	groupTestResult = MultiSets_GroupTest(Gmat=Gmat, obj.model=obj.model, obj_cc=obj_cc, 
                                                y=y, X=X, tauVec=tauVec, traitType=traitType, 
                                                function_group_marker_list=group_info_list[[i]], MAF_cutoff=MAF_cutoff, 
                                                function_group_test=function_group_test,
                                        cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, 
                                        cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude, 
                                        ratioVec=ratioVec, 
                                        kernel=kernel, method=method, 
                                        weights.beta.rare=weights.beta.rare,  weights.beta.common=weights.beta.common,weightMAFcutoff=weightMAFcutoff,
                                        r.corr=r.corr, sparseSigma=sparseSigma, IsSingleVarinGroupTest=IsSingleVarinGroupTest, 
                                        markerIDs=markerIDs, IsSparse=IsSparse, weights_specified=weights_specified,
                                        method_to_CollapseUltraRare = method_to_CollapseUltraRare,
                                        MACCutoff_to_CollapseUltraRare = MACCutoff_to_CollapseUltraRare, 
                                        DosageCutoff_for_UltraRarePresence = DosageCutoff_for_UltraRarePresence)
          
    	
    	out_df=Get_Results_DF(groupTestResult, geneID)
    	gene_base_test_df=out_df$gene_base_test_df
    	gene_base_test_df_A= rbind(gene_base_test_df_A, gene_base_test_df)
    	 
      	OUT_Filename<-SAIGEOutputFile
      	write.table(gene_base_test_df, file = OUT_Filename, col.names=!Append1, row.names=FALSE, quote=FALSE, append=Append1)
      	
      	if(IsSingleVarinGroupTest){
      	
      		single_test_df=out_df$single_test_df
      		single_test_df_A= rbind(single_test_df_A, single_test_df)   
      		OUT_Filename_Single<-sprintf("%s.single",SAIGEOutputFile )
      		if(!is.null(single_test_df_A)){
       			write.table(single_test_df, file = OUT_Filename_Single, col.names=!Append1, row.names=FALSE, quote=FALSE, append=Append1)
      		}
      	}
      	
      	Append1=TRUE
	}
	if (dosageFileType == "bgen") {
        closetestGenoFile_bgenDosage()
    }
    else if (dosageFileType == "vcf") {
        closetestGenoFile_vcfDosage()
    }
    summary(warnings())
    
    re_out = list(gene_base_test_df_A=gene_base_test_df_A, single_test_df_A=single_test_df_A)
    return(re_out)
}
