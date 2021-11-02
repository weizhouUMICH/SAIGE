bgenFile = "";
bgenFileIndex = "";
vcfFile = "./input/genotype_10markers.vcf.gz";
vcfFileIndex = "./input/genotype_10markers.vcf.gz.tbi";
vcfField = "GT";
savFile = "";
savFileIndex = "";
sampleFile = "./input/samplelist.txt";
idstoExcludeFile = "";
idstoIncludeFile = "";
rangestoExcludeFile = "";
rangestoIncludeFile = "";
chrom = "1";
start = 1;
end = 250000000;
minMAC = 0.5;
minMAF = 0;
maxMAFforGroupTest = 0.01;
minInfo = 0;
GMMATmodelFile = "./output/example_quantitative.rda";
varianceRatioFile = "./output/example_quantitative_cate.varianceRatio.txt";
SPAcutoff=2;
SAIGEOutputFile = "./test_multiset/example_quantitative.SAIGE.gene.txt";
numLinesOutput = 1;
IsSparse=TRUE;
LOCO=FALSE;
sparseSigmaFile="./output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx";
groupFile="./test_multiset/group_multiSets.txt";
kernel="linear.weighted";
method="optimal.adj";
weights.beta.rare = c(1,25);
weightMAFcutoff = 0.01;
weightsIncludeinGroupFile=FALSE;
r.corr=0;
IsSingleVarinGroupTest = TRUE;
cateVarRatioMinMACVecExclude=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5);
cateVarRatioMaxMACVecInclude=c(1.5,2.5,3.5,4.5,5.5,10.5,20.5);
dosageZerodCutoff = 0.2;
method_to_CollapseUltraRare="absence_or_presence";
MACCutoff_to_CollapseUltraRare = 10;
DosageCutoff_for_UltraRarePresence = 0.5;
function_group_test =c("lof", "missense");
MAF_cutoff=c(0.001,0.01);


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
ngroup<-length(group_info_list)
for(i in 1:ngroup){

  geneID = group_info_list[[i]][["all"]]$gene
  
  # Possible bug
  # Seems like there is a bug when read the same set twice
  # number of markers to read differ by the order
  
  markerIDs = sort(group_info_list[[i]][["all"]]$markerID)
  markerIDs = group_info_list[[i]][["all"]]$markerID   
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
  single_test_df = out_df$single_test_df
  
  gene_base_test_df_A= rbind(gene_base_test_df_A, gene_base_test_df)
  single_test_df_A= rbind(single_test_df_A, single_test_df)    
  closetestGenoFile_vcfDosage()
  
}
if (dosageFileType == "bgen") {
  closetestGenoFile_bgenDosage()
} else if (dosageFileType == "vcf") {
  closetestGenoFile_vcfDosage()
}
summary(warnings())

re_out = list(gene_base_test_df_A=gene_base_test_df_A, single_test_df_A=single_test_df_A)
