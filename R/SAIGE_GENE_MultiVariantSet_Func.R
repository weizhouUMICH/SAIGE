


splitfun_weight = function(x) {
    return(strsplit(x, split = ";")[[1]][2])
}
splitfun_markerID = function(x) {
    return(strsplit(x, split = ";")[[1]][1])
}
    

Check_File_Exist<-function(file, filetype){
    if (!file.exists(file)) {
        stop("ERROR! ", filetype , file, " does not exsit\n")
    }
}

Check_Genotypes_and_Samples<-function(bgenFile, vcfFile, vcfField, vcfFileIndex, savFile, sampleFile, sampleID, chrom){

    # genotype data
    if (bgenFile != "") {
        Check_File_Exist(bgenFile, "bgenFile")
        dosageFileType = "bgen"
    }
    else if (vcfFile != "") {
        Check_File_Exist(vcfFile, "vcfFile")
        
        if (!grepl("\\.sav$", vcfFile) && !file.exists(paste(vcfFile, ".csi", sep = ""))) {
            stop("ERROR! vcfFileIndex ", paste(vcfFile, ".csi", sep = ""), " does not exist\n")
        }
        dosageFileType = "vcf"
        if (chrom == "") {
            stop("ERROR! chrom needs to be specified for the vcf file\n")
        }
    }
    else if (savFile != "") {
        Check_File_Exist(savFile, "savFile")
        vcfFile = savFile
        dosageFileType = "vcf"
    }
	
	  # Check sample file
    sampleListinDosage = NULL
    if (!file.exists(sampleFile)) {
        if (dosageFileType == "bgen") {
            stop("ERROR! The dosage file type is bgen but sampleFile ", sampleFile, " does not exsit\n")
        }
    }
    else {
        sampleListinDosage = data.frame(data.table:::fread(sampleFile, header = F, stringsAsFactors = FALSE, colClasses = c("character")))
        sampleListinDosage$IndexDose = seq(1, nrow(sampleListinDosage), by = 1)
        cat(nrow(sampleListinDosage), " sample IDs are found in sample file\n")
        colnames(sampleListinDosage)[1] = "IIDDose"
    }

	# SampleInModel 
    sampleInModel = NULL
    sampleInModel$IID = sampleID
    sampleInModel = data.frame(sampleInModel)
    sampleInModel$IndexInModel = seq(1, length(sampleInModel$IID), by = 1)
    cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")

	# When VCF
   	if (dosageFileType == "vcf") {
        vcffileopen = FALSE
        isVariant = SAIGE:::setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
        sampleListinDosage_vec = SAIGE:::getSampleIDlist_vcfMatrix()
        
        if (is.null(sampleListinDosage)) {
            sampleListinDosage = data.frame(IIDDose = sampleListinDosage_vec)
            sampleListinDosage$IndexDose = seq(1, nrow(sampleListinDosage), by = 1)
            cat(nrow(sampleListinDosage), " sample IDs are found in the vcf file\n")
        }
    }
    dataMerge = merge(sampleInModel, sampleListinDosage, by.x = "IID", by.y = "IIDDose")
    dataMerge_sort = dataMerge[with(dataMerge, order(IndexInModel)), ]

  	if (nrow(dataMerge_sort) < nrow(sampleInModel)) {
        stop("ERROR!", nrow(sampleInModel) - nrow(dataMerge_sort), " samples used in glmm model fit do not have dosages\n")
    } else {
        dataMerge_v2 = merge(dataMerge_sort, sampleListinDosage, by.x = "IID", by.y = "IIDDose", all.y = TRUE)
        dataMerge_v2_sort = dataMerge_v2[with(dataMerge_v2, order(IndexDose.y)), ]
        sampleIndex = dataMerge_v2_sort$IndexInModel
        N = sum(!is.na(sampleIndex))
        cat(N, " samples were used in fitting the NULL glmm model and are found in sample file\n")
        sampleIndex[is.na(sampleIndex)] = -10
        sampleIndex = sampleIndex - 1
        rm(dataMerge)
        rm(dataMerge_v2)
        rm(dataMerge_sort)
        rm(dataMerge_v2_sort)
    }
    
    rm(sampleInModel)


	re<-list(dosageFileType=dosageFileType, vcfFile=vcfFile, sampleIndex=sampleIndex)
	return(re)
	
}


Get_Variance_Ratio<-function(varianceRatioFile, sparseSigmaFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude){

	ratioVec = c(1)
	
    # check variance ratio
    if (!file.exists(varianceRatioFile)) {
        if (sparseSigmaFile == "") {
            stop("ERROR! varianceRatioFile ", varianceRatioFile, " does not exsit but sparseSigmaFile also does not exist \n")
        }
        else {
            cat("varianceRatioFile is not specified so variance ratio won't be used\n")
        }
        ratioVec = c(1)
    }
    else {
        varRatioData = data.frame(data.table:::fread(varianceRatioFile, header = F, stringsAsFactors = FALSE))
        ln = length(cateVarRatioMinMACVecExclude)
        hn = length(cateVarRatioMaxMACVecInclude)
        if (nrow(varRatioData) == 1) {
            stop("ERROR! To perform gene-based tests, categorical variance ratios are required\n")
        }
        else {
            ratioVec = varRatioData[, 1]
            nrv = length(ratioVec)
            if (nrv != ln) {
                stop("ERROR! The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude\n")
            }
            if (ln != (hn + 1)) {
                stop("ERROR! The length of cateVarRatioMaxMACVecInclude does not match with the lenght of cateVarRatioMinMACVecExclude (-1)\n")
            }
        }
        cat("variance Ratio is ", ratioVec, "\n")
    }
    
    return(ratioVec)
    
}


Get_Results_DF<-function(groupTestResult, geneID){
  
  if(is.null(groupTestResult)){
  	return(list(gene_base_test_df=NULL, single_test_df=NULL))
  }
  
  re_test_gene_base = groupTestResult$re_test_gene_base
  re_test_single = groupTestResult$re_test_single
  
  nSets<-length(re_test_gene_base)
  group.a<-rep("", nSets)
  cutoff.a<-rep("", nSets)
  outvecs<-NULL
  for(i in 1:nSets){
    group.a[i]<-re_test_gene_base[[i]]$group
    cutoff.a[i]<-re_test_gene_base[[i]]$cutoff
    outvecs = rbind(outvecs, re_test_gene_base[[i]]$outvec)
  }
  gene_base_test_df = data.frame(GeneID = geneID, group=group.a, cutoff=cutoff.a, outvecs=outvecs )
  
  #re$p.value,  re$m, re$BETA_Burden, re$SE_Burden, re$pval_Burden,  MACg, MAC_caseg, MAC_ctrlg)
  colnames(gene_base_test_df)<-c("GeneID", "FuncGroup", "Cutoff", "pval", "m",
                                 "pval_Burden", "MAC", "MAC_case", "MAC_ctrl")
  
  single_test_df = data.frame(GeneID=geneID, re_test_single$out_vecs_df )
  return(list(gene_base_test_df=gene_base_test_df, single_test_df=single_test_df))
}



