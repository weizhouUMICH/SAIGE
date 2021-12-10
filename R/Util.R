checkArgsBool = function(arg0, arg0name){
  if(!arg0 %in% c(TRUE, FALSE))
    stop(arg0name, " should be TRUE or FALSE.")
}

checkArgsNumeric = function(arg0, arg0name, minVal, maxVal){
  if(!is.numeric(arg0) | arg0 < minVal | arg0 > maxVal)
    stop(arg0name, " should be a numeric value ranging from 0 to 0.5.")
}



Check_File_Exist<-function(file, filetype){
    if (!file.exists(file)) {
        stop("ERROR! ", filetype , file, " does not exsit\n")
    }
}

Check_OutputFile_Create<-function(file){
  if (file.exists(file)) {
    file.remove(file)
  }
  if (!file.exists(file)) {
    file.create(file, showWarnings = TRUE)
  }
}	


create_resultHeader<-function(traitType,
			      isGroupTest,
			      isCondition,
			      IsOutputAFinCaseCtrl,
			      IsOutputNinCaseCtrl,
			      IsOutputHetHomCountsinCaseCtrl,
			      adjustCCratioinGroupTest,
			      IsOutputBETASEinBurdenTest,
			      IsOutputPvalueNAinGroupTestforBinary,
			      cateVarRatioMinMACVecExclude,
			      IsOutputMAFinCaseCtrlinGroupTest){

  dosageFilecolnamesSkip=c("CHR", "POS", "SNPID",
                "Allele1", "Allele2", "AC_Allele2", "AF_Allele2",
                "imputationInfo")
    
  if (traitType == "binary") {
        cat("It is a binary trait\n")
        if (!isGroupTest) {
            if (!isCondition) {
                resultHeader = c(dosageFilecolnamesSkip, "N",
                  "BETA", "SE", "Tstat", "p.value", "p.value.NA",
                  "Is.SPA.converge", "varT", "varTstar")
            }else{
                resultHeader = c(dosageFilecolnamesSkip, "N",
                  "BETA", "SE", "Tstat", "p.value", "p.value.NA",
                  "Is.SPA.converge", "varT", "varTstar", "Tstat_cond",
                  "p.value_cond", "varT_cond", "BETA_cond", "SE_cond")
            }

            if (IsOutputAFinCaseCtrl) {
                resultHeader = c(resultHeader, "AF.Cases", "AF.Controls")
            }
            if (IsOutputNinCaseCtrl) {
                resultHeader = c(resultHeader, "N.Cases", "N.Controls")
            }
            if (IsOutputHetHomCountsinCaseCtrl) {
                resultHeader = c(resultHeader, "homN_Allele2_cases",
                  "hetN_Allele2_cases", "homN_Allele2_ctrls",
                  "hetN_Allele2_ctrls")
            }
            write(resultHeader, file = SAIGEOutputFile, ncolumns = length(resultHeader))
      } #if (!isGroupTest) {
    }else if (traitType == "quantitative") {
	adjustCCratioinGroupTest = FALSE    
        cat("It is a quantitative trait\n")
        if (!isGroupTest) {
            if (!isCondition) {
                resultHeader = c(dosageFilecolnamesSkip, "N",
                  "BETA", "SE", "Tstat", "p.value", "varT", "varTstar")
            }
            else {
                resultHeader = c(dosageFilecolnamesSkip, "N",
                  "BETA", "SE", "Tstat", "p.value", "varT", "varTstar",
                  "Tstat_cond", "p.value_cond", "varT_cond",
                  "BETA_cond", "SE_cond")
            }
            write(resultHeader, file = SAIGEOutputFile, ncolumns = length(resultHeader))
       }
    }else {
        stop("ERROR! The type of the trait has to be either binary or quantitative\n")
    }


    if(isGroupTest){
      MACcateNumHeader = paste0("Nmarker_MACCate_", seq(1,
            length(cateVarRatioMinMACVecExclude)))

    if (!isCondition) {
      if (adjustCCratioinGroupTest) {
                resultHeader = c("Gene", "Pvalue", MACcateNumHeader,
                  "markerIDs", "markerAFs")
                #if (method == "optimal.adj") {
                  if (IsOutputBETASEinBurdenTest) {
                    resultHeader = c("Gene", "Pvalue", MACcateNumHeader,
                      "markerIDs", "markerAFs", "Pvalue_Burden",
                      "Pvalue_SKAT", "BETA_Burden", "SE_Burden")
                  }else{
                    resultHeader = c("Gene", "Pvalue", MACcateNumHeader,
                      "markerIDs", "markerAFs", "Pvalue_Burden",
                      "Pvalue_SKAT")
                  }
                #}
    }

              if (IsOutputPvalueNAinGroupTestforBinary) {
                if (!adjustCCratioinGroupTest) {
                  resultHeader = c("Gene", "Pvalue", MACcateNumHeader,
                    "markerIDs", "markerAFs")
                  if (method == "optimal.adj") {
                    if (IsOutputBETASEinBurdenTest) {
                      resultHeader = c("Gene", "Pvalue", MACcateNumHeader,
                        "markerIDs", "markerAFs", "Pvalue_Burden",
                        "Pvalue_SKAT", "BETA_Burden", "SE_Burden")
                    }
                    else {
                      resultHeader = c("Gene", "Pvalue", MACcateNumHeader,
                        "markerIDs", "markerAFs", "Pvalue_Burden",
                        "Pvalue_SKAT")
                    }
                  }
                }else {
                  resultHeader = c(resultHeader, "Pvalue.NA")
                  if (method == "optimal.adj") {
                    if (IsOutputBETASEinBurdenTest) {
                      resultHeader = c(resultHeader, "Pvalue_Burden.NA",
                        "Pvalue_SKAT.NA", "BETA_Burden.NA", "SE_Burden.NA")
                    }
                    else {
                      resultHeader = c(resultHeader, "Pvalue_Burden.NA",
                        "Pvalue_SKAT.NA")
                    }
                  }
                }
            }    

   }else{
            if (adjustCCratioinGroupTest) {
                resultHeader = c("Gene", "Pvalue", "Pvalue_cond",
                  MACcateNumHeader, "markerIDs", "markerAFs")
                if (method == "optimal.adj") {
                  if (IsOutputBETASEinBurdenTest) {
                    resultHeader = c("Gene", "Pvalue", "Pvalue_cond",
                      MACcateNumHeader, "markerIDs", "markerAFs",
                      "Pvalue_Burden", "Pvalue_Burden_cond",
                      "Pvalue_SKAT", "Pvalue_SKAT_cond", "BETA_Burden",
                      "SE_Burden", "BETA_Burden_cond", "SE_Burden_cond")
                  }
                  else {
                    resultHeader = c("Gene", "Pvalue", "Pvalue_cond",
                      MACcateNumHeader, "markerIDs", "markerAFs",
                      "Pvalue_Burden", "Pvalue_Burden_cond",
                      "Pvalue_SKAT", "Pvalue_SKAT_cond")
                  }
                }
            }

	               if (IsOutputPvalueNAinGroupTestforBinary) {
                if (!adjustCCratioinGroupTest) {
                  resultHeader = c("Gene", "Pvalue", "Pvalue_cond",
                    MACcateNumHeader, "markerIDs", "markerAFs")
                  if (method == "optimal.adj") {
                    if (IsOutputBETASEinBurdenTest) {
                      resultHeader = c("Gene", "Pvalue", "Pvalue_cond",
                        MACcateNumHeader, "markerIDs", "markerAFs",
                        "Pvalue_Burden", "Pvalue_Burden_cond",
                        "Pvalue_SKAT", "Pvalue_SKAT_cond", "BETA_Burden",
                        "SE_Burden", "BETA_Burden_cond", "SE_Burden_cond")
                    }
                    else {
                      resultHeader = c("Gene", "Pvalue", "Pvalue_cond",
                        MACcateNumHeader, "markerIDs", "markerAFs",
                        "Pvalue_Burden", "Pvalue_Burden_cond",
                        "Pvalue_SKAT", "Pvalue_SKAT_cond")
                    }
                  }
                }else{
                  resultHeader = c(resultHeader, "Pvalue.NA",
                    "Pvalue.NA_cond")
                  #if (method == "optimal.adj") {
                    if (IsOutputBETASEinBurdenTest) {
                      resultHeader = c(resultHeader, "Pvalue_Burden.NA",
                        "Pvalue_Burden.NA_cond", "Pvalue_SKAT.NA",
                        "Pvalue_SKAT.NA_cond", "BETA_Burden.NA",
                        "SE_Burden.NA", "BETA_Burden.NA_cond",
                        "SE_Burden.NA_cond")
                    }
                    else {
                      resultHeader = c(resultHeader, "Pvalue_Burden.NA",
                        "Pvalue_Burden.NA_cond", "Pvalue_SKAT.NA",
                        "Pvalue_SKAT.NA_cond")
                    }
                 # }
                }
            }
    }

   if (IsOutputMAFinCaseCtrlinGroupTest) {
       resultHeader = c(resultHeader, "MAF_in_cases", "MAF_in_controls")
   }
   write(resultHeader, file = SAIGEOutputFile, ncolumns = length(resultHeader))


   }#if(isGroupTest)
}








# https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
#' Extract sample identifiers from BGEN file (for BGEN v1.2)
#'
#' Extract sample identifiers from BGEN file (for BGEN v1.2)
#'
#' @param bgenFile a character of BGEN file.
#' @examples
#'
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "GRAB")
#' getSampleIDsFromBGEN(BGENFile)
#' @export
getSampleIDsFromBGEN = function(bgenFile)
{
  if(!checkIfSampleIDsExist(bgenFile))
    stop("The BGEN file does not include sample identifiers. Please refer to ?checkIfSampleIDsExist and ?GRAB.ReadGeno for more details")
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH + 4)
  N = readBin(con, n = 1, what = "integer", size = 4)  # number of samples
  samplesInGeno = rep(0, N)

  # cycle for all samples to extract IDs
  for(i in 1:N){
    LS = readBin(con, n = 1, what = "integer", size = 2)
    sample = readChar(con, nchars = LS)
    samplesInGeno[i] = sample
  }

  # close connection
  close(con)

  return(samplesInGeno)
}

#' Check if sample identifiers are stored in a BGEN file (for BGEN v1.2)
#'
#' Check if sample identifiers are stored in a BGEN file (for BGEN v1.2)
#'
#' @param bgenFile a character of BGEN file. Sometimes, BGEN file does not include sample IDs. This information can be extracted from BGEN file. Please refer to https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for more details.
#' @examples
#'
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "GRAB")
#' checkIfSampleIDsExist(BGENFile)
#' @export
checkIfSampleIDsExist = function(bgenFile)
{
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH - 4)
  header = rawToBits(readBin(con, n = 4, what = "raw", size = 1, signed = FALSE))
  close(con)
  return(header[32] == 01)
}

Check_GenoFile<-function(GenoFile,
                         GenoFileIndex = NULL,
			 control = NULL){
  if(missing(GenoFile))
    stop("Argument 'GenoFile' is required.")

  Check_File_Exist(GenoFile, "GenoFile")
  #Check_File_Exist(GenoFileIndex, "GenoFileIndex")

  checkControl.ReadGeno(control) 
  AlleleOrder<-control$AlleleOrder
  chrom<-control$chrom
  vcfField<-control$vcfField
  sampleFile<-control$sampleFile

  if(grepl("\\.sav$", GenoFile) | grepl("\\.vcf$", GenoFile) | grepl("\\.vcf.gz$", GenoFile){

    dosageFileType <- "VCF"
    cat("genotypes/dosages file type is VCF\n")
    if (is.null(chrom)) {
      stop("ERROR! chrom needs to be specified for the VCF file\n")
    }

    if(is.null(GenoFileIndex)){
      GenoFileIndex = paste0(GenoFile, ".csi") 
      cat("'GenoFileIndex' is not given,", GenoFileIndex," will be used for the index file\n")
    }

    Check_File_Exist(GenoFileIndex, "GenoFileIndex")

    if(is.null(vcfField)) {
      vcfField <- "DS"
      cat("vcfField is not specified, so DS will be read as vcfField from the vcf file\n")
    }else{
      cat(vcfField, " will be read as vcfField from the vcf file\n")
    }
    samplesInGeno = getSampleIDsFromVCF(GenoFile, GenoFileIndex, vcfField)

    GenoFileList = list(vcfFile = GenoFile, vcfFileIndex = GenoFileIndex, vcfField = vcfField, chrom = chrom)


  }else if(grepl("\\.bgen$", GenoFile){
    dosageFileType <- "BGEN"
    cat("genotypes/dosages file type is BGEN\n")

    if(is.null(AlleleOrder)) {
      AlleleOrder <- "ref-first"
      cat("AlleleOrder is not specified, so bgen will be read as ref-first\n")
    }else{
      cat("bgen will be read as ",AlleleOrder,"\n")
    }	    

    if(is.null(GenoFileIndex)){
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bgen.bgi' file
      GenoFileIndex = c(gsub("bgen$", "bgen.bgi", GenoFile),
                        gsub("bgen$", "bgen.samples", GenoFile))
    }

    if(length(GenoFileIndex) != 1 & length(GenoFileIndex) != 2)
      stop("For genotype input of BGEN format , 'GenoFileIndex' should be of length 1 or 2. Check 'Details' section in '?GRAB.ReadGeno' for more information.")

    Check_File_Exist(GenoFileIndex[1], "bgiFile")

    if(length(GenoFileIndex) == 1){
      cat("No sample file for bgen is specified. Sample IDs will be read from the bgen file\n")
      samplesInGeno = getSampleIDsFromBGEN(bgenFile)
    }

    if(length(GenoFileIndex) == 2){
      sampleFile = GenoFileIndex[2]
      cat("Sample file for bgen is specified. Sample IDs will be read from the sample file\n")
      Check_File_Exist(sampleFile, "sampleFile")
      sampleData = data.table::fread(sampleFile, header=T, colClasses = c("character"))
      if(toupper(colnames(sampleData)[1]) != "BGEN_SAMPLE"){	      
        stop("The header of the first column in bgen.samples file should be 'BGEN_SAMPLE'.")
        samplesInGeno = as.character(sampleData[,1])
      }
    }

    GenoFileList = list(bgenFile = GenoFile, bgenFileIndex = GenoFileIndex[1], AlleleOrder = AlleleOrder, chrom=chrom)

  }else if(grepl("\\.bed$", GenoFile){
        if(is.null(AlleleOrder)) AlleleOrder = "alt-first"

    dosageFileType <- "PLINK"
    cat("genotypes/dosages file type is PLINK\n")	
    if(is.null(GenoFileIndex)){
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bim' and 'fam' files
      GenoFileIndex = c(gsub("bed$", "bim", GenoFile),
                        gsub("bed$", "fam", GenoFile))
    }

    if(length(GenoFileIndex) != 2)
      stop("If Plink format is used, argument 'GenoFileIndex' should be 'NULL' or a character vector of c(bimFile, famFile).")

    bimFile = GenoFileIndex[1]
    famFile = GenoFileIndex[2]
    bedFile = GenoFile
    Check_File_Exist(bedFile, "bedFile")
    Check_File_Exist(bimFile, "bimFile")
    Check_File_Exist(famFile, "famFile")



    if(is.null(AlleleOrder)) {
      AlleleOrder <- "ref-first"
      cat("AlleleOrder is not specified, so bgen will be read as ref-first\n")
    }else{
      cat("bed file will be read as ",AlleleOrder,"\n")
    }

    sampleInfo = data.table::fread(famFile)
    samplesInGeno = sampleInfo$V2
    GenoFileList = list(bedFile = bedFile, bimFile = bimFile, famFile = famFile, AlleleOrder = AlleleOrder, chrom=chrom)
  }
  cat(length(samplesInGeno)," samples are found in the dosages/genotypes file\n")
  return(samplesInGeno = samplesInGeno, dosageFileType = dosageFileType, GenoFileList = GenoFileList)  
}

  #genoFilelist <- Check_GenoFile(GenoFile = GenoFile, GenoFileIndex = GenoFileIndex, control = control) 
  #sampleIndex <- Check_GenotypeSamples(genoFilelist$samplesInGeno, sampleInModel) 
  #dosageFileType = genoFilelist$dosageFileType

set_Geno = function(dosageFileType,
                    genoFilelist,
		    sampleInModel,
                    control = NULL){

	
  if(dosageFileType == "VCF"){
      setVCFobjInCPP(genoFilelist$GenoFileList$vcfFile, genoFilelist$GenoFileList$vcfFileIndex, genoFilelist$GenoFileList$vcfField, genoFilelist$GenoFileList$chrom, , F, T, sampleInModel)
      markerInfo=NULL	
  }else if(dosageFileType == "BGEN"){
      setBGENobjInCPP(genoFilelist$GenoFileList$bgenFile, genoFilelist$GenoFileList$bgenFileIndex, genoFilelist$GenoFileList$samplesInGeno, sampleInModel, F, T, genoFilelist$GenoFileList$AlleleOrder)
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), genoFilelist$GenoFileList$bgenFileIndex)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    bgiData = dplyr::tbl(db_con, "Variant")
    bgiData = as.data.frame(bgiData)

    if(AlleleOrder == "alt-first")
      markerInfo = bgiData[,c(1,2,3,6,5,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    if(AlleleOrder == "ref-first")
      markerInfo = bgiData[,c(1,2,3,5,6,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html

    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT","genoIndex")	  
  }else if(dosageFileType == "PLINK"){
    setPLINKobjInCPP(genoFilelist$GenoFileList$bimFile, genoFilelist$GenoFileList$famFile, genoFilelist$GenoFileList$bedFile, sampleInModel, genoFilelist$GenoFileList$AlleleOrder)
    markerInfo = data.table::fread(genoFilelist$GenoFileList$bimFile, header = F)
    markerInfo = as.data.frame(markerInfo)

    if(AlleleOrder == "alt-first")
      markerInfo = markerInfo[,c(1,4,2,6,5)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    if(AlleleOrder == "ref-first")
      markerInfo = markerInfo[,c(1,4,2,5,6)]  # https://www.cog-genomics.org/plink/2.0/formats#bim

    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex = 1:nrow(markerInfo) - 1	
  }

  ##query
  Files = c("IDsToIncludeFile", "IDsToExcludeFile", "RangesToIncludeFile", "RangesToExcludeFile")


  if(dosageFileType == "VCF"){
    if(!is.null(control$IDsToIncludeFile) | !is.null(control$IDsToExcludeFile) | !is.null(control$RangesToExcludeFile) | !is.null(control$RangesToIncludeFile)){	  
      stop("We currenlty do not support the query function for VCF files. Please use chrom, start, end for query\n")
    }
  }

  anyInclude = FALSE
  anyExclude = FALSE

  markersInclude = c()
  markersExclude = c()

  if(!is.null(control$IDsToIncludeFile)){

    for(i in 1:nrow(RangesToInclude)){
      CHROM1 = RangesToInclude$CHROM[i]
      START = RangesToInclude$START[i]
      END = RangesToInclude$END[i]
      posRows = with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
      if(length(posRows) != 0)
        markersInclude = c(markersInclude, markerInfo$ID[posRows])
    }
    anyInclude = TRUE
  }

  if(!is.null(control$IDsToExcludeFile)){
    if(anyInclude)
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile'.")
    IDsToExclude = data.table::fread(control$IDsToExcludeFile,
                                     header = F, colClasses = c("character"))
    if(ncol(IDsToExclude) != 1)
      stop("IDsToExcludeFile should include one column.")
    IDsToExclude = IDsToExclude[,1]

    posRows = which(markerInfo$ID %in% IDsToExclude)
    if(length(posRows) != 0)
      markersExclude = c(markersExclude, markerInfo$ID[posRows])
    anyExclude = TRUE
  }

  if(!is.null(control$RangesToExcludeFile)){
    if(anyInclude)
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile'.")

    RangesToExclude = data.table::fread(control$RangesToExcludeFile,
                                        header = F, colClasses = c("character", "numeric", "numeric"))
    if(ncol(RangesToExclude) != 3)
      stop("RangesToExcludeFile should include three columns.")

    colnames(RangesToExclude) = c("CHROM","START","END")

    for(i in 1:nrow(RangesToExclude)){
      CHROM1 = RangesToExclude$CHROM[i]
      START = RangesToExclude$START[i]
      END = RangesToExclude$END[i]
      if(length(posRows) != 0)
        markersExclude = c(markersExclude, markerInfo$ID[posRows])
    }
    anyExclude = TRUE
  }

  markersInclude = unique(markersInclude)
  markersExclude = unique(markersExclude)

  # return genotype
  print(paste("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data."))

  if(anyInclude)
    markerInfo = subset(markerInfo, ID %in% markersInclude)

  if(anyExclude)
    markerInfo = subset(markerInfo, !ID %in% markersExclude)

  anyQueue = anyInclude | anyExclude

  genoList = list(markerInfo = markerInfo, anyQueue = anyQueue)
  return(genoList)
}

