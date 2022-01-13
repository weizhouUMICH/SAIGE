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
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "SAIGE")
#' getSampleIDsFromBGEN(BGENFile)
#' @export
getSampleIDsFromBGEN = function(bgenFile)
{
  if(!checkIfSampleIDsExist(bgenFile))
    stop("The BGEN file does not include sample identifiers. Please refer to ?checkIfSampleIDsExist and ?SAIGE.ReadGeno for more details")
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
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "SAIGE")
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

  if(grepl("\\.sav$", GenoFile) | grepl("\\.vcf$", GenoFile) | grepl("\\.vcf.gz$", GenoFile)){

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


  }else if(grepl("\\.bgen$", GenoFile)){
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
      stop("For genotype input of BGEN format , 'GenoFileIndex' should be of length 1 or 2. Check 'Details' section in '?SAIGE.ReadGeno' for more information.")

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

  }else if(grepl("\\.bed$", GenoFile)){
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

checkOutputFile = function(OutputFile,
                           OutputFileIndex,
                           AnalysisType,      ## "Marker" or "Region"
                           nEachChunk)
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  # message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."

  ## an R list of output
  if(missing(OutputFile))
    stop("Argument of 'OutputFile' is required.")

  if(file.exists(OutputFile)){
    if(!file.exists(OutputFileIndex)){
      stop(paste0("'OutputFile' of '", OutputFile,"' has existed.
                  Please use another 'OutputFile' or remove the existing one."))
    }
    else{
      outIndexData = read.table(OutputFileIndex, header = F, sep="\t")

      if(outIndexData[1,1] != message1 | outIndexData[2,1] != message2 | outIndexData[3,1] != message3)
        stop(paste0("'OutputFileIndex' of '", OutputFileIndex, "' is not as expected.
                    Probably, it has been modified by user, which is not permitted.
                    Please remove the existing files of 'OutputFile' and 'OutputFileIndex'."))

      lastMessage = outIndexData[nrow(outIndexData), 1]
      if(lastMessage == message5){
        End = TRUE
        indexChunk = outIndexData[nrow(outIndexData)-1, 1];
        indexChunk = as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        cat("Based on 'OutputFile' and 'OutputFileIndex', the analysis has been completed for the toal", indexChunk, "chunks.\n")
      }else{
        End = FALSE;
        indexChunk = lastMessage;
        indexChunk = as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        cat("Based on 'OutputFile' and 'OutputFileIndex', we restart the analysis from the", indexChunk+1, "chunk.\n")
      }
    }
    Start = FALSE
  }else{
    Start = TRUE;
    End = FALSE;
    indexChunk = 0;
    if(AnalysisType == "Marker"){
      indexChunk = 1;
    }	    
  }

  returnList = list(Start = Start, End = End, indexChunk = indexChunk)
  return(returnList)
}



writeOutputFile = function(Output,
                           OutputFile,
                           OutputFileIndex,
                           AnalysisType,
                           nEachChunk,
                           indexChunk,
                           Start,   # TRUE or FALSE, to indicate is the 'Output' is the first one to save into 'OutputFile'
                           End)     # TRUE or FALSE, to indicate is the 'Output' is the last one to save into 'OutputFile'
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."

  n1 = length(Output)
  n2 = length(OutputFile)
  print("write Output 1")

  if(n1 != n2)
    stop("length(Output) != length(OutputFile)")

  if(n1 != 0){
    for(i in 1:n1){
      if(Start){
        write.table(Output[[i]], OutputFile[[i]], quote = F, sep = "\t", append = F, col.names = T, row.names = F)
      }else{
        write.table(Output[[i]], OutputFile[[i]], quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
    }
  }
  print("write Output 2")
  if(Start)
    write.table(c(message1, message2, message3), OutputFileIndex,
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)

  write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)

  if(End)
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
}
