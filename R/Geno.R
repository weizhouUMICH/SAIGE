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

checkGenoInput = function(bgenFile = "",
                 bgenFileIndex = "",
                 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
                 savFileIndex = "",
                 sampleFile = "",
                 bedFile="",
                 bimFile="",
                 famFile=""){
   
    # genotype data
    if (bgenFile != "") {
        Check_File_Exist(bgenFile, "bgenFile")
        Check_File_Exist(bgenFileIndex, "bgenFileIndex")
        dosageFileType = "bgen"
    }else if (vcfFile != "") {
        Check_File_Exist(vcfFile, "vcfFile")

        if (!grepl("\\.sav$", vcfFile) && !file.exists(paste(vcfFile, ".csi", sep = ""))) {
            stop("ERROR! vcfFileIndex ", paste(vcfFile, ".csi", sep = ""), " does not exist\n")
        }
        dosageFileType = "vcf"
    }else if (savFile != "") {
        Check_File_Exist(savFile, "savFile")
	Check_File_Exist(savFileIndex, "savFileIndex")
	vcfFile = savFile
        dosageFileType = "vcf"
    }else if(bedFile != ""){
	Check_File_Exist(bedFile, "bedFile")
	if(bimFile == ""){
		bimFile = gsub("bed$", "bim", bedFile)
	}
	if(famFile == ""){
		famFile = gsub("bed$", "fam", bedFile)
	}	

	Check_File_Exist(bimFile, "bimFile")
	Check_File_Exist(famFile, "famFile")
	dosageFileType = "plink"
    }	    
    
    cat("dosageFile type is ", dosageFileType, "\n")

    if(dosageFileType == "vcf"){
        if (chrom == "") {
            stop("ERROR! chrom needs to be specified for the vcf file\n")
        }
	if(vcfField != "DS" & vcfField != "GT"){
	    stop("vcfField has to be DS or GT\n")	
	}
    }

    return(dosageFileType)
}	

setGenoInput = function(bgenFile = "",
                 bgenFileIndex = "",
                 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
                 savFileIndex = "",
                 sampleFile = "",
                 bedFile="",
                 bimFile="",
                 famFile="",
		 idstoExcludeFile = "",
                 idstoIncludeFile = "",
                 rangestoExcludeFile = "",
                 rangestoIncludeFile = "",
                 chrom = "",
                 start = 1,
                 end = 250000000, 
		 AlleleOrder = NULL)
{

  dosageFileType = checkGenoInput(bgenFile = bgenFile,
                 bgenFileIndex = bgenFileIndex,
                 vcfFile = vcfFile,
                 vcfFileIndex = vcfFileIndex,
                 vcfField = vcfField,
                 savFile = savFile,
                 savFileIndex = savFileIndex,
                 bedFile = bedFile,
                 bimFile = bimFile,
                 famFile = famFile)


  ########## ----------  Plink format ---------- ##########
  
  if(dosageFileType == "plink"){
    if(is.null(AlleleOrder)) AlleleOrder = "alt-first"

    cat("allle order in the plink file is ", AlleleOrder, ".\n")

    if(bimFile == ""){
    	bimFile = gsub("bed$", "bim", bedFile)
    }
    if(famFile == ""){
    	famFile = gsub("bed$", "fam", bedFile)
    } 
    markerInfo = data.table::fread(bimFile, header = F)
    markerInfo = as.data.frame(markerInfo)
    
    if(AlleleOrder == "alt-first")
      markerInfo = markerInfo[,c(1,4,2,6,5)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    if(AlleleOrder == "ref-first")
      markerInfo = markerInfo[,c(1,4,2,5,6)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex = 1:nrow(markerInfo) - 1  # -1 is to convert 'R' to 'C++' 
    
    sampleInfo = data.table::fread(famFile)
    samplesInGeno = sampleInfo$V2
    #SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
      
    setPLINKobjInCPP(bimFile, famFile, bedFile, AlleleOrder)
  }
  
  ########## ----------  BGEN format ---------- ##########
  
  if(dosageFileType == "bgen"){
    
    
    if(is.null(AlleleOrder)) AlleleOrder = "ref-first"
	
    if(!checkIfSampleIDsExist(bgenFile)){
	print("Sample IDs were not found in the bgen file.")
	Check_File_Exist(sampleFile)
	sampleData = data.table::fread(sampleFile, header=F, colClasses = c("character"))
        samplesInGeno = as.character(sampleData[,1])
    }else{
	samplesInGeno = getSampleIDsFromBGEN(bgenFile)
		    
    }    
    
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgenFileIndex)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    bgiData = dplyr::tbl(db_con, "Variant")
    bgiData = as.data.frame(bgiData)
    
    if(AlleleOrder == "alt-first")
      markerInfo = bgiData[,c(1,2,3,6,5,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    if(AlleleOrder == "ref-first")
      markerInfo = bgiData[,c(1,2,3,5,6,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT","genoIndex")
    
    setBGENobjInCPP(bgenFile, bgenFileIndex, samplesInGeno, AlleleOrder)
  }
  

  if(dosageFileType == "vcf"){
    setVCFobjInCPP(vcfFile, vcfFileIndex, vcfField, chrom, start, end)
    markerInfo = NULL
  }


  Files = c("IDsToIncludeFile", "IDsToExcludeFile", "RangesToIncludeFile", "RangesToExcludeFile")
 
  if(dosageFileType == "vcf"){
    if(!is.null(IDsToIncludeFile) | !is.null(IDsToExcludeFile) | !is.null(RangesToExcludeFile) | !is.null(RangesToIncludeFile)){
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

  genoList = list(dosageFileType = dosageFileType, markerInfo = markerInfo, anyQueue = anyQueue)
  return(genoList)
}




updateSampleIDs = function(SampleIDs, samplesInGeno)
{
  if(is.null(SampleIDs)){
    print("Since 'SampleIDs' not specified, we use all samples in 'GenoFile'.")
    SampleIDs = samplesInGeno;
  }
  
  if(any(!SampleIDs %in% samplesInGeno))
    stop("At least one sample from 'SampleIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  
  return(SampleIDs)
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



