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
                 famFile="", 
		 sampleInModel = NULL){
   
    if(is.null(sampleInModel)){
    	stop("sampleInModel is not specified.")
    }	    
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
        #if (chrom == "") {
        #    stop("ERROR! chrom needs to be specified for the vcf file\n")
        #}
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
                 idstoIncludeFile = "",
                 rangestoIncludeFile = "",
                 chrom = "",
		 AlleleOrder = NULL,
		 sampleInModel = NULL)
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
                 famFile = famFile,
		 sampleInModel = sampleInModel)

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
    markerInfo$ID = paste0(markerInfo$CHROM,":", markerInfo$POS ,"_", markerInfo$REF, "/", markerInfo$ALT) 
    setPLINKobjInCPP(bimFile, famFile, bedFile, sampleInModel, AlleleOrder)
  }
  
  ########## ----------  BGEN format ---------- ##########
  
  if(dosageFileType == "bgen"){
    
    
    if(is.null(AlleleOrder)) AlleleOrder = "ref-first"
    
    if(sampleFile != "" | !checkIfSampleIDsExist(bgenFile)){
	print("Sample IDs were not found in the bgen file.")
	Check_File_Exist(sampleFile)
	sampleData = data.table::fread(sampleFile, header=F, colClasses = c("character"))
        samplesInGeno = as.character(sampleData[,1])
    }else{
	samplesInGeno = getSampleIDsFromBGEN(bgenFile)
 	print(samplesInGeno[1:100])		    
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
    markerInfo$ID = paste0(markerInfo$CHROM,":", markerInfo$POS ,"_", markerInfo$REF, "/", markerInfo$ALT)    
    setBGENobjInCPP(bgenFile, bgenFileIndex, t_SampleInBgen = samplesInGeno, t_SampleInModel = sampleInModel, AlleleOrder)
  }
  

  if(dosageFileType == "vcf"){
    if(idstoIncludeFile != "" & rangestoIncludeFile != ""){
      stop("We currently do not support both 'idstoIncludeFile' and 'rangestoIncludeFile' at the same time for vcf files\n")
    }
    if(chrom==""){
      stop("chrom needs to be specified for VCF/BCF/SAV input\n")
    }
  }


  #Files = c("idstoIncludeFile", "idstoExcludeFile", "rangestoIncludeFile", "rangestoExcludeFile")

  anyInclude = FALSE
  #anyExclude = FALSE

  markersInclude = c()
  #markersExclude = c()
  IDsToInclude = NULL
  RangesToInclude = NULL
  if(idstoIncludeFile != ""){
    IDsToInclude = data.table::fread(idstoIncludeFile, header = F, colClasses = c("character"), data.table=F)
    if(ncol(IDsToInclude) != 1)
      stop("'idstoIncludeFile' of ", idstoIncludeFile, " should only include one column.")
    IDsToInclude = IDsToInclude[,1]

    posRows = which(markerInfo$ID %in% IDsToInclude)
    if(length(posRows) != 0)
      markersInclude = c(markersInclude, markerInfo$ID[posRows])
    anyInclude = TRUE
  }

  if(rangestoIncludeFile != ""){
    RangesToInclude = data.table::fread(rangestoIncludeFile, header = F, colClasses = c("character", "numeric", "numeric"), data.table=F)
    if(ncol(RangesToInclude) != 3)
      stop("rangestoIncludeFile should only include three columns.")

    colnames(RangesToInclude) = c("CHROM", "START", "END")

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

if(FALSE){

  if(!is.null(idstoExcludeFile)){
    if(anyInclude)
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    IDsToExclude = data.table::fread(idstoExcludeFile, header = F, colClasses = c("character"))
    if(ncol(IDsToExclude) != 1)
      stop("idstoExcludeFile should only include one column.")
    IDsToExclude = IDsToExclude[,1]

    posRows = which(markerInfo$ID %in% IDsToExclude)
    if(length(posRows) != 0)
      markersExclude = c(markersExclude, markerInfo$ID[posRows])
    anyExclude = TRUE
  }

  if(!is.null(rangestoExcludeFile)){
    if(anyInclude)
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    RangesToExclude = data.table::fread(rangestoExcludeFile, header = F, colClasses = c("character", "numeric", "numeric"))
    if(ncol(RangesToExclude) != 3)
      stop("rangestoExcludeFile should only include three columns.")

    colnames(RangesToExclude) = c("CHROM", "START", "END")

    for(i in 1:nrow(RangesToExclude)){
      CHROM1 = RangesToExclude$CHROM[i]
      START = RangesToExclude$START[i]
      END = RangesToExclude$END[i]
      posRows = with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
      if(length(posRows) != 0)
        markersExclude = c(markersExclude, markerInfo$ID[posRows])
    }
    anyExclude = TRUE
  }
}

  markersInclude = unique(markersInclude)
#  markersExclude = unique(markersExclude)

  # return genotype
  #cat("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data.\n")

  if(anyInclude)
    markerInfo = subset(markerInfo, ID %in% markersInclude)

#  if(anyExclude)
#    markerInfo = subset(markerInfo, !ID %in% markersExclude)

#  anyQueue = anyInclude | anyExclude

  if(dosageFileType == "vcf"){
    setVCFobjInCPP(vcfFile, vcfFileIndex, vcfField, t_SampleInModel = sampleInModel)


    if(!is.null(IDsToInclude)){
      SNPlist = paste(c("set1", IDsToInclude), collapse = "\t")
      in_chrom="fake_chrom"
      in_beg_pd=1
      in_end_pd=200000000
      set_iterator_inVcf(SNPlist, in_chrom, in_beg_pd, in_end_pd)
    }

    if(!is.null(RangesToInclude)){
      if(length(CHROM1) > 1){
        stop("We do not support query with multiple regions for vcf file. Please only include one region in the ", rangestoIncludeFile, "\n")
      }else{
	inSNPlist=""
        in_chrom=CHROM1[1]
  	in_beg_pd=START[1]
	in_end_pd=END[1]	
        set_iterator_inVcf(inSNPlist, in_chrom, in_beg_pd, in_end_pd)
      }
    } 
    markerInfo = NULL
  }
  #genoList = list(genoType = genoType, markerInfo = markerInfo, SampleIDs = SampleIDs, AlleleOrder = AlleleOrder, GenoFile = GenoFile, GenoFileIndex = GenoFileIndex, anyQueue = anyQueue)
  #genoList = list(dosageFileType = dosageFileType, markerInfo = markerInfo, anyQueue = anyQueue, genoType = dosageFileType)
  genoList = list(dosageFileType = dosageFileType, markerInfo = markerInfo, genoType = dosageFileType)
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


extract_genoIndex_condition = function(condition, markerInfo, genoType){
   if(condition != ""){
       	condition_original = unlist(strsplit(condition, ","))
	#if(!is.null(weight_cond)){
	#	weight_original = unlist(strsplit(weight_cond, ","))
	#}
	conditionDat = data.frame(SNP = condition_original, condIndex = seq(1,length(condition_original)))
   	if(genoType != "vcf"){
		markerInfo_conditionDat = merge(conditionDat, markerInfo, by.x="SNP", by.y="ID", sort = F)
		markerInfo_conditionDat = markerInfo_conditionDat[order(markerInfo_conditionDat$condIndex), ]	
		#markerInfo_conditionDat = markerInfo_conditionDat[which()]
       		#posInd = which(markerInfo$ID %in% condition_original)
		#if(length(posInd) == length(condition_original)){
		if(nrow(markerInfo_conditionDat) == length(condition_original)){
			cond_genoIndex = markerInfo_conditionDat$genoIndex
       			#cond_genoIndex = genoIndex[posInd]
		}else{

			stop(length(condition_original)-length(posInd), " conditioning markers are not found in the geno file. Please Check.\n")	
		}	
       }else{
	        condition_group_line = paste(c("condition", condition_original), collapse = "\t")	
		set_iterator_inVcf(condition_group_line, "1", 1, 200000000)
      		cond_genoIndex = rep(-1, length(condition_original)) 
       }
   }else{
	stop("condition is empty!")
   }
   return(cond_genoIndex)
}	
