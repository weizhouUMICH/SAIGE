options(stringsAsFactors=F)
#' Run single variant score tests with SPA based on the logistic mixed model.
#'
#' @param dosageFile character. Path to dosage file. Each line contains dosages for a marker to be tested
#' @param dosageFileNrowSkip integer(>=0). Number of lines to be skiped in the dosage file.
#' @param dosageFileNcolSkip integer(>=0). Number of columns to be skiped in the dosage file
#' @param dosageFilecolnamesSkip vector of characters. The column names of the skipped columns. Default: c("SNPID", "CHR", "POS", "Allele0", "Allele1")
#' @param dosageFileChrCol string. The column name for the chromosome column. Must be in the dosageFilecolnamesSkip. Required If LOCO = TRUE and chrom ="".  
#' @param bgenFile character. Path to bgen file. Currently version 1.2 with 8 bit compression is supported
#' @param bgenFileIndex character. Path to the .bgi file (index of the bgen file)
#' @param vcfFile character. Path to vcf file
#' @param vcfFileIndex character. Path to index for vcf file by tabix, ".tbi" by "tabix -p vcf file.vcf.gz"
#' @param vcfField character. genotype field in vcf file to use. "DS" for dosages or "GT" for genotypes. By default, "DS".
#' @param savFile character. Path to sav file
#' @param savFileIndex character. Path to index for sav file .s1r
#' @param idstoExcludeFile character. Path to the file containing variant ids to be excluded from the bgen or vcf file
#' @param idstoIncludeFile character. Path to the file containing variant ids to be included from the bgen or vcf file
#' @param rangestoExcludeFile character. Path to the file containing genome regions to be excluded from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header 
#' @param rangestoIncludeFile character. Path to the file containing genome regions to be included from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header 
#' @param chrom character. string for the chromosome to include from vcf file. Required for vcf file. If LOCO is specified, providing chrom will save computation cost
#' @param start numeric. start genome position to include from vcf file. 
#' @param end numeric. end genome position to include from vcf file. 
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 1. The higher threshold between minMAC and minMAF will be used
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used
#' @param maxMAFforGroupTest numeric. Maximum minor allele frequency of markers to test in group test. By default 1.
#' @param minInfo numeric. Minimum imputation info of markers to test (in bgen file)
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the dosage, vcf, sav, or bgen file with NO header
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param Cutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param IsSparse logical. Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for dichotomous traits. By default, TRUE 
#' @param numLinesOutput numeric. Output results for how many marker each time.    
#' @param SAIGEOutputFile character. Path to the output file containing the SPAGMMAT test results
#' @param IsOutputAFinCaseCtrl logical. Whether to output allele frequency in cases and controls. By default, FALSE
#' @param groupFile character. Path to the group file containing one column "GeneID", and one column for ID of the tested genetic markers delimited by comma. This file is space-delimited can only work with the bgen,sav,and vcf format 
#' @param condition. For conditional analysis. Genetic marker ids (chr:pos_ref/alt) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A. Note that currently conditional analysis is only for vcf/sav input.
#' @return SAIGEOutputFile
#' @export
SKATtest_usingSKATLib = function(dosageFile = "",
                 dosageFileNrowSkip = 0, 
                 dosageFileNcolSkip = 0,
                 dosageFilecolnamesSkip = c("SNPID", "CHR", "POS", "Allele0", "Allele1"),
		 dosageFileChrCol = "CHR",   ##for LOCO
		 bgenFile = "",
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
		 maxMAFforGroupTest = 1,
        	 minInfo = 0,
                 SKATmodelFile = "", 
                 SAIGEOutputFile = "",
		 numLinesOutput = 10000, 
		 IsOutputAFinCaseCtrl=FALSE,
		 groupFile="",
		 condition="",
		 kernel="linear.weighted",
                 method="optimal.adj",
                 weights.beta=c(1,25),
                 r.corr=0
){


  #check and read files
  #output file
  if(!file.exists(SAIGEOutputFile)){
    file.create(SAIGEOutputFile, showWarnings = TRUE)
  }

  #file for the glmm null model
  if(!file.exists(SKATmodelFile)){
    stop("ERROR! SKATmodelFile ", SKATmodelFile, " does not exsit\n")
  }else{
    load(SKATmodelFile) #out.obj     
    sampleInModel = NULL
    sampleInModel$IID = out.obj$sampleID
    sampleInModel = data.frame(sampleInModel)
    sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
    cat(nrow(sampleInModel), " samples have been used to fit the null model\n")
    traitType = out.obj$traitType 
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
    }
  }

  ##Needs to check the number of columns and the number of samples in sample file
  if(dosageFile != ""){
    if(!file.exists(dosageFile)){
      stop("ERROR! dosageFile ", dosageFile, " does not exsit\n")
    }else{
      if(dosageFileNrowSkip < 0 | dosageFileNcolSkip < 0){
        stop("ERROR! dosageFileNrowSkip or dosageFileNcolSkip can't be less than zero\n")
      }
    }
    dosageFileType = "plain"

  }else if(bgenFile != ""){ 
    if(!file.exists(bgenFile)){
      stop("ERROR! bgenFile ", bgenFile, " does not exsit\n")
    }
    if(!file.exists(bgenFileIndex)){
      stop("ERROR! bgenFileIndex ", bgenFileIndex, " does not exsit\n")
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
      stop("ERROR! savFile ", savFile, " does not exist\n")
    }else{
      vcfFile = savFile	
    }

    if(!file.exists(savFileIndex)){
      stop("ERROR! savFileIndex ", savFileIndex, " does not exist\n")
    }else{
      vcfFileIndex = savFileIndex
    }	
    dosageFileType = "vcf"
  }

  if(dosageFileType != "plain"){
    if(!file.exists(groupFile)){
      stop("ERROR! groupFile ", groupFile, " does not exist\n")
    }
  }


  if(condition != ""){
    isCondition = TRUE
  }else{
    isCondition = FALSE
  }


  if(isCondition){
    conditionlist = paste(c("condMarkers",unlist(strsplit(condition,","))),collapse="\t")
    cat("conditionlist is ", conditionlist, "\n")
    if(dosageFileType == "vcf"){
      setMAFcutoffs(0, 0.5)
      isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
      SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N)
      Gx_cond = getGenoOfGene_vcf(conditionlist, minInfo)
    }else if(dosageFileType == "bgen"){
      Gx_cond = getGenoOfGene_bgen(bgenFile,bgenFileIndex, conditionlist, testMinMAF, 0.5)
    }else{
      cat("WARNING: conditional analysis can only work for dosageFileType vcf, sav or bgen\n")
    }
    #print(Gx_cond)
    cat("conditioning on ", unlist(Gx_cond$markerIDs), "\n")
    #G0 = Gx_cond$dosages
     cntMarker = Gx_cond$cnt
     if(cntMarker > 0){
	if(dosageFileType == "vcf"){
	 dosage_cond = Matrix:::sparseMatrix(i = as.vector(Gx_cond$iIndex), j = as.vector(Gx_cond$jIndex), x = as.vector(Gx_cond$dosages), symmetric = FALSE, dims = c(N, Gx_cond$cnt))
	 dosage_cond = as.matrix(dosage_cond)
	}else if(dosageFileType == "bgen"){
          dosage_cond = matrix(Gx_cond$dosages, byrow=F, ncol = cntMarker)
       }else{
         stop("ERROR: conditional analysis can only work for dosageFileType vcf, sav or bgen\n")
      }
}
    print(dim(dosage_cond))
  }


  #determine minimum MAF for markers to be tested
  if(minMAC == 0){minMAC = 1} ##01-19-2018
  cat("minMAC: ",minMAC,"\n")
  cat("minMAF: ",minMAF,"\n")
  minMAFBasedOnMAC = minMAC/(2*N) 
  testMinMAF = max(minMAFBasedOnMAC, minMAF) 
  cat("Minimum MAF of markers to be testd is ", testMinMAF, "\n")
  
  if(dosageFileType == "vcf"){ setMAFcutoffs(testMinMAF, maxMAFforGroupTest) }

  ##############START TEST########################
  startTime = as.numeric(Sys.time())  # start time of the SPAGMMAT tests
  cat("Analysis started at ", startTime, "Seconds\n")

  if(file.exists(SAIGEOutputFile)){file.remove(SAIGEOutputFile)}

#  sampleIndex = sampleIndex - 1

if(dosageFileType == "plain"){
  isCondition = FALSE
}else if (dosageFileType == "bgen"){
  SetSampleIdx(sampleIndex, N)
}else if(dosageFileType == "vcf"){
  isVariant = setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
  SetSampleIdx_forGenetest_vcfDosage(sampleIndex, N) 
}
#cat("sampleIndex: ", sampleIndex, "\n")

if(traitType == "quantitative"){
    OUT = NULL
    cat("It is a quantitative trait\n")
    mth = 0

    resultHeader = c("Gene","Pvalue", "Is.converged", "markerNumber", "testedMarkerNumber","markerIDs", "markerAFs")

    if(method=="optimal.adj"){
        resultHeader = c("Gene","Pvalue", "Is.converged", "markerNumber", "testedMarkerNumber","markerIDs", "markerAFs", "Pvalue_Burden","Pvalue_SKAT")
    }	


    write(resultHeader,file = SAIGEOutputFile, ncolumns = length(resultHeader))
#OUT = rbind(OUT, c(geneID, skatTest$p.value, skatTest$param$n.marker, skatTest$param$n.marker.test, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))
    gf = file(groupFile, "r")
    while ( TRUE ) {
      marker_group_line = readLines(gf, n = 1)
      if(length(marker_group_line) == 0 ){
        break
      }else{	
        print(marker_group_line)
        geneID = strsplit(marker_group_line, split="\t")[[1]][1]
        if(dosageFileType == "vcf"){
          Gx = getGenoOfGene_vcf(marker_group_line, minInfo)
        }else if(dosageFileType == "bgen"){
          Gx = getGenoOfGene_bgen(bgenFile,bgenFileIndex,marker_group_line, testMinMAF, maxMAFforGroupTest, minInfo)          
        }
        cntMarker = Gx$cnt

        G0 = Gx$dosages
	print(G0)
#	cat("markerIDs: ", Gx$markerIDs, "\n")
#	cat("G0: ", G0, "\n")
        if(cntMarker > 0){
	if(dosageFileType == "bgen"){
		Gmat = matrix(G0, byrow=F, ncol = cntMarker)
		Gmat = as(Gmat, "sparseMatrix") 
	}else if(dosageFileType == "vcf"){
		Gmat = Matrix:::sparseMatrix(i = as.vector(Gx$iIndex), j = as.vector(Gx$jIndex), x = as.vector(Gx$dosages), symmetric = FALSE, dims = c(N, Gx$cnt))
		
	}else{
      		stop("ERROR: gene-based tests can only work for dosageFileType vcf, sav or bgen\n")
    	}

         #Gmat = matrix(G0, byrow=F, ncol = cntMarker)	  
	 cat("dim(Gmat): ", dim(Gmat), "\n")	
	 Gmat = as.matrix(Gmat)
	 #cat("Gmat[,1]: ", Gmat[,1], "\n")	
	 #cat("colSums(Gmat): ", colSums(Gmat), "\n")
         skatTest = SKAT:::SKAT(Gmat, out.obj, max_maf = 1, method=method, kernel = kernel, weights.beta = weights.beta, r.corr = r.corr, is_check_genotype=FALSE, is_dosage = TRUE)
	print(skatTest$param)
	if(!is.null(skatTest$param$Is_Converged)){
		Is_Converged = skatTest$param$Is_Converged
	}else{
		Is_Converged = NA
	}
	 
	outVec = c(geneID, skatTest$p.value, Is_Converged, skatTest$param$n.marker, skatTest$param$n.marker.test, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";"))
	if(method=="optimal.adj"){
                #rho = 1 for burden, 0 for skat
	    if(cntMarker > 1){	
                p.val.vec = skatTest$param$p.val.each
                if(!is.null(skatTest$param$rho)){
			rho.val.vec = skatTest$param$rho
                	outVec = c(outVec, p.val.vec[which(rho.val.vec == 1)], p.val.vec[which(rho.val.vec == 0)])
			print("optimal.adj")
		}else{
			outVec = c(outVec, skatTest$param$liu_pval, skatTest$param$liu_pval)
		}	
	    }else{
		outVec = c(outVec, skatTest$param$liu_pval, skatTest$param$liu_pval)
	    }	
	}
	OUT = rbind(OUT, outVec)

#	 OUT = rbind(OUT, c(geneID, skatTest$p.value, skatTest$param$Is_Converged, skatTest$param$n.marker, skatTest$param$n.marker.test, paste(Gx$markerIDs, collapse=";"), paste(Gx$markerAFs, collapse=";")))


          mth = mth + 1
          if(mth %% numLinesOutput == 0){
	    ptm <- proc.time()
            print(ptm)
            print(mth)
            OUT = as.data.frame(OUT)
            write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
            OUT = NULL
          }
        }    
      }#end of else for if(length(line) == 0 )
    } # end of while ( TRUE ) {
    if(!is.null(OUT)){
      OUT = as.data.frame(OUT)
      write.table(OUT, SAIGEOutputFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
      OUT = NULL
    }
  }else{
    stop("ERROR! The type of the trait has to be quantitative\n")
  }

  #close the dosage file after tests
  if(dosageFileType == "plain"){
    closetestGenoFile_plainDosage()  
  }else if (dosageFileType == "bgen"){
    closetestGenoFile_bgenDosage()
  }else if(dosageFileType == "vcf"){
    closetestGenoFile_vcfDosage()
  }


  endTime = as.numeric(Sys.time()) #end time of the SPAGMMAT tests
  cat("Analysis ended at ", endTime, "Seconds\n")
  tookTime = endTime - startTime
  cat("Analysis took ", tookTime, "Seconds\n")
  
}
