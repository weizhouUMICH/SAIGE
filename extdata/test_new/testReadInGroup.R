> setwd("~/projects/Dec2021/SAIGE/extdata/test_new")
> library(SAIGE, lib.loc="~/projects/Dec2021/install")
> RegionList = SAIGE.getRegionList(groupFile, annolist, markerInfo)
> groupFile="./group_plink.txt"
> annolist = c("lof","missense")
> markerInfo = read.table("markerInfo.txt", header=F)
> groupFile="./group_plink.txt"
> RegionList = SAIGE.getRegionList(groupFile, annolist, markerInfo)
Start extracting marker-level information from 'groupFile' of ./group_plink.txt ....
Error in `[<-`(`*tmp*`, which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]),  :


		source("/net/dumbo/home/zhowei/projects/Dec2021/SAIGE/R/SAIGE_GENE_MultiVariantSet_Group.R")
	       source("/net/dumbo/home/zhowei/projects/Dec2021/SAIGE/R/Util.R")
	       source("/net/dumbo/home/zhowei/projects/Dec2021/SAIGE/R/SAIGE_SPATest_Region.R")
	       groupFile="./group_plink.txt"
	       annolist = c("lof","missense")
	       markerInfo = read.table("markerInfo.txt", header=F)
		RegionList = SAIGE.getRegionList(groupFile, annolist, markerInfo)
