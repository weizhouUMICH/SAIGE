setwd("~/projects/Dec2021/SAIGE/extdata/test_new")
library(data.table)
data = fread("./group_plink.txt", header=F, data.table=F)
datanew = NULL
for(gene in unique(data[,1])){
	data2 = data[which(data[,1] == gene), c(3:ncol(data))]
	data3 = t(data2)
	for(anno in c("lof", "missense", "lof;missense")){
		genename = paste0(gene, "_", anno)
		annolist = unlist(strsplit(anno, split=";")[[1]])	
   		markerlist = data3[which(data3[,2] %in% annolist),1]
		datanew = rbind(datanew, c(genename, markerlist))
		write.table(datanew, "./group_old.txt", col.name=F, row.name=F, sep="\t", quote=F, append = T)
		datanew = NULL		
	}
}	
