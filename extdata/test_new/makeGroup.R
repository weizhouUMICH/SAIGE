> outfile="~/projects/Dec2021/SAIGE/extdata/test_new/group_plink.txt"
> a=c("GENE1","var1", paste0("rs",seq(1,50)))
> a
 [1] "GENE1" "var1"  "rs1"   "rs2"   "rs3"   "rs4"   "rs5"   "rs6"   "rs7"
[10] "rs8"   "rs9"   "rs10"  "rs11"  "rs12"  "rs13"  "rs14"  "rs15"  "rs16"
[19] "rs17"  "rs18"  "rs19"  "rs20"  "rs21"  "rs22"  "rs23"  "rs24"  "rs25"
[28] "rs26"  "rs27"  "rs28"  "rs29"  "rs30"  "rs31"  "rs32"  "rs33"  "rs34"
[37] "rs35"  "rs36"  "rs37"  "rs38"  "rs39"  "rs40"  "rs41"  "rs42"  "rs43"
[46] "rs44"  "rs45"  "rs46"  "rs47"  "rs48"  "rs49"  "rs50"
> write.table(a, outfile, quote=F, col.names=F, row.names=F)
> write.table(t(a), outfile, quote=F, col.names=F, row.names=F)
> b=c("GENE1","anno", rep("lof",10),rep("missense",20),rep("lof",20))
> write.table(t(b), outfile, quote=F, col.names=F, row.names=F, append=T)
> a=c("GENE2","var", paste0("rs",seq(51,100)))
> write.table(t(a), outfile, quote=F, col.names=F, row.names=F,append=T)
> b=c("GENE2","anno",rep("missense",30),rep("lof",20))
> write.table(t(b), outfile, quote=F, col.names=F, row.names=F, append=T)
