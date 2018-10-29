library(data.table)
data = fread("pheno_1000samples.txt_withdosages.txt", header=T)
data = data.frame(data)

set.seed(1)
n.fam = 100

kin100 = as.matrix(as.data.frame(read.table("/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/simulation_Sep16/ped5ColsForPedigree10_10fam.GRM", header=F, stringsAsFactors = FALSE)))
fam.kin = kin100[1:10,1:10]
N = n.fam*nrow(fam.kin)

kins<-diag(n.fam) %x% fam.kin
out.eigen<-eigen(fam.kin)
id<-which(out.eigen$values > mean(out.eigen$values)*10^-5)
factor<- t(out.eigen$vectors) * sqrt(out.eigen$values)
kin.chol<-diag(n.fam) %x% factor

rn = rnorm(N)
tau = 2
b<-t(kin.chol) %*% rn[1:N]
b = b * sqrt((tau))

Y= b + data$x1 + data$x2 + rnorm(N)
y_quantitative = qnorm((rank(Y, na.last="keep")-0.5)/sum(!is.na(Y)))

data1 = cbind(y_quantitative, data)
colnames(data1)[2] = "y_binary"
write.table(data1, "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt", col.names=T, row.names=F, quote=F)

