Beta_Weight<-function(MAF,weights.beta){

	n<-length(MAF)
	weights<-rep(0,n)	
	IDX_0<-which(MAF == 0)
	if(length(IDX_0) == n){
		stop("No polymorphic SNPs")
	} else if( length(IDX_0) == 0){
		weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
	} else {
		weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
	}

	
	#print(length(IDX_0))
	#print(weights[-IDX_0])
	return(weights)
	
}



SPA_ER_kernel_related<-function(G,obj,  obj.noK, Cutoff=2, Phi,  weight,VarRatio_Vec, mu.a){

	zscore.all_0<-matrix(rep(0, ncol(G)), ncol=ncol(G))
	zscore.all_1<-matrix(rep(0, ncol(G)), ncol=ncol(G))
	VarS=c()
	
	g.sum=0
	q.sum=0
	p.old=c()
	p.new=c()
	MAFsum=colSums(G)
	for (jj in 1:ncol(G)){
		n.g<-sum(G[,jj])
		NAset<-which(G[,jj]==0)

		G1<-G[,jj]  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G[,jj]) ####equal to G[,jj] in terms of score statistics.
		q<-(t(G1) %*% (obj.noK$y)) 
		g=G1 


		mu.qtemp=mu.a; g.qtemp=g   
		mu1 <- sum(mu.qtemp * g.qtemp)
		var1<-Phi[jj, jj]/weight[jj]^2

 		stat.qtemp<-(q - mu1)^2/var1
    		p_temp1<-pchisq(stat.qtemp, lower.tail = FALSE, df = 1)  
		p.old[jj]=p_temp1

		zscore.all_0[,jj]=(q-mu1)  ##sum(G[,jj]*(obj.noK$y-mu.a))  
		
		id1<-which(stat.qtemp > Cutoff^2) 
#		print("test5")	
		if (MAFsum[jj]<=10){
			if (length(id1)>0 ){
				G_temp=G[,jj]
				G_temp[which(G_temp<=0.2)]=0
#		print("test6")	
#		print(dim(as.matrix(G_temp)))	
#		print(obj$res)
				p_temp1=SKAT::SKATBinary(as.matrix(G_temp),obj, method.bin="Hybrid")$p.value
#		print("test7")	
										
			}
	
		}else {
			if (length( id1)>0){  				
    				p_temp1 = scoreTest_SPAGMMAT_binaryTrait(g, n.g, NAset, obj.noK$y, mu.a, varRatio=VarRatio_Vec[jj], Cutoff = Cutoff)$p.value
#		print("test8")	
			}
		}	
		
		p.new[jj]=p_temp1
		if (Phi[jj,jj]<=0){zscore.all_1[,jj]=0} else{
			zscore.all_1[,jj]=qnorm(p_temp1/2, mean = 0, sd =sqrt( Phi[jj,jj]),lower.tail = FALSE, log.p = FALSE)*sign(q-mu1)
		}
		if (p_temp1>0){
			VarS[jj]= zscore.all_0[,jj]^2/qchisq(p_temp1, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
		} else {
			VarS[jj]= zscore.all_0[,jj]^2/500 
		}
		if (p_temp1<1){        
			g.sum = g.sum + g.qtemp * weight[jj] 
        		q.sum = q.sum + q * weight[jj] 
		}
	}##for every col of G
#		print("test9")	
	outlist=list();
	
	outlist$zscore.all_0=zscore.all_0
	outlist$VarS=VarS
	outlist$mu=mu.qtemp
	outlist$g.sum=g.sum
	outlist$q.sum=q.sum
	outlist$p.old=p.old
	outlist$p.new=p.new
	outlist$zscore.all_1=zscore.all_1
	return(outlist) ;
}

Related_ER<-function(G, obj, obj.noK, ratioVec=ratioVec,sparseSigma, mac_cutoff, Cutoff=2, weights.beta=c(1,25)){
	if (length(G)==0) {stop("WARNING: no-variantion in the whole genotype matrix!\n")}
    	for (gi in 1:dim(G)[2]){
		temp_gi=which(G[,gi]==9 | G[,gi]==NA)
		if (length(temp_gi)>0){
			G[temp_gi,gi]=mean(G[-temp_gi,gi])
			cat("The missing values in column", gi," are imputed by the mean genotype value. \n")
    		}
    	}
	mu.a=       obj.noK$mu
	mu2.a=    obj.noK$V
    	MAF_0 = which(colSums(G)==0)
    	if (length(MAF_0)>0){
		cat("The following columns are removed due to no-variation: ", MAF_0,"\n")
		G=G[,-MAF_0]
	} 
	if (length(G)==0) {stop("WARNING: no-variantion in the whole genotype matrix!\n")}

	mac_cutoff=c(0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5)
	for (jj in 1:ncol(G)){
		n.g<-sum(G[,jj])
		if(n.g/(2*length(G[,jj]))>0.5)
		{
			G[,jj]<-2-G[,jj]
			n.g<-sum(G[,jj])
		}

	}

	MAF=colMeans(G)/2
	MAFsum=colSums(G)
	VarRatio_Vec=rep(0,length(MAFsum))
	for (G_k in 1:length(MAFsum)){
		for (mac_k in 1:(length(mac_cutoff)-1)){
			if (MAFsum[G_k]>=mac_cutoff[mac_k] & MAFsum[G_k]<mac_cutoff[mac_k+1]){
				VarRatio_Vec[G_k]=ratioVec[mac_k]
			} 
		}
		if (MAFsum[G_k]>=mac_cutoff[length(mac_cutoff)]){VarRatio_Vec[G_k]=ratioVec[length(mac_cutoff)]}

	}


        	
	mafcutoff=0.01  ###########1/sqrt(nrow(G_o) * 2)
	weight=rep(0,length(MAF))
	maf_temp=which(MAF<mafcutoff)##rare variants
	if (length(maf_temp)>0 & length(maf_temp)<length(MAF)){
		#weight[maf_temp]=Beta_Weight(MAF[maf_temp],c(1,25))
		weight[maf_temp]=SKAT:::Beta.Weights(MAF[maf_temp],weights.beta)
		#weight[-maf_temp]=Beta_Weight(MAF[-maf_temp],c(0.5,0.5))
		weight[-maf_temp]=SKAT:::Beta.Weights(MAF[-maf_temp],c(0.5,0.5))
		flag=1

	} else {
		#if (length(maf_temp)==0){weight=Beta_Weight(MAF,c(0.5,0.5));flag=2} ###flag 1 means both; 2 only common; 3 only rare;
		if (length(maf_temp)==0){weight=SKAT:::Beta.Weights(MAF,c(0.5,0.5));flag=2} ###flag 1 means both; 2 only common; 3 only rare;
		#if (length(maf_temp)==length(MAF)){weight=Beta_Weight(MAF,c(1,25));flag=3}
		if (length(maf_temp)==length(MAF)){weight=SKAT:::Beta.Weights(MAF,weights.beta);flag=3}
		
#		cat("MAF: ", MAF, "\n")
#		cat("maf_temp: ", maf_temp, "\n")
#		cat("flag: ", flag, "\n")
#	        cat("weight from Zhaocheng's code : ", weight, "\n")	


	}

	G_w=Matrix(t(t(G)*weight),sparse=TRUE)
	G1_tilde_Ps_G1_tilde = getCovM_nopcg(G1=G_w, G2=G_w, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2.a)
#	print("test3")

	if(length(VarRatio_Vec) > 1){
		VarRatio_12m=Matrix(diag(sqrt(VarRatio_Vec)), sparse=TRUE)
		#cat("VarRatio_12m dim :", dim(VarRatio_12m), "\n")
		#cat("G1_tilde_Ps_G1_tilde dim :", dim(G1_tilde_Ps_G1_tilde), "\n")
		Phi=as.matrix(VarRatio_12m %*% G1_tilde_Ps_G1_tilde %*% VarRatio_12m)
#	print("test4")
	}else{
		Phi = as.matrix(VarRatio_Vec[1] * G1_tilde_Ps_G1_tilde)
	}
	out_kernel=SPA_ER_kernel_related(G,obj, obj.noK, Cutoff=Cutoff, Phi, weight,VarRatio_Vec, mu.a);
#	print("test2")

	zscore.all_1=out_kernel$zscore.all_0* weight
	VarS=out_kernel$VarS*weight^2
	gc()

	r.all = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
	r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
	IDX<-which(r.all >= 0.999)
	if(length(IDX) > 0){
		r.all[IDX]<-0.999	
	}
#	print("test1")
	out=SKAT:::Met_SKAT_Get_Pvalue(Score=zscore.all_1, Phi=as.matrix(Phi), r.corr=r.all, method="optimal.adj",Score.Resampling=NULL)
#	cat("Score from Zhaocheng's code: ", zscore.all_1, "\n")	
#	print(Phi)
	
	list_myfun=list();
	list_myfun$p_skato_old=out$p.value
	rho.val.vec = out$param$rho
	#list_myfun$p_each_old=out$param$p.val.each
	list_myfun$p_each_old = c(out$param$p.val.each[which(rho.val.vec == 1)], out$param$p.val.each[which(rho.val.vec == 0)])

	VarS_org=diag(Phi)		
	vars_inf=which(VarS==Inf)
	if (length(vars_inf)>0){
		VarS[vars_inf] = 0
		zscore.all_1[vars_inf]=0
		Phi[vars_inf,]=0
		Phi[,vars_inf]=0
	}
	if(length(VarS) > 1){
	G2_adj_n=as.matrix(Phi)%*%diag(VarS/VarS_org)	
	}else{
	G2_adj_n=as.matrix(Phi)%*%(VarS/VarS_org)
	}
	mu =out_kernel$mu
	g.sum =out_kernel$g.sum
	q.sum=out_kernel$q.sum
	p.value_burden<-SPAtest:::Saddle_Prob(q.sum , mu=mu, g=g.sum, Cutoff=2,alpha=2.5*10^-6)$p.value


	v1=rep(1,dim(G2_adj_n)[1])
	VarQ=t(v1)%*%G2_adj_n %*%v1


	p.m<-dim(G)[2]
	Q_b=p.m^2 * rowMeans(zscore.all_1)^2

	VarQ_2=Q_b/qchisq(p.value_burden, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

	if (VarQ_2== 0) {r=1} else {r=VarQ/VarQ_2}
	r=min(r,1)
		

	out=SKAT:::Met_SKAT_Get_Pvalue(Score=zscore.all_1, Phi=as.matrix(G2_adj_n), r.corr=r.all, method="optimal.adj",Score.Resampling=NULL)		
	

	list_myfun$p_skato=out$p.value
	#list_myfun$p_each=out$param$p.val.each
	rho.val.vec = out$param$rho
        list_myfun$p_each = c(out$param$p.val.each[which(rho.val.vec == 1)], out$param$p.val.each[which(rho.val.vec == 0)])	

	out=SKAT:::Met_SKAT_Get_Pvalue(Score=zscore.all_1, Phi=as.matrix(G2_adj_n%*%diag(rep(1/r,dim(G2_adj_n)[2]))), r.corr=r.all, method="optimal.adj",Score.Resampling=NULL)	
	list_myfun$p_skato_2=out$p.value
	#list_myfun$p_each_2=out$param$p.val.each
	rho.val.vec = out$param$rho
        list_myfun$p_each_2 = c(out$param$p.val.each[which(rho.val.vec == 1)], out$param$p.val.each[which(rho.val.vec == 0)])


	list_myfun$r=r
		
	p_old=out_kernel$p.old
	p_new=out_kernel$p.new

		
	if (flag==2) {list_myfun$rare_n=0; list_myfun$common_n=length(MAF); list_myfun$rare_mac=0;list_myfun$common_mac=sum(G);}
	if (flag==1){list_myfun$rare_n=length(maf_temp); list_myfun$common_n=length(MAF)-length(maf_temp); list_myfun$rare_mac=sum(G[,maf_temp]);list_myfun$common_mac=sum(G[,-maf_temp]);}
	if (flag==3) {list_myfun$rare_n=length(MAF); list_myfun$common_n=0; list_myfun$rare_mac=sum(G);list_myfun$common_mac=0;}
	list_myfun$mac=MAFsum
	list_myfun$p_single_new=p_new
	list_myfun$p_single_old=p_old
	return (list_myfun);

}

