#obj is the rda. file output from SAIGE step 1
#G1 is genotypes for testing gene, which contains m markers
#G2_cond is G2 in the word document, genotypes for m_cond conditioning marker(s)
#G2_cond_es is beta_2_hat (effect size for the conditioning marker(s))
SAIGE_SKAT_withRatioVec  = function(G1, obj, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, G2_cond = NULL, G2_cond_es, kernel= "linear.weighted", method="optimal.adj", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
, r.corr=0, is_check_genotype=FALSE, is_dosage = TRUE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, sparseSigma = NULL, singleGClambda = 1, mu2 = NULL, adjustCCratioinGroupTest = FALSE, mu=NULL, IsOutputPvalueNAinGroupTestforBinary = FALSE){


#	xt <- proc.time()	
        #check the input genotype G1
        obj.noK = obj$obj.noK
        m = ncol(G1)
        n = nrow(G1)
	AF = colMeans(G1)/2
	
	flipindex = which(AF > 0.5)
	if(length(flipindex) > 0){
		G1[,flipindex] = 2 - G1[,flipindex]
		cat("Note the ", flipindex, "th variants were flipped to use dosages for the minor alleles in gene-based tests\n")
	}

	MAF = colMeans(G1)/2
		
	if(adjustCCratioinGroupTest){
          y = obj$obj.glm.null$y
          N = length(y)
#          obj.noK = obj$obj.noK
	  obj.noK$mu = mu
          obj.noK$res=y - obj.noK$mu 
	  obj.noK$V = mu2
#	xt0 <- proc.time()
          #obj_cc<-SKAT::SKAT_Null_Model(y ~ obj.noK$X1-1, out_type="D", Adjustment = FALSE)
	  obj_cc <- obj$obj_cc	

#	xt1 <- proc.time()
#	print("xt1-xt0")
#	print(xt1-xt0)
          obj_cc$mu=mu
          obj_cc$res=y-obj_cc$mu
          obj_cc$pi_1=obj_cc$mu*(1-obj_cc$mu)
          testtime = system.time({re=Related_ER(G1, obj_cc,  obj.noK, ratioVec=ratioVec, sparseSigma, mac_cutoff = cateVarRatioMinMACVecExclude,Cutoff=2, weights.beta = weights.beta, IsOutputPvalueNAinGroupTestforBinary=IsOutputPvalueNAinGroupTestforBinary)})
	  print("testtime")
	  print(testtime)
	#if(IsOutputPvalueNAinGroupTestforBinary){
	#  	re$p_skato_old2=re$p_skato_old  ##SKATO
	#	re$p_each_old2=Out_List$p_each_old ##SKAT and burden
	#}
	###With additional adjustment (Current robust approach)
	#re$p_skato_2_2=Out_List$p_skato_2 ##SKATO
	#re$p_each_2_2=Out_List$p_each_2 ##SKAT and burden




	}else{ #if(adjustCCratioinGroupTest){


        id_include<-1:n
        out.method<-SKAT:::SKAT_Check_Method(method,r.corr, n=n, m=m)
        method=out.method$method
        r.corr=out.method$r.corr
        IsMeta=out.method$IsMeta
        SKAT:::SKAT_Check_RCorr(kernel, r.corr)

	if(is.null(weights)){
		weights <- SKAT:::Beta.Weights(MAF, weights.beta)
	}
	#cat("weights : ", weights, "\n")	

	indexNeg = NULL
        #if more than 1 marker is left, continue the test
        if(m  >  0){

                #cbind G1 and G2_cond to estimate the variance ratio matrix (m+m_cond) x (m+m_cond)
                if(!is.null(G2_cond)){
                        m_cond = ncol(G2_cond)
                        Zall = cbind(G1, G2_cond)
                }else{
                        Zall = G1
                }



		MACvec_indVec_Zall = getCateVarRatio_indVec(Zall, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)
		rm(Zall)

		GratioMatrixall = getGratioMatrix(MACvec_indVec_Zall, ratioVec)
		if(!is.null(G2_cond)){
			MACvec_indVec = MACvec_indVec_Zall[1:m] 
		}else{
			MACvec_indVec = MACvec_indVec_Zall
		}

                ##summaize the number of markers falling in each MAC category
		markerNumbyMAC = NULL
		for(i in 1:length(cateVarRatioMinMACVecExclude)){
			markerNumbyMAC = c(markerNumbyMAC, sum(MACvec_indVec == i))
		}


                if (kernel == "linear.weighted") {
                        G1 = t(t(G1) * (weights))
                        #Z_tilde = t(t(Z_tilde) * (weights))
                }


                Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])

                if(!is.null(G2_cond)){
                        T2 = as.vector(t(G2_cond) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
                }



                #if no P is provides, use sparseSigma or identity Sigma
                if(is.null(obj$P)){

			G1_tilde_Ps_G1_tilde = getCovM_nopcg(G1=G1, G2=G1, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)

			#check if variance for any marker is negative, remove the variant
                        if(!is.null(G2_cond)){
                                G2_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G2_cond, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
                                G1_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G1, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
                                G2_tilde_Ps_G1_tilde = t(G1_tilde_Ps_G2_tilde) 

		#		print("G2_tilde_Ps_G2_tilde")
		#		print(G2_tilde_Ps_G2_tilde)

		#		print("G1_tilde_Ps_G2_tilde")
                 #               print(G1_tilde_Ps_G2_tilde)

		#		print("GratioMatrixall")
		#		print(GratioMatrixall)


                                G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (G1_tilde_Ps_G2_tilde*(GratioMatrixall[1:m,c((m+1):(m+m_cond))]))%*%(solve(G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])))
		#		print("G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv")
		#		print(G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv)
		#		print("Score")
		#		print(Score)
		#		print("T2")
		#		print(T2)

                                Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2
                                Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (G2_tilde_Ps_G1_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m]))
                                Phi_cond = as.matrix(Phi_cond)
                        }#if(!is.null(G2_cond)){
                        Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])

                }else{ #if(is.null(obj$P)){

			G1_tilde<- G1  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G1)
                        if(!is.null(G2_cond)){
				G2_cond_tilde<- G2_cond  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G2_cond)
                                #G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% solve(t(G2_cond_tilde) %*% (obj$P %*% G2_cond_tilde))
                                G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% getcovM(G2_cond_tilde, G2_cond_tilde, obj$P)

                                Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2

                                Phi_cond = t(G1_tilde) %*% (obj$P %*% G1_tilde) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (t(G2_cond_tilde) %*% (obj$P %*% G1_tilde))
               #                Phi_cond = t(G1_tilde) %*% (obj$P %*% G1_tilde) - (t(G1_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% solve(t(G2_cond_tilde) %*% (obj$P %*% G2_cond_tilde)) %*% (t(G2_cond_tilde) %*% (obj$P %*% G1_tilde))
                	}
                        Phi = t(G1_tilde) %*% (obj$P %*% G1_tilde)

                } #end of else if(is.null(obj$P)){


		#check if variance for each marker is negative, remove the variant
		indexNeg = which(diag(as.matrix(Phi)) <= (.Machine$double.xmin)^(1/4)) 
		#print(Phi)
		#print(diag(as.matrix(Phi)))
		#print(indexNeg)
		#cat("Phi: ", Phi, "\n")
		#cat("Score: ", Score, "\n")
		if(length(indexNeg) > 0){
			Phi = Phi[-indexNeg, -indexNeg]
			Score = Score[-indexNeg]
			if(!is.null(G2_cond)){
				Phi_cond = Phi_cond[-indexNeg, -indexNeg]
				Score_cond = Score_cond[-indexNeg]
			}
			MACvec_indVec = MACvec_indVec[-indexNeg]
			m = m - length(indexNeg)
			markerNumbyMAC = NULL
                	for(i in 1:length(cateVarRatioMinMACVecExclude)){
                        	markerNumbyMAC = c(markerNumbyMAC, sum(MACvec_indVec == i))
                	}
			cat("WARNING: ", indexNeg, " th marker(s) are excluded because of negative variance\n")		
		}

	if(length(Score) > 0){
		m_new = length(Score)
                #Perform the SKAT test
                if(!is.null(G2_cond)){
			
		#	print(Phi_cond)
			#cat("Phi_cond ", Phi_cond, "\n")
		#	cat("Score_cond ", Score_cond, "\n")
			if(m_new == 1){
                        	#if(sum(diag(Phi_cond) < 10^-60) > 0){
                        	if(sum(diag(Phi_cond) < (.Machine$double.xmin)^(1/4)) > 0){
					re_cond = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)

                        	}else{
	                       		re_cond = SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL)
                        	}
			}else{# if(m_new == 1){
				re_cond = SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL)
			}
                }


		
		if(singleGClambda == 1){
                  	Phi = Phi * singleGClambda
		  	Phi = as.matrix(Phi)
			if(m_new == 1){
                 		if(sum(diag(Phi) < (.Machine$double.xmin)^(1/4)) > 0){
                            		re = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)
				}else{	
                  	    		re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
				}
			}else{#if(m_new == 1)
				re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)	
				
         #cat("Score code: ", Score, "\n")
         #print(Phi)

	
			}


		}else{

			if(m_new == 1){
				if(sum(diag(Phi) < (.Machine$double.xmin)^(1/4)) > 0){
					re = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA, P_singlGCadj=NA, GCadjOut=NA)
				}else{	
					re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
					     #cat("Score code: ", Score, "\n")
        print(Phi)

				}
			}else{#if(m_new == 1){

				re =  SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
				     #cat("Score code: ", Score, "\n")
        print(Phi)
			}
			Phi = Phi * singleGClambda
			Phi = as.matrix(Phi)
			if(m_new == 1){	
				if(sum(diag(Phi) < (.Machine$double.xmin)^(1/4)) > 0){
					re$P_singlGCadj = 1
					re$GCadjOut = NA							
				}else{
					re_GCadj = SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)					
					re$P_singlGCadj = re_GCadj$p.value
					re$GCadjOut = re_GCadj
				}

			}else{

		  		re_GCadj = SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
		  		re$P_singlGCadj = re_GCadj$p.value
		  		re$GCadjOut = re_GCadj	

			}

			#}
		}

                if(!is.null(G2_cond)){
                        re$p.value.cond = re_cond$p.value
			re$condOut = re_cond
                }else{
                        re$p.value.cond = NA
                }
		m = length(Score)
	   }else{ #if(length(Score) > 0){
		 #else: no marker is left for test, m = 0
                re = list(p.value = NA, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA)
                #markerNumbyMAC = c(0,0,0,0,0,0)
                markerNumbyMAC = rep(0, length(cateVarRatioMinMACVecExclude))
		m = 0	
	   }	
         }else{ #if(m == 0)

                #else: no marker is left for test, m = 0
                re = list(p.value = NA, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA)
                #markerNumbyMAC = c(0,0,0,0,0,0)
                markerNumbyMAC = rep(0, length(cateVarRatioMinMACVecExclude))
		
        }

        re$IsMeta=TRUE
        re$markerNumbyMAC = markerNumbyMAC
	re$m = m
	re$indexNeg = indexNeg
}#if(adjustCCratioinGroupTest){ else

	#xt1 <- proc.time()
	#print(xt1 - xt)
	cat("re!\n")
        print(re)
        return(re)

}



getGratioMatrix = function(MACvec_indVec, ratioVec){

	numCate = length(ratioVec)
        #cat("MACvec_indVec: ", MACvec_indVec, "\n")

        indMatrix = contr.sum(numCate, contrasts = FALSE)

        GindMatrix = NULL
        for(i in MACvec_indVec){
          GindMatrix = rbind(GindMatrix, indMatrix[i,])
        }

        GratioVec = GindMatrix %*% matrix(ratioVec, ncol=1)
        #mxm
        GratioMatrix = sqrt(GratioVec) %*% t(sqrt(GratioVec))

        return(GratioMatrix)
}




# Function to calculate t(G1_tilde) %*% Sigma_Inverse %*% G1_tilde
# need to get SI_XXVX_inv first
getCovM_nopcg<-function(G1, G2, XV, XXVX_inv, sparseSigma=NULL, mu2){

        # XV<-obj.noK$XV; XXVX_inv<-obj.noK$XXVX_inv
        nSNP2<-ncol(G2)
        nSNP1<-ncol(G1)
        Mat<-matrix(0, nrow=nSNP1, ncol=nSNP2)
        XV_G1 = XV %*% G1

        if(!is.null(sparseSigma)){
        	XV_G2 = XV %*% G2
		SI_XXVX_inv = solve(sparseSigma, XXVX_inv ,sparse=TRUE)
                for(i in 1:nSNP2){
                        SI_G2<-solve(sparseSigma, G2[,i], sparse = TRUE)
                        SI_A_G2<-SI_XXVX_inv %*% XV_G2[,i]
                        A1<- t(G1) %*% (SI_G2 - SI_A_G2)
                        A2<- t(XXVX_inv) %*% (SI_G2 - SI_A_G2)
                        Mat[,i] <-(A1 - t(XV_G1) %*% A2)[,1]
                }
        }else{
                 G2 = G2 * mu2
        	 XV_G2 = XV %*% G2
                 SI_A_G2<-XXVX_inv %*% XV_G2
                 A1<- t(G1) %*% (G2 - SI_A_G2)
                 A2<- t(XXVX_inv) %*% (G2 - SI_A_G2)
                 Mat<-(A1 - t(XV_G1) %*% A2)
                 Mat<-as.matrix(Mat)

        }

        return(Mat)
}





getcovM = function(G1, G2, sparseSigma, mu2 = NULL){

  if(!is.null(sparseSigma)){
   pcginvSigma = NULL
   for(i in 1:ncol(G2)){
     c3<-pcg(sparseSigma, G2[,i])
     pcginvSigma<-cbind(pcginvSigma, c3)
   }
   covM = as.matrix(t(G1) %*% pcginvSigma)
  }else{
      if(!is.null(mu2)){
        G2 = G2 * mu2
      }
      covM = as.matrix(t(G1) %*% G2)
  }
   return(covM)
}


getCateVarRatio_indVec = function(G, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude){

	
	if(ncol(G) > 1){
        	MACvector = colSums(G)
		
	}else{
		MACvector = NULL
		MACvector = c(MACvector, sum(as.vector(G[,1])) )
	}

	MACvector[which(MACvector > nrow(G))] = 2*nrow(G) - MACvector[which(MACvector > nrow(G))]

	#cat("MACvector: ", MACvector, "\n")
        #print(length(MACvector))
        MACvec_indVec = rep(0, length(MACvector))
	#cat("here1 MACvec_indVec: ", MACvec_indVec, "\n")		
	#cat("cateVarRatioMinMACVecExclude: ", cateVarRatioMinMACVecExclude, "\n")
	#cat("cateVarRatioMaxMACVecInclude: ", cateVarRatioMaxMACVecInclude, "\n")
  	numCate = length(cateVarRatioMinMACVecExclude)
#	cat("numCate: ", numCate, "\n")
#	cat("MACvector: ", MACvector, "\n")
    	for(i in 1:(numCate-1)){
		if(i == 1){
			MACvecIndex = which(MACvector >= cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
		}else{
			MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
		}
		if(length(MACvecIndex) > 0){
			MACvec_indVec[MACvecIndex] = i
		}

    	}
#	cat("here2 MACvec_indVec: ", MACvec_indVec, "\n")		
#	cat("here2 length(cateVarRatioMaxMACVecInclude): ", length(cateVarRatioMaxMACVecInclude), "\n")		

    	if(length(cateVarRatioMaxMACVecInclude) == (numCate-1)){
		MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate])
    	}else{
		MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate] & MACvector <= cateVarRatioMaxMACVecInclude[numCate])	
    	}

	if(length(MACvecIndex) > 0){
        	MACvec_indVec[MACvecIndex] = numCate
        }
#	cat("here3 MACvec_indVec: ", MACvec_indVec, "\n")		

        return(MACvec_indVec)
}

###
getVarRatio = function(G, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec){
  if(length(ratioVec) == 1 & ncol(as.matrix(G)) == 1){
    return(ratioVec[1])
  }else{
	G = as.matrix(G)
	if(length(ratioVec) == 1){
	  ratioVec = c(ratioVec, rep(ratioVec[1], ncol(as.matrix(G))))
	}
        if(ncol(G) > 1){
                MACvector = colSums(G)

        }else{
                MACvector = NULL
                MACvector = c(MACvector, sum(as.vector(G[,1])) )
        }

        MACvector[which(MACvector > nrow(G))] = 2*nrow(G) - MACvector[which(MACvector > nrow(G))]

        #cat("MACvector: ", MACvector, "\n")
        #print(length(MACvector))
        MACvec_indVec = rep(0, length(MACvector))
        #cat("here1 MACvec_indVec: ", MACvec_indVec, "\n")
        #cat("cateVarRatioMinMACVecExclude: ", cateVarRatioMinMACVecExclude, "\n")
        #cat("cateVarRatioMaxMACVecInclude: ", cateVarRatioMaxMACVecInclude, "\n")
        numCate = length(cateVarRatioMinMACVecExclude)
#       cat("numCate: ", numCate, "\n")
#       cat("MACvector: ", MACvector, "\n")
        for(i in 1:(numCate-1)){
		if(i == 1){
			MACvecIndex = which(MACvector >= cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
		}else{
                	MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[i] & MACvector <= cateVarRatioMaxMACVecInclude[i])
		}
                if(length(MACvecIndex) > 0){
                        MACvec_indVec[MACvecIndex] = i
                }

        }
#       cat("here2 MACvec_indVec: ", MACvec_indVec, "\n")
#       cat("here2 length(cateVarRatioMaxMACVecInclude): ", length(cateVarRatioMaxMACVecInclude), "\n")

        if(length(cateVarRatioMaxMACVecInclude) == (numCate-1)){
                MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate])
        }else{
                MACvecIndex = which(MACvector > cateVarRatioMinMACVecExclude[numCate] & MACvector <= cateVarRatioMaxMACVecInclude[numCate])
        }

        if(length(MACvecIndex) > 0){
                MACvec_indVec[MACvecIndex] = numCate
        }
#       cat("here3 MACvec_indVec: ", MACvec_indVec, "\n")


	GratioMat = getGratioMatrix(MACvec_indVec, ratioVec)
        #return(MACvec_indVec)
        return(GratioMat)
  }
}

