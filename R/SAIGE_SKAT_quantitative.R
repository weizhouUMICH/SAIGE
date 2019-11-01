#obj is the rda. file output from SAIGE step 1
#G1 is genotypes for testing gene, which contains m markers
#G2_cond is G2 in the word document, genotypes for m_cond conditioning marker(s)
#G2_cond_es is beta_2_hat (effect size for the conditioning marker(s))
SAIGE_SKAT_withRatioVec  = function(G1, obj, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, G2_cond = NULL, G2_cond_es, kernel= "linear.weighted", method="optimal.adj", weights.beta.rare=c(1,25), weights.beta.common=c(0.5,0.5), weightMAFcutoff = 0.01,impute.method="fixed", r.corr=0, is_check_genotype=FALSE, is_dosage = TRUE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, sparseSigma = NULL, mu2 = NULL, adjustCCratioinGroupTest = FALSE, mu=NULL, IsOutputPvalueNAinGroupTestforBinary = FALSE){

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
        id_include<-1:n
        out.method<-SKAT:::SKAT_Check_Method(method,r.corr, n=n, m=m)
        method=out.method$method
        r.corr=out.method$r.corr
        IsMeta=out.method$IsMeta
        SKAT:::SKAT_Check_RCorr(kernel, r.corr)

       if(!is.null(G2_cond)){
         AF_G2 = colMeans(G2_cond)/2
         flipindex_G2 = which(AF_G2 > 0.5)
         if(length(flipindex_G2) > 0){
                G2_cond[,flipindex_G2] = 2 - G2_cond[,flipindex_G2]
                cat("Note the ", flipindex, "th variants of conditioing variants were flipped to use dosages for the minor alleles in gene-based tests\n")
         }
         MAF_G2_cond = colMeans(G2_cond)/2
         MAF = c(MAF, MAF_G2_cond)
       }

	#print(MAF)
	if(length(MAF) > 1){
		weights=rep(0,length(MAF))
		index1 = which(MAF<=weightMAFcutoff)
		if(length(index1) > 0) {weights[which(MAF<=weightMAFcutoff)] = SKAT:::Beta.Weights(MAF[which(MAF<=weightMAFcutoff)],weights.beta.rare)}
		index2 = which(MAF>weightMAFcutoff)
		if(length(index2) > 0) {weights[which(MAF>weightMAFcutoff)] = SKAT:::Beta.Weights(MAF[which(MAF>weightMAFcutoff)],weights.beta.common)}	
	}else{
		if(MAF<=weightMAFcutoff){
			weights = SKAT:::Beta.Weights(MAF,weights.beta.rare)
		}else{
			weights = SKAT:::Beta.Weights(MAF,weights.beta.common)
			
		}
	}


	cat("weights : ", weights, "\n")	
	indexNeg = NULL

        G1_org = G1
        #if more than 1 marker is left, continue the test
        if(m  >  0){
                if(!is.null(G2_cond)){
			if(adjustCCratioinGroupTest){
				G2_cond_org = G2_cond
				G1_org = G1
			}
                        m_cond = ncol(G2_cond)
                        Zall = cbind(G1, G2_cond)
                }else{
			if(adjustCCratioinGroupTest){
                                G1_org = G1
                        }
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
                        G1 = t(t(G1) * (weights[1:m]))
			if(!is.null(G2_cond)){
				G2_cond = t(t(G2_cond) * weights[(m+1):(m+m_cond)])
			}
                }


		Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])

                if(!is.null(G2_cond)){
                        T2 = as.vector(t(G2_cond) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
                }


		if(IsOutputPvalueNAinGroupTestforBinary){
                	#if no P is provides, use sparseSigma or identity Sigma
                	if(is.null(obj$P)){
				G1_tilde_Ps_G1_tilde = getCovM_nopcg(G1=G1, G2=G1, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)

				#check if variance for any marker is negative, remove the variant
                        	if(!is.null(G2_cond)){
                                	G2_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G2_cond, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
                                	G1_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G1, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
					Phi12 = G1_tilde_Ps_G2_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])
                                	G2_tilde_Ps_G1_tilde = t(G1_tilde_Ps_G2_tilde) 
					Phi2 = G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])
                                	G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (Phi12)%*%(solve(Phi2))

                                	Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2
                                	Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (t(Phi12))
                                	Phi_cond = as.matrix(Phi_cond)
                        	}#if(!is.null(G2_cond)){
                        	Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])

                	}
		}

		if(adjustCCratioinGroupTest){
			y = obj$obj.glm.null$y	
			obj_cc <- obj$obj_cc
			obj_cc$mu=mu
			obj_cc$res=y-obj_cc$mu
			obj_cc$pi_1=obj_cc$mu*(1-obj_cc$mu)

                	if(is.null(obj$P)){

				if(!IsOutputPvalueNAinGroupTestforBinary){
					G1_tilde_Ps_G1_tilde = getCovM_nopcg(G1=G1, G2=G1, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
                        		Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])
				}
				Phi_ccadj = SPA_ER_kernel_related_Phiadj(G1_org, obj_cc, obj.noK, Cutoff=2, Phi, weights[1:m], VarRatio_Vec = as.vector(GratioMatrixall[1:m,1]), mu)

				#check if variance for any marker is negative, remove the variant
                        	if(!is.null(G2_cond)){
					if(!IsOutputPvalueNAinGroupTestforBinary){
                                		G2_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G2_cond, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
						Phi2 = G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])
                                		G1_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G1, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
						Phi12 = G1_tilde_Ps_G2_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])
                                		G2_tilde_Ps_G1_tilde = t(G1_tilde_Ps_G2_tilde)
						G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (Phi12)%*%(solve(Phi2)) 

                                		Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2
                                		Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (t(Phi12))
                                		Phi_cond = as.matrix(Phi_cond)
					}


					Phi2_ccadj = SPA_ER_kernel_related_Phiadj(G2_cond_org, obj_cc, obj.noK, Cutoff=2, Phi2, weights[((m+1):(m+m_cond))], VarRatio_Vec = as.vector(GratioMatrixall[c((m+1):(m+m_cond)),1]), mu)			
					Phi12_ccadj_val = Phi_ccadj$scaleFactor %*% Phi12 %*% Phi2_ccadj$scaleFactor
					G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv_ccadj = Phi12_ccadj_val%*%solve(Phi2_ccadj$val)
					Score_cond_ccadj = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv_ccadj %*% T2
					Phi_cond_ccadj = Phi_ccadj$val - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv_ccadj %*% t(Phi12_ccadj_val)

                        	}#if(!is.null(G2_cond)){

			}
		}


		indexNeg = which(diag(as.matrix(Phi)) <= (.Machine$double.xmin)^(1/4)) 
                #Score = as.vector(t(G1) %*% matrix(obj$residuals, ncol=1))/as.numeric(obj$theta[1])
		if(length(indexNeg) > 0){
                	Phi = Phi[-indexNeg, -indexNeg]
                        Score = Score[-indexNeg]
                                #if(!is.null(G2_cond)){
                                #        Phi_cond = Phi_cond[-indexNeg, -indexNeg]
                                #        Score_cond = Score_cond[-indexNeg]
                                #}
                        MACvec_indVec = MACvec_indVec[-indexNeg]
                        m = m - length(indexNeg)
                        markerNumbyMAC = NULL
                        for(i in 1:length(cateVarRatioMinMACVecExclude)){
                        	markerNumbyMAC = c(markerNumbyMAC, sum(MACvec_indVec == i))
                        }
                        cat("WARNING: ", indexNeg, " th marker(s) are excluded because of negative variance\n")
                }


		

		if(IsOutputPvalueNAinGroupTestforBinary){
			#check if variance for each marker is negative, remove the variant
			if(length(indexNeg) > 0){
			# 	Phi = Phi[-indexNeg, -indexNeg]
			#	Score = Score[-indexNeg]
				if(!is.null(G2_cond)){
					T2 = T2[-indexNeg]
					Phi_cond = Phi_cond[-indexNeg, -indexNeg]
					Score_cond = Score_cond[-indexNeg]
				}
			}

		}


		if(adjustCCratioinGroupTest){
                        if(length(indexNeg) > 0){
                                Phi_ccadj = Phi_ccadj[-indexNeg, -indexNeg]
                                if(!is.null(G2_cond)){
                                        Phi_cond_ccadj = Phi_cond_ccadj[-indexNeg, -indexNeg]
                                        Score_cond_ccadj = Score_cond_ccadj[-indexNeg]
                                }
                        }

		}

		if(length(Score) > 0){
			m_new = length(Score)

			if(IsOutputPvalueNAinGroupTestforBinary){
                        	if(m_new == 1){
                                                #if(sum(diag(Phi_cond) < 10^-60) > 0){
                                	if(sum(diag(Phi) < (.Machine$double.xmin)^(1/4)) > 0){
                                        	re = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)

                                        }else{
                                                re = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL))
						if(class(re) == "try-error"){
							re = list(p.value = NA)
						}	
                                        }
                         	}else{# if(m_new == 1){

					
                                	re = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
)
					if(class(re) == "try-error"){
                                                        re = list(p.value = NA)
                                                }					

                                }


                         }


			if(adjustCCratioinGroupTest){
				if(m_new == 1){
                                                #if(sum(diag(Phi_cond) < 10^-60) > 0){
                                        if(sum(diag(Phi_ccadj$val) < (.Machine$double.xmin)^(1/4)) > 0){

						if(!IsOutputPvalueNAinGroupTestforBinary){
							Out_ccadj = list(p.value = 1, param = NULL)
							#re=list(p.value_cc = 1, param.ccadj=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)
							re = list(Out_ccadj = Out_ccadj)

						}else{
							Out_ccadj = list(p.value = 1)
							re$Out_ccadj = Out_ccadj

						}

                                        }else{
						if(!IsOutputPvalueNAinGroupTestforBinary){
							re = list()
						}

                                                        re_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi_ccadj$val, r.corr=r.corr, method=method, Score.Resampling=NULL))
							if(class(re_ccadj) == "try-error"){
								re_ccadj = list(p.value = NA)
							}	
							re$Out_ccadj = re_ccadj
                                        }
                                }else{# if(m_new == 1){
					if(!IsOutputPvalueNAinGroupTestforBinary){
                                        	re = list()
                                        }
                                                re_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi_ccadj$val, r.corr=r.corr, method=method, Score.Resampling=NULL))
						if(class(re_ccadj) == "try-error"){
							re_ccadj = list(p.value=NA)
						}	
						re$Out_ccadj = re_ccadj

                                }
			}




                	#Perform the SKAT test
                	if(!is.null(G2_cond)){
				if(IsOutputPvalueNAinGroupTestforBinary){			
					if(m_new == 1){
                        			#if(sum(diag(Phi_cond) < 10^-60) > 0){
                        			if(sum(diag(Phi_cond) < (.Machine$double.xmin)^(1/4)) > 0){
							re_cond = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)

                        			}else{
	                       				re_cond = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL))
							  if(class(re_cond) == "try-error"){
                                                        re_cond = list(p.value=NA)
                                                	}
                        			}
					}else{# if(m_new == 1){
						re_cond = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL))
						if(class(re_cond) == "try-error"){
                                                        re_cond = list(p.value=NA)
                                                 }
					}
					re$condOut = re_cond

				}



				if(adjustCCratioinGroupTest){
					if(m_new == 1){
                                                #if(sum(diag(Phi_cond) < 10^-60) > 0){
                                                if(sum(diag(Phi_cond_ccadj) < (.Machine$double.xmin)^(1/4)) > 0){
                                                        re_cond_ccadj = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)

                                                }else{
                                                        re_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=r.corr, method=method, Score.Resampling=NULL))
							if(class(re_cond_ccadj) == "try-error"){
                                                        re_cond_ccadj = list(p.value=NA)
                                                 	}
                                                }
                                        }else{# if(m_new == 1){
                                                re_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=r.corr, method=method, Score.Resampling=NULL))
						if(class(re_cond_ccadj) == "try-error"){
                                                        re_cond_ccadj = list(p.value=NA)
                                                }
                                        }

					re$condOut_ccadj = re_cond_ccadj

				}
                	}	
		
			m = length(Score)
	   	}else{ #if(length(Score) > 0){
		 		#else: no marker is left for test, m = 0
			
                	re = list(p.value = NA, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA, p.value.cond.ccadj=NA)
			re =list()	
                	markerNumbyMAC = rep(0, length(cateVarRatioMinMACVecExclude))
			m = 0	
	  	}


	
         }else{ #if(m == 0)

                #else: no marker is left for test, m = 0
                re = list(p.value = NA, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA, p.value.cond=NA, p.value.cond_ccadj=NA)
                #markerNumbyMAC = c(0,0,0,0,0,0)
                markerNumbyMAC = rep(0, length(cateVarRatioMinMACVecExclude))
		
        }

        re$IsMeta=TRUE
        re$markerNumbyMAC = markerNumbyMAC
	re$m = m
	re$indexNeg = indexNeg

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

