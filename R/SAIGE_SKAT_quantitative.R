#obj is the rda. file output from SAIGE step 1
#G1 is genotypes for testing gene, which contains m markers
#G2_cond is G2 in the word document, genotypes for m_cond conditioning marker(s)
#G2_cond_es is beta_2_hat (effect size for the conditioning marker(s))
SAIGE_SKAT_withRatioVec  = function(G1, obj, y, X, tauVec, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, ratioVec, G2_cond = NULL, G2_cond_es, kernel= "linear.weighted", method="optimal.adj", weights.beta.rare=c(1,25), weights.beta.common=c(0.5,0.5), weightMAFcutoff = 0.01,impute.method="fixed", r.corr=0, is_check_genotype=FALSE, is_dosage = TRUE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID = NULL, sparseSigma = NULL, mu2 = NULL, adjustCCratioinGroupTest = FALSE, mu=NULL, IsOutputPvalueNAinGroupTestforBinary = FALSE, weights_specified = NULL, weights_for_G2_cond = NULL, weightsIncludeinGroupFile = FALSE, IsOutputBETASEinBurdenTest=FALSE,  method_to_CollapseUltraRare = "absence_or_presence",  MACCutoff_to_CollapseUltraRare = 10, DosageCutoff_for_UltraRarePresence = 0.5){
	#offset = obj$offset
        #check the input genotype G1
        obj.noK = obj$obj.noK
	obj.noK$y = y
        m = ncol(G1)
        n = nrow(G1)
	MACvec = colSums(G1)
	MAF = MACvec/(2*n)
        AF = MAF
	#AF = colMeans(G1)/2	
	#flipindex = which(AF > 0.5)
	#if(length(flipindex) > 0){
	#	G1[,flipindex] = 2 - G1[,flipindex]
	#	MACvec[flipindex] = 2*n - MACvec[flipindex]
	#	cat("Note the ", flipindex, "th variants were flipped to use dosages for the minor alleles in gene-based tests\n")
	#}

	#MAF = colMeans(G1)/2


        macle10Index = which(MACvec <= MACCutoff_to_CollapseUltraRare)
        if(method_to_CollapseUltraRare != "" & length(macle10Index) > 0){
                #Gnew = rowSums(Gmat[,macle10Index,drop=F])/(length(macle10Index))
		G1rare=G1[,macle10Index, drop=F]
                if(method_to_CollapseUltraRare == "absence_or_presence"){
                	Gnew = rowSums(G1rare)
                        Gnew[which(Gnew >= DosageCutoff_for_UltraRarePresence)] = 1
                }else if(method_to_CollapseUltraRare == "sum_geno"){ #####NOT active
      
			##determine the weights of ultra rare variants
			MAFle10 = MAF[macle10Index]
			if(!weightsIncludeinGroupFile){
                		if(length(MAFle10) > 1){
                        		weights_MAFle10=rep(0,length(MAFle10))
                        		index1 = which(MAFle10<=weightMAFcutoff)
                        		if(length(index1) > 0) {weights_MAFle10[which(MAFle10<=weightMAFcutoff)] = SKAT:::Beta.Weights(MAFle10[which(MAFle10<=weightMAFcutoff)],weights.beta.rare)}
                        		index2 = which(MAFle10>weightMAFcutoff)
                        		if(length(index2) > 0) {weights_MAFle10[which(MAFle10>weightMAFcutoff)] = SKAT:::Beta.Weights(MAFle10[which(MAFle10>weightMAFcutoff)],weights.beta.common)}
                		}else{
                        		if(MAFle10<=weightMAFcutoff){
                                		weights_MAFle10 = SKAT:::Beta.Weights(MAFle10,weights.beta.rare)
                        		}else{
                                		weights_MAFle10 = SKAT:::Beta.Weights(MAFle10,weights.beta.common)

                        		}
                		}
        		}else{
                		weights_MAFle10 = weights_specified[macle10Index]
                		cat("weights is specified in the group file for ultra rare variants.\n")
                        }
			Gnew = rep(0, n)
			for (i in 1:length(MAFle10)){
				Gnew = Gnew + weights_MAFle10[i] * G1rare[,i]
			}	
                } #####NOT active

		if(length(macle10Index) < m){
                        G1 = cbind(Gnew, G1[,-macle10Index, drop=F])
                }else{
                        G1_sub = NULL
                        G1 = cbind(Gnew, G1_sub)
                }

		m = ncol(G1)
        	MACvec = colSums(G1)
        	MAF = MACvec/(2*n)
        	AF = MAF
        }

        id_include<-1:n
        out.method<-SKAT:::SKAT_Check_Method(method,r.corr, n=n, m=m)
        method=out.method$method
        r.corr=out.method$r.corr
        IsMeta=out.method$IsMeta
        SKAT:::SKAT_Check_RCorr(kernel, r.corr)

       if(!is.null(G2_cond)){
         #AF_G2 = colMeans(G2_cond)/2
	 MACvec_cond = colSums(G2_cond)
	 AF_G2 = MACvec_cond/(2*n) 
         flipindex_G2 = which(AF_G2 > 0.5)
         if(length(flipindex_G2) > 0){
                G2_cond[,flipindex_G2] = 2 - G2_cond[,flipindex_G2]
		MACvec_cond[flipindex_G2] = 2*n - MACvec_cond[flipindex_G2]
                cat("Note the ", flipindex, "th variants of conditioing variants were flipped to use dosages for the minor alleles in gene-based tests\n")
         }
         MAF_G2_cond = colMeans(G2_cond)/2
         MAF = c(MAF, MAF_G2_cond)
       }

	#print(MAF)
	if(!weightsIncludeinGroupFile){
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
	}else{
		weights = weights_specified
		cat("weights is specified in the group file.\n")

		if(!is.null(G2_cond)){
			if(!is.null(weights_for_G2_cond)){
				weights = c(weights, weights_for_G2_cond)
			}else{
				stop("weights is not specified for the conditioning marker(s)\n")
			} 
		}
	}

	cat("weights: ", weights, "\n")	
	indexNeg = NULL
        print("DEBUG1")
        G1_org = G1
        #if more than 1 marker is left, continue the test
        if(m  >  0){
                if(!is.null(G2_cond)){
			if(adjustCCratioinGroupTest){
				G2_cond_org = G2_cond
				G1_org = G1
			}
                        m_cond = ncol(G2_cond)
                        #Zall = cbind(G1, G2_cond)
			MACvec = c(MACvec, MACvec_cond)
                }else{
			if(adjustCCratioinGroupTest){
                                G1_org = G1
                        }
                        #Zall = G1
                }
		
		#MACvec_indVec_Zall = getCateVarRatio_indVec(Zall, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)
		#rm(Zall)
		print("DEBUG2")
		MACvec_indVec_Zall = getCateVarRatio_indVec(MACvector=MACvec, cateVarRatioMinMACVecExclude=cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude=cateVarRatioMaxMACVecInclude)
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
					#Phi12 = G1_tilde_Ps_G2_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])
					Phi12 = G1_tilde_Ps_G2_tilde * (GratioMatrixall[1:m, c((m+1):(m+m_cond))])
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
			#y = obj$y	
			obj_cc <- SKAT::SKAT_Null_Model(y ~ X-1, out_type="D", Adjustment = FALSE) 
			obj_cc$mu=mu
			obj_cc$res=y-obj_cc$mu
			obj_cc$pi_1=obj_cc$mu*(1-obj_cc$mu)

                	if(is.null(obj$P)){
                print("DEBUG3")

				if(!IsOutputPvalueNAinGroupTestforBinary){
					G1_tilde_Ps_G1_tilde = getCovM_nopcg(G1=G1, G2=G1, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
                print("DEBUG4")
                        		Phi = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m])
				}
				Phi_ccadj = SPA_ER_kernel_related_Phiadj(G1_org, obj_cc, obj.noK, Cutoff=2, Phi, weights[1:m], VarRatio_Vec = as.vector(GratioMatrixall[1:m,1]), mu, sparseSigma)

				#check if variance for any marker is negative, remove the variant
                        	if(!is.null(G2_cond)){
					if(!IsOutputPvalueNAinGroupTestforBinary){
                                		G2_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G2_cond, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
						Phi2 = G2_tilde_Ps_G2_tilde*(GratioMatrixall[c((m+1):(m+m_cond)),c((m+1):(m+m_cond))])
                                		G1_tilde_Ps_G2_tilde = getCovM_nopcg(G1=G1, G2=G2_cond, XV=obj.noK$XV, XXVX_inv=obj.noK$XXVX_inv, sparseSigma = sparseSigma, mu2 = mu2)
						#Phi12 = G1_tilde_Ps_G2_tilde * (GratioMatrixall[c((m+1):(m+m_cond)), 1:m])
						Phi12 = G1_tilde_Ps_G2_tilde * (GratioMatrixall[1:m, c((m+1):(m+m_cond))])
                                		G2_tilde_Ps_G1_tilde = t(G1_tilde_Ps_G2_tilde)
						G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv = (Phi12)%*%(solve(Phi2)) 

                                		Score_cond = Score - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% T2
                                		Phi_cond = G1_tilde_Ps_G1_tilde*(GratioMatrixall[1:m,1:m]) - G1_tilde_P_G2_tilde_G2_tilde_P_G2_tilde_inv %*% (t(Phi12))
                                		Phi_cond = as.matrix(Phi_cond)
					}


					Phi2_ccadj = SPA_ER_kernel_related_Phiadj(G2_cond_org, obj_cc, obj.noK, Cutoff=2, Phi2, weights[((m+1):(m+m_cond))], VarRatio_Vec = as.vector(GratioMatrixall[c((m+1):(m+m_cond)),1]), mu, sparseSigma)			
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
			print("MACvec_indVec")	
			print(MACvec_indVec)
                        MACvec_indVec = MACvec_indVec[-indexNeg]
                        m = m - length(indexNeg)
                        markerNumbyMAC = NULL
                        for(i in 1:length(cateVarRatioMinMACVecExclude)){
                        	markerNumbyMAC = c(markerNumbyMAC, sum(MACvec_indVec == i))
                        }
                        cat("WARNING: ", indexNeg, " th marker(s) are excluded because of negative variance\n")
			print("MACvec_indVec")	
			print(MACvec_indVec)

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
				#print("Phi_ccadj")
                        	#print(Phi_ccadj)
                                Phi_ccadj$scaleFactor = Phi_ccadj$scaleFactor[-indexNeg, -indexNeg]
                                Phi_ccadj$val = Phi_ccadj$val[-indexNeg, -indexNeg]
                                if(!is.null(G2_cond)){
                                        Phi_cond_ccadj$val = Phi_cond_ccadj$val[-indexNeg, -indexNeg]
                                        Phi_cond_ccadj$scaleFactor = Phi_cond_ccadj$scaleFactor[-indexNeg, -indexNeg]
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
							re_btemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=1, method=method, Score.Resampling=NULL))
							re_stemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=0, method=method, Score.Resampling=NULL))
							if(class(re_btemp) == "try-error" | class(re_stemp) == "try-error"){	
								re = list(p.value = NA)
							}else{
								re = list(p.value = 2*min(re_btemp$p.value, re_stemp$p.value, 0.5), param = list())
								re$param = list(p.val.each = c(re_btemp$p.value, re_stemp$p.value), rho=c(1,0))

							}
						}	
                                        }
                         	}else{# if(m_new == 1){

					
                                	re = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=NULL)
)
					if(class(re) == "try-error"){
                                        	re_btemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=1, method=method, Score.Resampling=NULL)) 
                                                re_stemp = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=0, method=method, Score.Resampling=NULL)) 

                                                if(class(re_btemp) == "try-error" | class(re_stemp) == "try-error"){
                                                	re = list(p.value = NA)
                                                }else{
                                                	re = list(p.value = 2*min(re_btemp$p.value, re_stemp$p.value, 0.5), param=list())
							re$param = list(p.val.each = c(re_btemp$p.value, re_stemp$p.value), rho=c(1,0))


                                                }
                                        }					

                                }
				if(IsOutputBETASEinBurdenTest){
					re$Phi_sum = sum(Phi)
					re$Score_sum = sum(Score)
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
								re_btemp_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi_ccadj$val, r.corr=1, method=method, Score.Resampling=NULL)) 
								re_stemp_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi_ccadj$val, r.corr=0, method=method, Score.Resampling=NULL)) 

                                                        	if(class(re_btemp_ccadj) == "try-error" | class(re_stemp_ccadj) == "try-error"){
                                                                	re_ccadj = list(p.value = NA)
                                                        	}else{
                                                                	re_ccadj = list(p.value = 2*min(re_btemp_ccadj$p.value, re_stemp_ccadj$p.value, 0.5), param=list())
									re_ccadj$param = list(p.val.each = c(re_btemp_ccadj$p.value, re_stemp_ccadj$p.value), rho=c(1,0))
                                                        	}

							}

	
							re$Out_ccadj = re_ccadj
                                        }
                                }else{# if(m_new == 1){
					if(!IsOutputPvalueNAinGroupTestforBinary){
                                        	re = list()
                                        }
                                                re_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi_ccadj$val, r.corr=r.corr, method=method, Score.Resampling=NULL))
						if(class(re_ccadj) == "try-error"){
							re_btemp_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi_ccadj$val, r.corr=1, method=method, Score.Resampling=NULL))
                                                	re_stemp_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi_ccadj$val, r.corr=0, method=method, Score.Resampling=NULL))

                                                        if(class(re_btemp_ccadj) == "try-error" | class(re_stemp_ccadj) == "try-error"){
                                                        	re_ccadj = list(p.value = NA)
                                                        }else{
                                                                re_ccadj = list(p.value = 2*min(re_btemp_ccadj$p.value, re_stemp_ccadj$p.value, 0.5), param=list())
								re_ccadj$param = list(p.val.each = c(re_btemp_ccadj$p.value, re_stemp_ccadj$p.value), rho=c(1,0))

                                                        }

						}	
						re$Out_ccadj = re_ccadj

                                }

				if(IsOutputBETASEinBurdenTest){
					re$Phi_ccadj_sum = sum(Phi_ccadj$val)
					re$Score_sum = sum(Score)
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
                                                        	re_btemp_cond = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=1, method=method, Score.Resampling=NULL)) 
								re_stemp_cond = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=0, method=method, Score.Resampling=NULL))

                                                        	if(class(re_btemp_cond) == "try-error" | class(re_stemp_cond) == "try-error"){
                                                                	re_cond = list(p.value = NA)
                                                        	}else{
                                                                	re_cond = list(p.value = 2*min(re_btemp_cond$p.value, re_stemp_cond$p.value, 0.5), param=list())
									re_cond$param = list(p.val.each = c(re_btemp_cond$p.value, re_stemp_cond$p.value), rho=c(1,0))


                                                        	}
                                                	}
                        			}
					}else{# if(m_new == 1){
						re_cond = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=r.corr, method=method, Score.Resampling=NULL))
						if(class(re_cond) == "try-error"){
							re_btemp_cond = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=1, method=method, Score.Resampling=NULL))
                                                        re_stemp_cond = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond, Phi=Phi_cond, r.corr=0, method=method, Score.Resampling=NULL))

                                                        if(class(re_btemp_cond) == "try-error" | class(re_stemp_cond) == "try-error"){
                                                        	re_cond = list(p.value = NA)
                                                        }else{
                                                                re_cond = list(p.value = 2*min(re_btemp_cond$p.value, re_stemp_cond$p.value, 0.5), param=list())
								re_cond$param = list(p.val.each = c(re_btemp_cond$p.value, re_stemp_cond$p.value), rho=c(1,0))
                                                        }

                                                 }
					}
					re$condOut = re_cond
					if(IsOutputBETASEinBurdenTest){
						re$Score_cond_sum = sum(Score_cond)
						re$Phi_cond_sum = sum(Phi_cond)
					}

				}



				if(adjustCCratioinGroupTest){
					if(m_new == 1){
                                                #if(sum(diag(Phi_cond) < 10^-60) > 0){
                                                if(sum(diag(Phi_cond_ccadj) < (.Machine$double.xmin)^(1/4)) > 0){
                                                        re_cond_ccadj = list(p.value = 1, param=NA, p.value.resampling=NA, pval.zero.msg=NA, Q=NA)

                                                }else{
                                                        re_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=r.corr, method=method, Score.Resampling=NULL))
							if(class(re_cond_ccadj) == "try-error"){
								re_btemp_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=1, method=method, Score.Resampling=NULL))
                                                        	re_stemp_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=0, method=method, Score.Resampling=NULL))

                                                        	if(class(re_btemp_cond_ccadj) == "try-error" | class(re_stemp_cond_ccadj) == "try-error"){
                                                                	re_cond_ccadj = list(p.value = NA)
                                                        	}else{
                                                                	re_cond_ccadj = list(p.value = 2*min(re_btemp_cond_ccadj$p.value, re_stemp_cond_ccadj$p.value, 0.5), param=list())
									re_cond_ccadj$param = list(p.val.each = c(re_btemp_cond_ccadj$p.value, re_stemp_cond_ccadj$p.value), rho=c(1,0))
                                                        	}

                                                 	}
                                                }
                                        }else{# if(m_new == 1){
                                                re_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=r.corr, method=method, Score.Resampling=NULL))
							if(class(re_cond_ccadj) == "try-error"){
                                                                re_btemp_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=1, method=method, Score.Resampling=NULL))
                                                                re_stemp_cond_ccadj = try(SKAT:::Met_SKAT_Get_Pvalue(Score=Score_cond_ccadj, Phi=Phi_cond_ccadj, r.corr=0, method=method, Score.Resampling=NULL))

                                                                if(class(re_btemp_cond_ccadj) == "try-error" | class(re_stemp_cond_ccadj) == "try-error"){
                                                                        re_cond_ccadj = list(p.value = NA)
                                                                }else{
                                                                        re_cond_ccadj = list(p.value = 2*min(re_btemp_cond_ccadj$p.value, re_stemp_cond_ccadj$p.value, 0.5), param = list())
									re_cond_ccadj$param = list(p.val.each = c(re_btemp_cond_ccadj$p.value, re_stemp_cond_ccadj$p.value), rho=c(1,0)) 
                                                                }
                                                        }	
                                        }
					re$condOut_ccadj = re_cond_ccadj
					if(IsOutputBETASEinBurdenTest){
						re$Score_cond_ccadj_sum = sum(Score_cond_ccadj)
						re$Phi_cond_ccadj_sum = sum(Phi_cond_ccadj)
					}	

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

getCateVarRatio_indVec = function(MACvector = NULL, G = NULL, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude){

     if(!is.null(G)){
        if(ncol(G) > 1){
                MACvector = colSums(G)

        }else{
                MACvector = NULL
                MACvector = c(MACvector, sum(as.vector(G[,1])) )
        }

        MACvector[which(MACvector > nrow(G))] = 2*nrow(G) - MACvector[which(MACvector > nrow(G))]
     }

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

        return(MACvec_indVec)
}
