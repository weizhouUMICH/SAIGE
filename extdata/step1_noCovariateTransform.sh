seedNum=$1
nindep=$2
nfam=$3
tau=$4



#Rscript /net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation_08_2018/jobs/SAIGE_SKATO/step1/jobs/step1_fitNULLGLMM_unit.R \
#for unfinished jobs from 0803, use differnt initial values

#Rscript /net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation_08_2018/jobs/SAIGE_SKATO/step1/jobs/step1_fitNULLGLMM_noCovariateTransform_test.R \
Rscript	step1_fitNULLGLMM.R	\
        --plinkFile=/net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation/geno/plink_saige_step1/output/nfam_${nfam}_nindep_${nindep}_step1_includeMoreRareVariants_poly \
        --phenoFile=/net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation_08_2018/pheno/output/seed_${seedNum}_nfam_${nfam}_nindep_${nindep}_tau_${tau}.pheno_noCov_invnormY \
        --phenoCol=Y \
        --sampleIDColinphenoFile=IND_ID \
        --traitType=quantitative \
        --nThreads=32 \
        --LOCO=FALSE \
        --IsSparseKin=TRUE      \
        --isCateVarianceRatio=TRUE	\
	--invNormalize=FALSE	\
	--outputPrefix=./_noCovariateTransform_noCov \
	--traceCVcutoff=0.00001	
	

#--outputPrefix=/net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation_08_2018/jobs/SAIGE_SKATO/step1/output/seed_${seedNum}_nfam_${nfam}_nindep_${nindep}_tau_${tau}_noCovariateTransform_noCov \
