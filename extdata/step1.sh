seedNum=$2
sampleSize=$1
trait=$3
#ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned-VCFInds_5000_5.bim
#/net/hunt/disk2/zhowei/project/SAIGE_SKAT/computationCost/pheno/output/randomSample_1e+05_seedNum_5.pheno
genoPath=/net/hunt/disk2/zhowei/project/SAIGE_SKAT/computationCost/geno/output/
phenoPath=/net/hunt/disk2/zhowei/project/SAIGE_SKAT/computationCost/pheno/output/
outPath=/net/hunt/disk2/zhowei/project/SAIGE_SKAT/computationCost/geno/SAIGE/step1/output/
#Rscript /net/hunt/disk2/zhowei/project/SAIGE_SKAT/computationCost/geno/SAIGE/step1/jobs/step1_fitNULLGLMM.R	\
Rscript step1_fitNULLGLMM.R	\
	--plinkFile=${genoPath}ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned-VCFInds_${sampleSize}_${seedNum}	\
	--phenoFile=${phenoPath}randomSample_${sampleSize}_seedNum_${seedNum}.pheno	\
	--phenoCol=${trait} \
	--covarColList=Sex,Age_visit0,PC1.wb,PC2.wb,PC3.wb,PC4.wb \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
	--outputPrefix=./${trait}_${sampleSize}_${seedNum}	\
        --nThreads=16 \
        --LOCO=FALSE \
        --IsSparseKin=TRUE      \
        --isCateVarianceRatio=TRUE	\
	--invNormalize=TRUE	\
	--cateVarRatioIndexVec=1,1,1,1,1,1,1,1	\
	--cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5	\
	--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,20.5	\
	--skipModelFitting=FALSE	\
