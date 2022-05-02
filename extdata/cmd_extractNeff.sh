##This script will create two files, which will be specified in the extractNeff.R with --sparseGRMFile and --sparseGRMSampleIDFile

#--sparseGRMFile=sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx
#--sparseGRMSampleIDFile=sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt


Rscript createSparseGRM.R       \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --nThreads=4  \
        --outputPrefix=./output/sparseGRM       \
        --numRandomMarkerforSparseKin=2000      \
        --relatednessCutoff=0.05

#Next extract Nglmm using extractNeff.R
#This script does not generate any output file. The screen outputcan be redirected to the log file. Here I use Nglmm.log. At the end of the Nglmm.log, the value of Nglmm can be found.   

Rscript extractNglmm.R	\
	--useSparseGRMtoFitNULL=TRUE	\
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --nThreads=4    \
        --minMAFforGRM=0.01	\
	--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx	\
	--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt &> Nglmm.log


Rscript extractNglmm.R   \
        --useSparseGRMtoFitNULL=TRUE    \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative        \
        --nThreads=4    \
        --minMAFforGRM=0.01     \
        --sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx  \
        --sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt &> Nglmm_quantitative.log
