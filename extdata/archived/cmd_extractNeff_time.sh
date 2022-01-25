/usr/bin/time -o /net/hunt/disk2/zhowei/project/SAIGE_SKAT/simulation/jobs/log/run.nomosix_1_saige.job.03142018.runinfo.txt -v Rscript /net/hunt/disk2/zhowei/project/SAIGE_SKAT/script/GeneBaseTest_Simu_saige.R --seedNum=1 --nfam=1000


/usr/bin/time -o ./temp.runinfo.txt -v Rscript createSparseGRM.R       \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --nThreads=4  \
        --outputPrefix=./output/sparseGRM       \
        --numRandomMarkerforSparseKin=2000      \
        --relatednessCutoff=0.05
