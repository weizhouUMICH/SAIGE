FROM ubuntu:16.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt-get update && \
    apt-get install curl make gcc g++ cmake gfortran libreadline-dev \
    libz-dev libbz2-dev liblzma-dev libpcre3-dev libssl-dev libcurl4-openssl-dev \
    libopenblas-dev default-jre unzip libboost-all-dev \
    libpng-dev libcairo2-dev tabix --yes && \
    apt-get clean && \
    curl -O https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz && \
    tar xvzf R-3.5.1.tar.gz && \
    cd R-3.5.1 && \
    ./configure --with-x=no --with-blas="-lopenblas" && \
    make && \
    mkdir -p /usr/local/lib/R/lib && \
    make install && \
    cd .. && \
    rm -rf R-3.5.1*

RUN apt-get install python3-pip --yes && \
    pip3 install cget

RUN apt install wget
#RUN src_branch=master  && \
#	repo_src_url=https://github.com/weizhouUMICH/SAIGE && \
#	git clone --depth 0 -b $src_branch $repo_src_url


#RUN wget https://github.com/weizhouUMICH/SAIGE/archive/master.zip && \
#        unzip master.zip && \
#        mv SAIGE-master SAIGE 
RUN apt-get update

RUN apt-get install -y libxml2-dev


#RUN Rscript SAIGE/extdata/install_packages.R
ADD install_packages.R /usr/local/bin/


RUN Rscript /usr/local/bin/install_packages.R


#RUN Rscript -e 'devtools::install_github("leeshawn/SKAT")'
#RUN echo "update 0.37 new"

#RUN Rscript -e 'devtools::install_github("leeshawn/SPAtest")'
#RUN echo "update 0.43.2 new2"

RUN Rscript -e 'devtools::install_github("leeshawn/MetaSKAT")'
#RUN echo "update 0.44.2"
#RUN echo "update 0.39.2_new" &&  wget https://github.com/weizhouUMICH/SAIGE/archive/0.39.2.zip && \
#        unzip 0.39.2.zip && \
#        mv SAIGE-0.39.2 SAIGE
#RUN echo "update 0.44.2.b hello" &&  wget https://github.com/weizhouUMICH/SAIGE/archive/master.zip && \
#	unzip master.zip && \
#	mv SAIGE-master SAIGE
#RUN echo "update3"

RUN export MKL_NUM_THREADS=1
RUN export MKL_DYNAMIC=false
RUN export OMP_NUM_THREADS=1
RUN export OMP_DYNAMIC=false

RUN apt-get update

RUN apt-get install -y git
RUN echo "hello world13"

#ADD install_packages_2.R /usr/local/bin/
#RUN Rscript /usr/local/bin/install_packages_2.R

#ENV SRC_BRANCH master
#ENV SRC_BRANCH ScoreSPARcpp_v0.44.3
#ENV SRC_BRANCH 0.43.1i
#RUN echo "hello world16"
#ADD install_packages_3.R /usr/local/bin/
#RUN Rscript /usr/local/bin/install_packages_3.R

ENV SRC_BRANCH new_dev2
#ENV SRC_BRANCH master
ENV REPO_SRC_URL https://github.com/weizhouUMICH/SAIGE

RUN echo "hello world19" 

RUN rm -rf SAIGE
#RUN echo "hello world562"

RUN git clone --depth 1 -b $SRC_BRANCH $REPO_SRC_URL


RUN mv SAIGE/extdata/install_packages.R /usr/local/bin/
RUN Rscript /usr/local/bin/install_packages.R

RUN R CMD INSTALL SAIGE


RUN Rscript SAIGE/extdata/createSparseGRM.R --help

RUN Rscript SAIGE/extdata/step1_fitNULLGLMM.R --help
RUN echo "hello world56"

RUN Rscript SAIGE/extdata/step2_SPAtests.R --help
RUN echo "hello world5"



#RUN Rscript SAIGE/extdata/extractNglmm.R --help
#RUN echo "hello world5"

#RUN Rscript SAIGE/extdata/install_packages.R \
#    && R CMD INSTALL SAIGE \
#    && mv SAIGE/extdata/step1_fitNULLGLMM.R SAIGE/extdata/step2_SPAtests.R \
#    SAIGE/extdata/createSparseGRM.R /usr/local/bin/ \

#RUN chmod a+x /usr/local/bin/step1_fitNULLGLMM.R
#RUN chmod a+x /usr/local/bin/step2_SPAtests.R


#ENV PATH "/usr/local/bin/:$PATH"

#ADD SAIGE/extdata/step1_fitNULLGLMM.R SAIGE/extdata/step2_SPAtests.R SAIGE/extdata/createSparseGRM.R /usr/local/bin/
#RUN mv SAIGE/extdata/extractNglmm.R  SAIGE/extdata/step1_fitNULLGLMM.R SAIGE/extdata/step2_SPAtests.R SAIGE/extdata/createSparseGRM.R /usr/local/bin/ 

RUN mv SAIGE/extdata/step1_fitNULLGLMM.R SAIGE/extdata/step2_SPAtests.R SAIGE/extdata/createSparseGRM.R /usr/local/bin/

RUN chmod a+x /usr/local/bin/step1_fitNULLGLMM.R
RUN chmod a+x /usr/local/bin/step2_SPAtests.R
RUN chmod a+x /usr/local/bin/createSparseGRM.R
#RUN chmod a+x /usr/local/bin/extractNglmm.R

#RUN Rscript /usr/local/bin/createSparseGRM.R  --help

#RUN Rscript /usr/local/bin/step1_fitNULLGLMM.R --help

#RUN Rscript /usr/local/bin/step2_SPAtests.R --help

RUN createSparseGRM.R  --help

RUN step1_fitNULLGLMM.R --help

RUN step2_SPAtests.R --help

#RUN extractNglmm.R --help
#RUN rm -r SAIGE
