FROM ubuntu:16.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV PATH "/usr/local/bin/:$PATH"

RUN apt-get update && \
    # Install apt dependencies
    apt-get install --yes --no-install-recommends \
    cmake \
    curl \
    default-jre-headless \
    g++ \
    gcc \
    gfortran \
    libboost-all-dev \
    libbz2-dev \
    libcairo2-dev  \
    libcurl4-openssl-dev \
    liblzma-dev \
    libopenblas-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libssl-dev \
    libxml2-dev \
    libz-dev \
    make \
    python3-dev \
    python3-pip \
    python3-setuptools \
    tabix && \
    rm -rf /var/lib/apt/lists/* && \
    # Install Python dependencies
    python3 -m pip install pip --upgrade && \
    python3 -m pip install cget && \
    # Install R
    curl -O https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz && \
    tar xvzf R-3.5.1.tar.gz && \
    cd R-3.5.1 && \
    ./configure --with-x=no --with-blas="-lopenblas" && \
    make && \
    mkdir -p /usr/local/lib/R/lib && \
    make install && \
    cd .. && \
    rm -rf R-3.5.1*

COPY . SAIGE

# Install SAIGE and R dependencies
RUN Rscript SAIGE/extdata/install_packages.R && \
    R CMD INSTALL SAIGE && \
    cp SAIGE/extdata/createSparseGRM.R /usr/local/bin/createSparseGRM.R && \
    chmod +x /usr/local/bin/createSparseGRM.R && \
    cp SAIGE/extdata/step1_fitNULLGLMM.R /usr/local/bin/step1_fitNULLGLMM.R && \
    chmod +x /usr/local/bin/step1_fitNULLGLMM.R && \
    cp SAIGE/extdata/step2_SPAtests.R /usr/local/bin/step2_SPAtests.R && \
    chmod +x /usr/local/bin/step2_SPAtests.R && \
    rm -rf SAIGE

# Display SAIGE help
RUN createSparseGRM.R  --help && \
    step1_fitNULLGLMM.R --help && \
    step2_SPAtests.R --help
