# getting bas image

FROM ubuntu:20.04

# versions -- bcl2fastq2 has reached end of life so version is hard-coded
ARG DEBIAN_FRONTEND=noninteractive
ARG STAR_VERSION=2.7.10b
ARG CELL_RANGER_VERSION=7.2.0
ARG BCL_CONVERT_VERSION=4.1.5-2.el7
ARG FASTQC_VERSION=0.12.1

# ===================================
# Set maintainer
LABEL maintainer="Andre Kurlovs <andre.kurlovs@sanofi.com>" \
      description="ToBridge portion of CellBridge Pipeline"
LABEL version="1.1.0"

RUN apt-get update && \
        apt-get install -y \
        apt-utils \
        ed \
        default-jre \
        less \
        locales \
        vim-tiny \
        wget \
        ca-certificates \
        apt-transport-https \
        curl \
        cmake \
        git \
        python3.9 \
        python3-pip \
        bzip2 \
        gzip \
        zip \
        tar \
        vim \
        nano \
        acl \
        software-properties-common \
        fonts-liberation \
        dirmngr \
        gcc \
        alien \
        libssl-dev \
        libcurl4-openssl-dev \
        libhdf5-dev \
        libxml2-dev \
        libnlopt-dev \
        libicu-dev \
        libjpeg62 \
        libgeos-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/*
 
# get R       
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    r-base \
    r-base-core
       
RUN pip3 install \
  setuptools \
  numpy \
  pandas \
  scipy \
  matplotlib \
  scikit-learn \
  RSeQC \
  bedparse \
  multiqc \
  pysam \
  psutil

RUN mkdir -p /opt/tobridge
WORKDIR /opt/tobridge

COPY cellranger-${CELL_RANGER_VERSION}.tar.gz /opt/tobridge/
RUN tar xzvf cellranger-${CELL_RANGER_VERSION}.tar.gz

COPY bcl2fastq2-v2-20-0-linux-x86-64.zip /opt/tobridge/
RUN unzip bcl2fastq2-v2-20-0-linux-x86-64.zip
RUN alien -i bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm

# genomes - get the human and mouse from 10x, custom will be provided on their own
RUN mkdir -p /opt/tobridge/genomes

WORKDIR /opt/tobridge/

# get fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip
RUN unzip fastqc_v${FASTQC_VERSION}.zip

# install bcl convert
COPY bcl-convert-${BCL_CONVERT_VERSION}.x86_64.rpm /opt/tobridge/
RUN alien -i bcl-convert-${BCL_CONVERT_VERSION}.x86_64.rpm

# get star
RUN wget https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz
RUN tar xzvf ${STAR_VERSION}.tar.gz
WORKDIR /opt/tobridge/STAR-${STAR_VERSION}/source
RUN make STAR
WORKDIR /opt/tobridge/

RUN cp -r cellranger-${CELL_RANGER_VERSION}/lib/python/cellranger/barcodes barcodes
COPY chemistry.csv /opt/tobridge/chemistry.csv
RUN gunzip barcodes/*.gz

COPY Base.py /opt/tobridge/
RUN chmod +x /opt/tobridge/Base.py
RUN ln /opt/tobridge/Base.py /usr/local/bin/tobridge

#this part is for converting indices from 10x to bcl convert friendly
#can be skipped if not needed
RUN mkdir -p /opt/tobridge/more_python
ADD index_csvs /opt/tobridge/index_csvs
COPY base_functions.py /opt/tobridge/more_python/base_functions.py
COPY convert_sample_sheet.py /opt/tobridge/more_python/convert_sample_sheet.py
ENV PYTHONPATH="/opt/tobridge/more_python:${PYTHONPATH}"

# The VOLUME instruction and the -v option to docker run tell docker to store 
# files in a directory on the host instead of in the container’s file system.
VOLUME /data
WORKDIR /data
