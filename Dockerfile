FROM ubuntu:20.04

################## METADATA ######################
LABEL about.summary="Collection of tools for the BSF atacseq_pipeline"
LABEL about.home="https://www.biomedical-sequencing.org"
LABEL about.tags="Epigenetics, ATAC-seq"
LABEL about.maintainer="Bekir Erguener @ CeMM/BSF"
LABEL version="0.2"
################## MAINTAINER ####################

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y --no-install-recommends \
	g++ \
	gcc \
	make \
	autoconf \
	wget \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	parallel \
	pigz \
	alien \
	zip \
	unzip \
	python3.8 \
	python3.8-dev \
	python3-pip \
	bedtools \
	bowtie2 \
	&& apt-get -y clean  && apt-get -y autoclean  && apt-get -y autoremove

RUN pip3 install deeptools==3.5.0 \
		macs2==2.2.7.1

# Compile BWA
WORKDIR /root/source
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
	&& tar -xf bwa-0.7.17.tar.bz2
WORKDIR /root/source/bwa-0.7.17
RUN make

# Compile bcftools
WORKDIR /root/source
RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 \
	&& tar -xf bcftools-1.11.tar.bz2 
WORKDIR /root/source/bcftools-1.11
RUN ./configure \
	&& make

# Compile samtools
WORKDIR /root/source
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
	&& tar -xf samtools-1.11.tar.bz2
WORKDIR /root/source/samtools-1.11
RUN ./configure  --without-curses \
	&& make

# Compile samblaster
WORKDIR /root/source
RUN wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.24/samblaster-v.0.1.24.tar.gz \
	&& tar -xf samblaster-v.0.1.24.tar.gz
WORKDIR /root/source/samblaster-v.0.1.24
RUN make

# Compile fastp
WORKDIR /root/source
RUN wget https://github.com/OpenGene/fastp/archive/v0.20.1.tar.gz \
	&& tar -xf v0.20.1.tar.gz
WORKDIR /root/source/fastp-0.20.1
RUN make

RUN cp /root/source/bwa-0.7.17/bwa \
	/root/source/bcftools-1.11/bcftools \
	/root/source/samtools-1.11/samtools \
	/root/source/samblaster-v.0.1.24/samblaster \
	/root/source/fastp-0.20.1/fastp /usr/local/bin/

WORKDIR /homer
RUN wget http://homer.ucsd.edu/homer/configureHomer.pl \
	&& perl configureHomer.pl -install \	
	&& perl configureHomer.pl -install hg38 mm10
ENV PATH /homer/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends \
	git

RUN pip3 install multiqc

RUN git clone https://github.com/berguner/atacseq_pipeline.git && \
	python3 atacseq_pipeline/multiqc_atacseq/setup.py install

CMD ["/bin/bash"]

