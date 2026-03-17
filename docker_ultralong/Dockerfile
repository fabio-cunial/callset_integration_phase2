FROM continuumio/miniconda3:24.1.2-0
MAINTAINER Fabio Cunial
ARG work_dir=/callset_integration
WORKDIR ${work_dir}


# OS
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        build-essential \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        python2 \
        python-dev \
        bsdmainutils \
        default-jre \
        default-jdk \
        unzip \
        python3-pycurl \
        python3-edlib 
        

# GSUTIL
ARG gcloud_version=544.0.0
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
ARG htslib_version=1.21
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
ARG samtools_version=1.21
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
ARG bcftools_version=1.21
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && bcftools --help

# BEDTOOLS
ARG bedtools_version=2.31.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools.static \
    && mv bedtools.static bedtools \
    && chmod +x bedtools \
    && mv bedtools /usr/local/bin \
    && bedtools --help

# TRUVARI
ARG truvari_version=5.4.0
RUN pip install truvari==${truvari_version} \
    && truvari --help

# SNIFFLES
# Remark: the latest version is currently 2.7.3, but it throws errors on our
# ultralong VCFs. E.g.:
#
# Sniffles2 Error: Error parsing input VCF: Line 1:
# 'pysam.libcbcf.VariantRecord' object has no attribute 'strip' (Fatal error,
# exiting.)
#
# The following is currently the latest release that doesn't throw errors.
ARG sniffles_version=2.7.3
RUN python -m venv /opt/sniffles_env
RUN /opt/sniffles_env/bin/pip install sniffles==${sniffles_version}

# CUTEFC
ARG cutefc_version=1.0.2
RUN python -m venv /opt/cutefc_env
RUN wget https://github.com/Meltpinkg/cuteFC/archive/refs/tags/cuteFC-v${cutefc_version}.tar.gz \
    && tar -xzf cuteFC-v${cutefc_version}.tar.gz \
    && rm -f cuteFC-v${cutefc_version}.tar.gz \
    && cd cuteFC-cuteFC-v${cutefc_version} \
    && /opt/cutefc_env/bin/python setup.py install

# LRCALLER
ARG lrcaller_version=1.0
RUN wget https://github.com/DecodeGenetics/LRcaller/releases/download/v${lrcaller_version}/lrcaller-${lrcaller_version}_linux-x86-64.tgz \
    && tar -xzvf lrcaller-${lrcaller_version}_linux-x86-64.tgz \
    && chmod +x lrcaller-${lrcaller_version}/bin/* \
    && mv lrcaller-${lrcaller_version}/bin/* /usr/local/bin \
    && lrcaller --help || echo "LRcaller returns an error code when called"

# AD HOC SCRIPTS
RUN pip install pandas numpy matplotlib seaborn cyvcf2 pysam scipy 
COPY *.java .
RUN javac *.java
COPY *.r .
COPY *.py .
