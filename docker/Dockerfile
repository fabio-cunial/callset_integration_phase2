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
ARG gcloud_version=517.0.0
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# AWS
RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip \
    && unzip awscli-exe-linux-x86_64.zip \
    && ./aws/install \
    && aws --version

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

# SEQKIT
ARG seqkit_version=2.10.0
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz \
    && tar -xzvf seqkit_linux_amd64.tar.gz \
    && rm -f seqkit_linux_amd64.tar.gz \
    && chmod +x seqkit \
    && ./seqkit --help

# Latest version of kanpig
ARG kanpig_version=1.0.2
RUN wget https://github.com/ACEnglish/kanpig/releases/download/v${kanpig_version}/kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && tar xzf kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && rm -f kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && mv ./kanpig-v${kanpig_version}-x86_64-unknown-linux-musl/kanpig . \
    && chmod +x ./kanpig \
    && ./kanpig --help

# TRUVARI
ARG truvari_version=5.2.0
RUN pip3 install truvari==${truvari_version} \
    && truvari --help

# BEDTOOLS
ARG bedtools_version=2.31.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools.static \
    && mv bedtools.static bedtools \
    && chmod +x bedtools \
    && mv bedtools /usr/local/bin \
    && bedtools --help

# PIGZ
RUN apt-get install -y --no-install-recommends pigz

# SEQKIT
ARG seqkit_version=2.8.2
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz \
    && tar -xzvf seqkit_linux_amd64.tar.gz \
    && rm -f seqkit_linux_amd64.tar.gz \
    && chmod +x seqkit \
    && ./seqkit --help
    
# MINIMAP2. We need v2.15 because it is the one wrapped by pbmm2 v1.10.0 that 
# is used in AoU Phase 2 production.
ARG minimap2_version=2.15
RUN wget https://github.com/lh3/minimap2/releases/download/v${minimap2_version}/minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && tar xjf minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && rm -f minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && mv ./minimap2-${minimap2_version}_x64-linux/ ./minimap2/ \
    && ./minimap2/minimap2 -h && echo 0 || echo 1

# VCFDIST
ARG vcfdist_version=2.6.3
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
RUN wget https://github.com/TimD1/vcfdist/archive/refs/tags/v${vcfdist_version}.tar.gz \
    && tar xzf v${vcfdist_version}.tar.gz \
    && rm -f v${vcfdist_version}.tar.gz \
    && cd vcfdist-${vcfdist_version}/src \
    && make -j 16 \
    && make install \
    && vcfdist --help

# AD HOC SCRIPTS
COPY *.java .
RUN javac *.java
COPY *.py .
