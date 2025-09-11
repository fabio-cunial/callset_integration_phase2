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

# Kanpig 1.1.0 supports --squish
ARG kanpig_version=1.1.0
RUN wget https://github.com/ACEnglish/kanpig/releases/download/v${kanpig_version}/kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && tar xzf kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && rm -f kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && mv ./kanpig-v${kanpig_version}-x86_64-unknown-linux-musl/kanpig . \
    && chmod +x ./kanpig \
    && ./kanpig --help

# Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs > script.sh \
    && chmod +x script.sh \
    && ./script.sh -y
    
# Custom Kanpig
ARG commit_id=262b89805a3c32b2a1b4beb99b3414b080a2b3ac
RUN mkdir kanpig-custom \
    && cd kanpig-custom \
    && git clone https://github.com/fabio-cunial/kanpig.git \
    && cd kanpig \
    && git checkout ${commit_id} \
    && $HOME/.cargo/bin/cargo build --release \
    && mv ./target/release/kanpig ../../kanpig-tmp \
    && cd ../.. \
    && rm -rf kanpig-custom/ \
    && mv kanpig-tmp kanpig-custom \
    && ./kanpig-custom gt --help

# TRUVARI
ARG truvari_version=5.2.0
RUN pip3 install truvari==${truvari_version} \
    && truvari --help

# Beta-binomial
COPY ./pyro-kanpig-environment.yml .
COPY ./genotype-beta-binomial-mixture.py .
RUN conda env create -f ./pyro-kanpig-environment.yml

# AD HOC SCRIPTS
COPY *.java .
RUN javac *.java
