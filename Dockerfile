#################
# 1st build stage - ubuntu based image
#################
FROM ubuntu:22.04 as stage_1

# Set DEBIAN_FRONTEND to noninteractive
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN apt-get update && apt-get install -y \
    build-essential python3.10 python3-pip python3-venv \
    libc6-dev git libncurses5-dev default-jre \
    libbz2-dev liblzma-dev zlib1g-dev wget vim-tiny \
    libxml2-dev libcurl4-openssl-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

RUN apt-get clean && rm -rf /var/lib/apt/lists/*


# RUN apt-get update && \
#     apt-get install -y build-essential python3.10 python3-pip python3-venv && \
#     apt-get install -y libc6-dev git libncurses5-dev default-jre && \
#     apt-get install -y libbz2-dev liblzma-dev zlib1g-dev wget vim-tiny && \
#     apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev && \
#     apt-get install -y libfontconfig1-dev libharfbuzz-dev libfribidi-dev && \
#     apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev



#################
# 2nd build stage - R & ribowaltz
#################
FROM stage_1 as stage_2


# Add R repository
RUN apt-get update && \
    apt-get install -y software-properties-common dirmngr --no-install-recommends && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && \
    apt-get install -y r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install ribowaltz
RUN R -e "install.packages('devtools', dependencies = TRUE); \
         library('devtools'); \
         devtools::install_github('LabTranslationalArchitectomics/riboWaltz@v1.2.0', dependencies = TRUE);"


#################
# 3rd build stage - orfmine external software dependencies (bowtie2, gffread, hisat2, samtools, fastQC)
#################
FROM stage_2 as stage_3

# copy software binaries (bowtie2, gffread, hisat2, blastp, makeblastdb) to /usr/local/bin
COPY softwares_dependencies/bin/ /usr/local/bin

# install samtools
WORKDIR /opt/samtools
COPY softwares_dependencies/samtools-1.16.1 .
RUN ./configure --prefix="/usr/local" && make all all-htslib && make install install-htslib

# install fastqc
WORKDIR /opt/FastQC
COPY softwares_dependencies/FastQC/ .
RUN ln -s /opt/FastQC/fastqc /usr/local/bin


#################
# 3rd build stage - orfmine python dependencies
#################
FROM stage_3 as stage_4

# create a user and go in home
RUN adduser orfuser

# go in /home
WORKDIR /home/orfuser/orfmine

# create a virtual environment
RUN python3 -m venv env-orfmine
ENV VIRTUAL_ENV /home/orfuser/orfmine/env-orfmine/bin

# # Make sure we use the virtualenv:
ENV PATH ${VIRTUAL_ENV}:${PATH}

# Copy the requirements file into the container
COPY requirements.txt .

# Install the Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt


#################
# final stage - orfmine package
#################
FROM stage_4

# add ORFmine main package and setup.py
COPY orfmine ./orfmine
COPY setup.py ./

# create ini file for optional softwares related to orfold (iupred & tango)
RUN printf "[EXTERNAL_SOFTWARE]\niupred = \"/opt/iupred2a\"\ntango = \"/opt/tango\"\n" > ./softwares.ini

# install ORFmine python libraries & dependencies 
RUN pip3 install -e .

# # make current workdir content owns by user
# RUN chown -R orfuser:orfuser ./


# create /inputs and /outputs directries with relevant user permissions
RUN mkdir /input /output && \
    chown orfuser:orfuser /input && \ 
    chown orfuser:orfuser /output && \
    chmod 755 /input /output

WORKDIR /input

# log as user himself
USER orfuser