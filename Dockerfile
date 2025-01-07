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


#################
# 2nd build stage - R 
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
    

#################
# 3rd build stage - Conda environment setup + STAR installation
#################
FROM stage_2 as stage_3

# Install Miniconda
WORKDIR /tmp
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Set path to conda
ENV PATH /opt/conda/bin:$PATH

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz && \
    tar -xzvf 2.7.10a.tar.gz && \
    rm 2.7.10a.tar.gz && \
    cd STAR-2.7.10a/source && \
    make STAR && \
    mv STAR /usr/local/bin/

# Clean up STAR build files
RUN rm -rf /tmp/STAR-2.7.10a

# Copy environment.yml to the container
COPY ORFmine_env.yml /tmp/ORFmine_env.yml

# Create Conda environment
RUN conda env create -f /tmp/ORFmine_env.yml && \
    conda clean -afy


#################
# 4th build stage - orfmine python dependencies
#################
FROM stage_3 as stage_4

# create a user and go in home
RUN adduser orfuser

# go in /home
WORKDIR /home/orfuser/orfmine

# create a virtual environment
RUN python3.10 -m venv env-orfmine
ENV VIRTUAL_ENV /home/orfuser/orfmine/env-orfmine/bin

# Make sure we use the virtualenv
ENV PATH ${VIRTUAL_ENV}:${PATH}

# Copy the requirements file into the container
COPY requirements.txt .

# Install the Python dependencies
RUN pip3 install --no-cache-dir numpy
RUN pip3 install --no-cache-dir -r requirements.txt


#################
# final stage - orfmine package
#################
FROM stage_4

# add ORFmine main package and setup.py
COPY orfmine ./orfmine
COPY setup.py ./

# create ini file for optional softwares related to orfold (iupred & tango)
#RUN printf "[EXTERNAL_SOFTWARE]\niupred = \"/opt/iupred2a\"\ntango = \"/opt/tango\"\n" > ./softwares.ini

# install ORFmine python libraries & dependencies 
RUN pip3 install -e .



# Activate Conda environment by default
ENV CONDA_DEFAULT_ENV=ORFmine_env
ENV PATH="/opt/conda/envs/ORFmine_env/bin:$PATH"


# create /inputs and /outputs directories with relevant user permissions
RUN mkdir /input /output && \
    chown orfuser:orfuser /input && \ 
    chown orfuser:orfuser /output && \
    chmod 755 /input /output

WORKDIR /input

# log as user himself
USER orfuser
