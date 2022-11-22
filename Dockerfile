#################
# 1st build stage
#################

FROM python:3.9-slim-buster as stage_1

RUN apt-get update && \
    apt-get install -y build-essential  && \
    apt-get install -y libc6-dev && \
    apt-get install -y git && \
    apt-get install -y libncurses5-dev && \
    apt-get install -y default-jre && \
    apt-get install -y libbz2-dev && \
    apt-get install -y liblzma-dev && \
    apt-get install -y zlib1g-dev

# create a user and go in home
RUN adduser orfmine
WORKDIR /home/orfmine

# create a virtual environment
RUN python -m venv .env

# create env variables for virtual env
ENV VIRTUAL_ENV /home/orfmine/.env/bin

# Make sure we use the virtualenv:
ENV PATH ${VIRTUAL_ENV}:${PATH} 

# copy orfribo required software binaries in VIRTUAL_ENV
COPY softwares_dependencies/bin/ ${VIRTUAL_ENV}

# install samtools
WORKDIR /home/orfmine/samtools
COPY softwares_dependencies/samtools-1.16.1 ./
RUN ./configure --prefix="/home/orfmine/.env" && make && make install



#################
# 2nd build stage
#################
FROM stage_1 as stage_2

RUN apt-get update  && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libssl-dev && \
    apt-get install -y libfontconfig1-dev  && \
    apt-get install -y libharfbuzz-dev && \
    apt-get install -y libfribidi-dev  && \
    apt-get install -y libfreetype6-dev  && \
    apt-get install -y libpng-dev && \
    apt-get install -y libtiff5-dev  && \
    apt-get install -y libjpeg-dev  && \
    apt-get install -y r-base

RUN R -e "install.packages('devtools', dependencies = TRUE);library('devtools');devtools::install_github('LabTranslationalArchitectomics/riboWaltz@v1.2.0', dependencies = TRUE);"

WORKDIR /home/orfmine


# USER orfmine

#################
# 3rd build stage
#################
FROM stage_2

# create a user
# RUN adduser orfmine
# COPY --from=env_riboseq /home/orfmine /home/orfmine

# Make sure we use the virtualenv:
# ENV PATH="/home/orfmine/.env/bin:$PATH"

# go in home
WORKDIR /home/orfmine

# add ORFmine packages and setup.py
COPY packages ./packages
COPY setup.py ./

# install ORFmine python libraries & dependencies 
RUN pip3 install -e .

COPY softwares.ini ./

# make current workdir content owns by user
RUN chown -R orfmine:orfmine ./

# edit .bashrc
RUN sed -i "s/#alias/alias/g" .bashrc

# log as user himself
USER orfmine