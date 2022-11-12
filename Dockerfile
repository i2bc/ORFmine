FROM ubuntu:18.04

LABEL version="0.8.7"

ENV LC_ALL=C

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y wget git-all \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -u -p /bin/miniconda3/ \
    && rm Miniconda3-latest-Linux-x86_64.sh

ADD . /ORFmine/
RUN rm -r /ORFmine/docs

RUN chmod +x /ORFmine/orfribo/RiboDoc_BAM2Reads/RiboDoc_BAM2Reads.sh /ORFmine/install.sh

RUN mv /root/.bashrc /bin/ && ln -s /bin/.bashrc /root/.bashrc && \
    mv /root/.conda /bin/ && ln -s /bin/.conda /root/.conda && \
    mv /root/.profile /bin/ && ln -s /bin/.profile /root/.profile && \
    ln -s /root/ /bin/

ENV PATH=/bin/miniconda3/bin:$PATH
ENV PATH=/bin/miniconda3/envs/ORFmine_env/bin:$PATH
ENV PATH=$PATH:/bin

SHELL ["/bin/bash", "-c"]

RUN conda install mamba -n base -c conda-forge
RUN mamba env create -f /ORFmine/ORFmine_env.yml
RUN mamba init
RUN mamba init bash
RUN printf "mamba activate ORFmine_env\n" >> /bin/.bashrc \
    && mamba info --envs \
    && mamba list --name ORFmine_env

RUN R -e "library('devtools');devtools::install_github('LabTranslationalArchitectomics/riboWaltz@v1.2.0', dependencies = FALSE);"

RUN /ORFmine/install.sh ORFmine_env

RUN printf '#!/bin/bash\nbash /ORFmine/orfribo/RiboDoc_BAM2Reads/RiboDoc_BAM2Reads.sh "$@"' > /usr/bin/orfribo && \
    chmod +x /usr/bin/orfribo

WORKDIR /workdir/

CMD ["echo","Please use an interactive shell for the container"]
