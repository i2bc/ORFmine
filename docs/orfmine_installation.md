## Installation

### Overview

ORFmine is a package that consists of different programs: ORFtrack, ORFribo, ORFold and ORFdate.
They can be used together or independently with the exception of ORFribo, which needs the ORFtrack's output file as input to be functional.


### Requirements

First of all, Docker<sup><a href="#references">1</a></sup> or Singularity<sup><a href="#references">2</a></sup> must be present in order to get the image containing all the scripts and tools needed for the analysis. 

Singularity might be prefered as it does not need super user rights, which is interesting for the use of ORFmine on a cluster or on a lab's computer. Singularity installation for UNIX-based OS is explained [here](https://singularity-tutorial.github.io/01-installation/).

Docker Engine is available on different OS like MacOS and Windows10 through Docker Desktop and as a static binary installation for a variety of Linux platforms (Docker Desktop for Linux is now available on some distributions). Everything is available [here](https://docs.docker.com/engine/install/).

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        For Windows, WSL2 and Ubuntu from Microsoft store applications usually are needed too.
    </p>
</div>

For reproducibility purposes, all programs and dependencies used by ORFmine are listed [here](./dependencies.md).


### Pull the ORFmine docker image

All the scripts and tools needed for the analysis were put and installed in a docker image stored on Dockerhub. You can pull and use it from a terminal as a docker image by writing :

``` bash
docker pull lopesi2bc/orfmine:latest
```

or as a singularity image :

```bash
# this will build a singularity image named orfmine_latest.sif that will be located in YOUR_PATH (to adapt)
singularity build YOUR_PATH/orfmine_latest.sif docker://lopesi2bc/orfmine:latest
```

This step might take about 10-20 minutes depending on your computer.

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        If you have any error, it might come from a permissions problem so you should try using these commands with sudo as prefix.  
    </p>
</div>


The quick installation is now complete! Please have a look [here](./orfmine_quickstart.md) to start your container.

<br>


## References

1. Merkel, Dirk. "Docker: lightweight linux containers for consistent development and deployment." Linux j 239.2 (2014)
2. Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459. Published 2017 May 11. doi:10.1371/journal.pone.0177459
