# Overview

OMICS studies attribute a new role to the noncoding genome in the production of novel peptides. The widespread transcription
of noncoding regions and the pervasive translation of the resulting RNAs offer a vast reservoir of novel peptides to the organisms.

ORFmine<sup><a href="#references">1</a>, <a href="#references">2</a></sup> is an open-source package that aims at extracting, annotating, and characterizing the sequence and structural properties of
all Open Reading Frames (ORFs) of a genome (including coding and noncoding sequences) along with their translation activity. ORFmine consists of several independent programs that can be used together or independently:
- ORFtrack
- ORFold
- ORFribo
- ORFdate

## Built with
- python 3.6
- miniconda 3
- pyHCA <sup><a href="#references">1</a></sup>
- R 
- bash
- Docker

All programs and dependencies are listed [here](docs/dependencies.md). 


# Getting started

## Prerequisites

- [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://singularity-tutorial.github.io/01-installation/)
- [IUPred](https://iupred2a.elte.hu/download_new) <sup><a href="#references">2, 3, 4</a></sup>  (optional)
- [Tango](http://tango.crg.es) <sup><a href="#references">5, 6, 7</a></sup> (optional)


## Installation

Simply pull the ORFmine image from Dockerhub.

For docker:
```bash
# pull the ORFmine docker image from Dockerhub
docker pull annelopes94/orfmine:v0.8.7
```

For singularity:
```bash
# this will build a sif image named orfmine_v0.8.6.sif that will be located in YOUR_PATH (to adapt)
singularity build YOUR_PATH/orfmine_v0.8.6.sif docker://annelopes94/orfmine:v0.8.7
```

If you have any error, it might come from a permissions problem so you should try using these commands with sudo as prefix.  


# Usage

For usage examples, please check the [Quick start](https://i2bc.github.io/ORFmine/orfmine_quickstart.html) section of our documentation page.


# Documentation

Our full documentation is accessible [here](https://i2bc.github.io/ORFmine/)

# Issues

If you have suggestions to improve ORFmine or face technical issues, please post an issue [here](https://github.com/i2bc/ORFmine/issues)


# Contact

Anne Lopes - anne.lopes@i2bc.paris-saclay.fr


# Citing

Yet to come.


# Licence

The ORFmine project is under the MIT licence. Please check [here](LICENSE.md) for more details.


# References

1. Bitard-Feildel, T. & Callebaut, I. HCAtk and pyHCA: A Toolkit and Python API for the Hydrophobic Cluster Analysis of Protein Sequences. bioRxiv 249995 (2018).
2. Dosztanyi, Z., Csizmok, V., Tompa, P. & Simon, I. The pairwise energy content estimated from amino acid composition discriminates between folded and intrinsically unstructured proteins. Journal of molecular biology 347, 827–839 (2005).
3. Dosztányi, Z. Prediction of protein disorder based on IUPred. Protein Science 27, 331– 340 (2018).
4. Mészáros, B., Erdős, G. & Dosztányi, Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46, W329–W337 (2018).
5. Fernandez-Escamilla, A.-M., Rousseau, F., Schymkowitz, J. & Serrano, L. Prediction of sequence-dependent and mutational effects on the aggregation of peptides and proteins. Nature biotechnology 22, 1302–1306 (2004).
6. Linding, R., Schymkowitz, J., Rousseau, F., Diella, F. & Serrano, L. A comparative study of the relationship between protein structure and β-aggregation in globular and intrinsically disordered proteins. Journal of molecular biology 342, 345–353 (2004). 
7. Rousseau, F., Schymkowitz, J. & Serrano, L. Protein aggregation and amyloidosis: confusion of the kinds? Current opinion in structural biology 16, 118–126 (2006).


`singularity shell --pwd /workdir/ -B workdir:/workdir/ -B database:/database/ -B fastq:/fastq/ -B ~/tango:/ORFmine/orfold_v1/orfold/softwares/ ../../ORFmine/orfmine_test.sif`