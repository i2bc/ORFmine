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


##### Download the latest ORFmine release
All the scripts and tools needed for the analysis were put and installed in a docker image stored on Dockerhub. You can pull and use it from a terminal as a docker image by writing :

``` bash
docker pull annelopes94/orfmine:v0.8.7
```

or as a singularity image :

```bash
# this will build a singularity image named orfmine_v0.8.7.sif that will be located in YOUR_PATH (to adapt)
singularity build YOUR_PATH/orfmine_v0.8.6.sif docker://annelopes94/orfmine:v0.8.7
```

This step might take about 10-20 minutes depending on your computer. 

If you have any error, it might come from a permissions problem so you should try using these commands with sudo as prefix.  

<br>

The quick installation is now complete! Please have a look [here](./orfmine_quickstart.md) to start your container.



<br><br>




## References

1. Merkel, Dirk. "Docker: lightweight linux containers for consistent development and deployment." Linux j 239.2 (2014)
2. Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459. Published 2017 May 11. doi:10.1371/journal.pone.0177459
3. Dosztanyi, Z., Csizmok, V., Tompa, P. & Simon, I. The pairwise energy content estimated from amino acid composition discriminates between folded and intrinsically unstructured proteins. Journal of molecular biology 347, 827–839 (2005).
4. Dosztányi, Z. Prediction of protein disorder based on IUPred. Protein Science 27, 331– 340 (2018).
5. Mészáros, B., Erdős, G. & Dosztányi, Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46, W329–W337 (2018).
6. Fernandez-Escamilla, A.-M., Rousseau, F., Schymkowitz, J. & Serrano, L. Prediction of sequence-dependent and mutational effects on the aggregation of peptides and proteins. Nature biotechnology 22, 1302–1306 (2004).
7. Linding, R., Schymkowitz, J., Rousseau, F., Diella, F. & Serrano, L. A comparative study of the relationship between protein structure and β-aggregation in globular and intrinsically disordered proteins. Journal of molecular biology 342, 345–353 (2004).
8. Rousseau, F., Schymkowitz, J. & Serrano, L. Protein aggregation and amyloidosis: confusion of the kinds? Current opinion in structural biology 16, 118–126 (2006).
