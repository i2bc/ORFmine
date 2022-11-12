## Installation


### 1. Overview
ORFmine is a package that consists of different programs: ORFtrack, ORFribo, ORFold and ORFdate.
They can be used together or independently with the exception of ORFribo, which needs the ORFtrack's output file as input to be functional.


### 2. Quick Installation (without Tango and IUPred)

##### Install Singularity or Docker
First of all, Docker[1] or Singularity[2] must be present in order to get the image containing all the scripts and tools needed for the analysis. 

Singularity might be prefered as it does not need super user rights, which is interesting for the use of ORFmine on a cluster or on a lab's computer. Singularity installation is explained [here](https://singularity-tutorial.github.io/01-installation/).

Docker Engine is available on different OS like MacOS and Windows10 through Docker Desktop and as a static binary installation for a variety of Linux platforms (Docker Desktop for Linux is now available on some distributions). Everything is available [here](https://docs.docker.com/engine/install/).

All programs and dependencies used by ORFmine are listed [here](./dependencies.md).

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        For Windows, WSL2 and Ubuntu from Microsoft store applications usually are needed too.
    </p>
</div>

##### Download the latest ORFmine release
All the scripts and tools needed for the analysis were put and installed in a docker image stored on Dockerhub. You can pull and use it from a terminal as a docker image by writing :

``` bash
docker pull DOCKERHUB_REPOSITORY_ORFMINE-name_of_repo/orfmine
```

or as a singularity image :

```bash
singularity build /path/to/your/image.sif docker://name_of_repo/orfmine
```

This step might take about 10-20 minutes depending on your computer. 

If you have any error, it might come from a permissions problem so you should try using these commands with sudo as prefix.  


<a name="general_install"></a>

If you are not interested in the calculation of the disorder and/or aggregation propensities with IUPred and/or Tango in ORFold, nothing more is required for the installation and you can skip the rest of this page (please note that HCA is included in the docker and you will be able to predict the foldability of your sequence(s) of interest). The quick installation is now complete! Please have a look [here](./orfmine_quickstart.md) to start your container.


### 3. Full installation (with Tango and IUPred)

If you want to include IUPred and Tango in ORFold, you have to download them on your own and build the docker/singularity image but do not worry, you will be guided. As usual with docker, you might need the admin privileges to use it.

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
    The calculation of the disorder or aggregation propensities are both optional and complementary to the HCA score. As a result, IUPred and Tango tools are not mandatory for ORFold to be operational. In addition, they are not necessarily coupled together. ORFold will properly be installed without them or even with only one of them.
    </p>
</div>

First of all, you have to contact the IUPred and/or Tango developers through the respective links to have access to their programs as these two softwares (links to [IUPred](https://iupred2a.elte.hu/download_new)[3][4][5] and [Tango](http://tango.crg.es)[6][7][8]).

Then, all scripts and information about virtual environments for ORFmine can be downloaded [here](https://github.com/i2bc/ORFmine) with :
``` bash
mkdir /path/to/directory/for/ORFmine/;
cd /path/to/directory/for/ORFmine/;
git clone https://github.com/i2bc/ORFmine
```

Once ORFmine is downloaded and you have access to the IUPred and Tango scripts, you must place those in a specific folder to be able to use them with ORFold.
From you directory made for ORFmine :

* Move the IUPred source code and data (provided by the developer):
    ``` bash
	mv iupred2a.py ORFmine/orfold_v1/orfold/softwares/
    ```
* Move Tango source code:
	* For linux:
        ``` bash
		mv tango_x86_64_release ORFmine/orfold_v1/orfold/softwares/
        ```
    * For MacOS:
        ``` bash
		mv tango2_3_1 ORFmine/orfold_v1/orfold/softwares/
        ```
    * For windows:
        ``` bash
		mv Tango.exe ORFmine/orfold_v1/orfold/softwares/
        ```

Still from you directory made for ORFmine, you can now build the docker image including IUPred and Tango (this takes several dozen of minutes) :
``` bash
docker build -f Dockerfile --tag orfmine:latest .
```

At this step,  the image is built and if you do not need to use singularity you can launch a container with docker as explained in the exemple [here](./orfmine_quickstart.md).

In case you prefer the use of singularity, you have to make the *.sif* image in order to launch it :
``` bash
singularity build orfmine_latest.sif docker-daemon://repo/orfmine:latest;
```

You may need to change the permissions for the use of your *.sif* image. If you do not know exactly how, you can do :
``` bash
sudo chmod ogu+rwx orfmine_latest.sif
```

Well done ! The full installation including TANGO and IUPRED is now complete !
You can launch a container with singularity as explained in the exemple [here](./orfmine_quickstart.md).


<a name="launch_install"></a>


<br><br><br>
#### References

1. Merkel, Dirk. "Docker: lightweight linux containers for consistent development and deployment." Linux j 239.2 (2014)
2. Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459. Published 2017 May 11. doi:10.1371/journal.pone.0177459
3. Dosztanyi, Z., Csizmok, V., Tompa, P. & Simon, I. The pairwise energy content estimated from amino acid composition discriminates between folded and intrinsically unstructured proteins. Journal of molecular biology 347, 827–839 (2005).
4. Dosztányi, Z. Prediction of protein disorder based on IUPred. Protein Science 27, 331– 340 (2018).
5. Mészáros, B., Erdős, G. & Dosztányi, Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46, W329–W337 (2018).
6. Fernandez-Escamilla, A.-M., Rousseau, F., Schymkowitz, J. & Serrano, L. Prediction of sequence-dependent and mutational effects on the aggregation of peptides and proteins. Nature biotechnology 22, 1302–1306 (2004).
7. Linding, R., Schymkowitz, J., Rousseau, F., Diella, F. & Serrano, L. A comparative study of the relationship between protein structure and β-aggregation in globular and intrinsically disordered proteins. Journal of molecular biology 342, 345–353 (2004).
8. Rousseau, F., Schymkowitz, J. & Serrano, L. Protein aggregation and amyloidosis: confusion of the kinds? Current opinion in structural biology 16, 118–126 (2006).
