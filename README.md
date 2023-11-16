# ORFmine

<div align="center">
  <img src="./docs/img/icons/ORFmine.png" width="80%"/>  
</div>

ORFmine is an open-source package that aims at extracting, annotating, and characterizing the sequence and structural properties of all Open Reading Frames (ORFs) of a genome, including coding as well as noncoding sequences, along with their translation activity. ORFmine consists of several independent programs that can be used together or independently:

- <i>**ORFtrack** searches for all possible ORFs longer than 60 nucleotides in the six frames of an input genome, and annotate them according to a set of genomic features.</i>
- <i>**ORFold** predicts the fold potential and the disorder and aggregation propensities of amino acid sequences.</i>
- <i>**ORFribo** probes the ORFs translation activity based on Ribosome Profiling data (Ribo-Seq).</i>
- <i>**ORFdate** estimates the ORFs evolutionary age based on phylostratigraphy information.</i>

More information can be found in the ORFmine [documentation](https://i2bc.github.io/ORFmine/).


## Requirements

ORFmine requires several dependencies and external softwares. To simplify installation, we offer a Docker image providing the complete environment required to use all of the ORFmine tools.

<details open>
<summary><h4>Minimal requirements for a container usage (recommended)</h4></summary>
To use the Docker image, you will need:

- Python >= 3.9
- ORFmine >= 2.0.0
- Docker or Singularity
</details>


<details>
<summary><h4>Minimal requirements for a local installation</h4></summary>
Alternatively, if you want to set up your environment for a local usage of ORFmine:

- Python >= 3.9
- blast >= 2.13
- bowtie2 == 2.5.0
- hisat2 == 2.2.1
- gffread == 0.12.7
- samtools == 1.16.1
- FastQC == 0.11.9
</details>


<details>
<summary><h4>Optional requirements</h4></summary>
Two other external softwares may be used for ORFold computations:

- Tango == 3.1 
- IUPred2A
</details>


## Recommendation

Before installing ORFmine, we strongly recommend to set up an Python isolated environment in order to avoid potential version conflicts between python libraries when working on different projects or different ORFmine versions.

Click in the section below for a short illustration on how to use an Python isolated environment.

<details style="margin-left: 32px">
<summary>How to use an isolated environment (recommended)</summary>
<br>
<p>
By using an isolated environment you will avoid potential version conflicts between python libraries when working on different projects. Some of the most popular tools to work with isolated python environments are [virtualenv](https://pypi.org/project/virtualenv/), [pyenv](https://pypi.org/project/pyenv/), [pipenv](https://pypi.org/project/pipenv/). 
</p>

Below is an example on how to use [virtualenv](https://pypi.org/project/virtualenv/).

#### 1. Install virtualenv
```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install virtualenv
python3 -m pip install virtualenv
```

#### 2. Create and activate an isolated environment
```bash
# create an isolated environment named 'orfmine_env' (to adapt)
virtualenv orfmine_env

# activate your isolated environment
source orfmine_env/bin/activate
```

Once activated, any python library you'll install using pip will be installed in this isolated environment, and python will only have access to these packages.

Once you're done working on your project, simply type `deactivate` to exit the environment.
</details>


## Installation

> :bell: **Note**
 The ORFmine package must be installed locally even if you plan to use the Docker image. This is because ORFmine includes a feature that simplifies the Docker usage, eliminating the need for complex volume mounting commands.


ORFmine can be accessed in different ways. Follow instructions described in option 1 or 2 if you're not interested in accessing/modifying the source code, otherwise prefer option 3. 

<details open>
<summary><h4>Option 1: from the archive (git not required)</h4></summary>

First download an archive of our latest release <a href="https://github.com/i2bc/ORFmine/releases/latest" target="_blank">here</a>.

```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install ORFmine vx.x.x
python3 -m pip install ORFmine-vx.x.x.zip # (or .tar.gz)
```
</details>


<details>
<summary><h4>Option 2: from the version control systems</h4></summary>

```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install ORFmine vx.x.x
python -m pip install -e git+https://github.com/i2bc/ORFmine.git@v2.0.0#egg=orfmine
```
</details>

<details>
<summary><h4>Option 3: from this project repository</h4></summary>

```bash
# clone ORFmine on your machine
git clone https://github.com/i2bc/ORFmine.git

# go in the ORFmine/ directory
cd ORFmine

# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install ORFmine in edition mode (useful for a development process)
python3 -m pip install -e .
```
</details>


## Documentation

All details about ORFtrack, ORFold, ORFdate, ORFribo and their usage can be found in the ORFmine [documentation](https://i2bc.github.io/ORFmine/).


## Licence

The ORFmine project is under the MIT licence. Please check [here](https://github.com/i2bc/ORFmine/blob/ORFmine_complete/LICENSE.md) for more details.


## Citing

If you use ORFmine for your research, please cite:
> Papadopoulos, C., Chevrollier, N., Lopes, A. Exploring the peptide potential of genomes. Meth. Mol. Biol. (2022)

> Papadopoulos, C., Arbes, H., Chevrollier, N., Blanchet, S., Cornu, D., Roginski, P., Rabier, C., Atia, S., Lespinet, O., Namy, O., Lopes, A. The Ribosome Profiling landscape of yeast reveals a high diversity in pervasive translation. bioRxiv (2023)

