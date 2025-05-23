# ORFmine

<div align="center">
  <img src="./docs/img/icons/ORFmine.png" width="80%"/>  
</div>

**ORFmine** is an open-source package designed to extract, annotate, and characterize the sequence and structural properties of all Open Reading Frames (ORFs) of a genome, including coding and noncoding sequences, along with their translation activity.

---

## Key Features

ORFmine includes several independent programs that can be used together or separately:

- **ORFtrack**: Searches for all possible ORFs (>60 nucleotides) in the six reading frames of a genome and annotates them based on genomic features.
- **ORFold**: Predicts the folding potential, disorder, and aggregation propensities of amino acid sequences.
- **ORFribo**: Analyzes ORFs' translation activity using Ribosome Profiling data (Ribo-Seq).
- **ORFdate**: Estimates the evolutionary age of ORFs using phylostratigraphy information.

More information is available in the [official documentation](https://i2bc.github.io/ORFmine/).

### Requirements

To get started, you will need:

- Python 3.9 (or 3.10)
- pip
- Either Docker **or** Singularity (Apptainer)

---

### Python Setup 

Check if Python 3.9 or 3.10 is already available on your system.

- **Already installed ?** Go to [Installation Options](#installation-options)

- **Not installed ?** follow one of the options below to create isolated python environment. 


#### Option A – Using virtualenv

```bash
   python3.9 -m pip install --upgrade pip
   python3.9 -m pip install virtualenv

   # Create environment
   virtualenv env-3.9
   
   # Activate python environment 
   source env-3.9/bin/activate

   # To deactivate the environment:
   deactivate
```

#### Option B – Using Conda

```bash
   # Create conda environment
   conda create --name env_python3.9 python=3.9

   # To activate the environement: 
   conda activate env_python3.9 
```

---

### Installation Options

ORFmine can be installed in several ways. Choose the option that best suits your needs: 
It will install the last version of orfmine 

#### Option 1 : From Pypi (Recommanded): 

```
    pip install orfmine 
```


#### Option 2: From a Local Repository

1. Clone the ORFmine repository:  

```
    git clone https://github.com/i2bc/ORFmine.git
```
2. Navigate to the cloned directory:  

```
    cd ORFmine
```

3. Install ORFmine in editable mode:  

```
    python3 -m pip install --upgrade pip
    python3 -m pip install -e .
```

---

**ORFmine relies on a suite of external tools in addition to its core Python scripts.
The commands above only install the Python components and a few dependencies. To ensure full functionality, including access to all required third-party software, you must complete the installation using the containerized environment provided (via Docker or Singularity). This guarantees a fully configured, portable, and reproducible setup.[Docker and Singularity build](#docker-and-singularity-build)**

---

**Note:**  Regardless of the installation method you choose, you must install the following tool to use ORFold:

```
   pip install git+https://github.com/T-B-F/pyHCA.git
```



## Docker and Singularity build

For containerized environments, ORFmine supports Docker and Singularity.
Make sure to have docker or singularity installed in your machine before building the image. 
For Docker, make sure you have root permissions. 

- **Docker**:  

```
    docker pull fadwa06/orfmine:3.0.1
```

- **Singularity (or Apptainer)**:  

```
   singularity build orfmine_latest.sif docker://fadwa06/orfmine:3.0.1
```

---

## Documentation

For detailed installation instructions, usage examples, and pipeline configurations, visit the [official documentation](https://i2bc.github.io/ORFmine/).

---

## License and Citation

- **License**: ORFmine is distributed under the MIT License.
- **Citation**: If you use ORFmine in your research, please cite the following works:

> Papadopoulos, C., Chevrollier, N., Lopes, A. Exploring the peptide potential of genomes. Meth. Mol. Biol. (2022).  
> Papadopoulos, C., et al. The Ribosome Profiling landscape of yeast reveals a high diversity in pervasive translation. bioRxiv (2023).
