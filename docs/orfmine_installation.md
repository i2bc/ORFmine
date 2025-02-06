# ORFmine Installation Guide

## Requirements

Before installing ORFmine, it is strongly recommended to set up an isolated Python environment to avoid potential library version conflicts across different projects.

### Using an Isolated Python Environment (Recommended)

Setting up an isolated Python environment (python >= 3.9) prevents library version conflicts. Here’s how to create and use one with `virtualenv`:

1. **Install virtualenv**:  

```
    python3.9 -m pip install --upgrade pip
    python3.9 -m pip install virtualenv
```

2. **Create and activate an isolated environment**:  

```
    virtualenv orfmine_env
    source orfmine_env/bin/activate
```

   To deactivate the environment:
```
   deactivate
```

Alternatively, ORFmine provides a Docker image for a fully configured environment.

---

## Installation Options

ORFmine can be installed in several ways. Choose the option that best suits your needs: 
It will install the last version of orfmine 

### Option 1 : From Pypi (Recommanded): 

```
    pip install orfmine 
```


### Option 2: From a Local Repository

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

**Note:** If you need to use ORFold, you must install the following tool:

```
   pip install git+https://github.com/T-B-F/pyHCA.git
```

## Docker and Singularity Usage

For containerized environments, ORFmine supports Docker and Singularity.
For Docker, make sure you have root permissions. 
During the first execution of the tool's modules, it will take a bit more time to retrieve and configure the Docker/Singularity image.

- **Docker**:  

```
    $package_name $args --docker
```

- **Singularity**:  

```
    $package_name $args --singularity
```


## Conda Usage (Not Recommanded)

If you are not familiar with Docker or Singularity, you can create a Conda environment.

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
4. Create Conda environement from yml file: 

```
conda env create -f ORFmine_env.yml
```

5. Activate the environment:  

```
conda activate ORFmine_env
```


---

## License

ORFmine is licensed under the MIT License. For more details, see the [LICENSE file](https://github.com/i2bc/ORFmine/blob/ORFmine_complete/LICENSE.md).

---

## Citation

If you use ORFmine in your research, please cite the following papers:

> Papadopoulos, C., Chevrollier, N., Lopes, A. Exploring the peptide potential of genomes. Meth. Mol. Biol. (2022)  
> Papadopoulos, C., Arbes, H., Chevrollier, N., Blanchet, S., Cornu, D., Roginski, P., Rabier, C., Atia, S., Lespinet, O., Namy, O., Lopes, A. The Ribosome Profiling landscape of yeast reveals a high diversity in pervasive translation. bioRxiv (2023)

