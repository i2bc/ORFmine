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

### Option 1: From an Archive (No Git Required)
1. Download the latest release archive from [here](https://github.com/i2bc/ORFmine/releases/latest).
2. Install ORFmine:  

```
    python3.9 -m pip install --upgrade pip
    python3.9 -m pip install ORFmine-vx.x.x.zip
```

### Option 2: From Version Control (Git)
Install directly from the GitHub repository:  

```
   python3.9 -m pip install --upgrade pip
   python3.9 -m pip install -e git+https://github.com/i2bc/ORFmine.git@v2.0.0#egg=orfmine
```

### Option 3: From a Local Repository
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

## Docker and Singularity Usage

For containerized environments, ORFmine supports Docker and Singularity:
For Docker, make sure you have root permissions 

- **Docker**:  

```
    $package_name $args --docker
```

- **Singularity**:  

```
    $package_name $args --singularity
```

---

## Conda Environment Usage (Not recommanded)

If you are using conda, you can create a Conda environment using the `ORFmine_env.yml` file:  

```
conda env create -f ORFmine_env.yml
```

Activate the environment:  
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

