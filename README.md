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

---

## Installation

### Requirements

To use ORFmine, the following versions are recommended:

- **Python** >= 3.9
- **ORFmine** >= 2.0.0
- **Docker** or **Singularity** (for containerized usage)

We recommend using an isolated Python environment to avoid version conflicts between libraries. See the section below for details.

---

### Setting Up an Isolated Python Environment (Recommended)

To create an isolated environment:

```bash
python3.9 -m pip install --upgrade pip
python3.9 -m pip install virtualenv
virtualenv orfmine_env
source orfmine_env/bin/activate
```

To deactivate the environment:

```bash
deactivate
```

---

### Installation Options

> **Note**: ORFmine must be installed locally even if you plan to use the Docker image to simplify its usage.

#### 1. Install via GitHub Releases

Download the latest release from [GitHub](https://github.com/i2bc/ORFmine/releases/latest):

```bash
python3.9 -m pip install ORFmine-vx.x.x.zip
```

#### 2. Install from the GitHub Repository

Install directly from the repository:

```bash
python3.9 -m pip install -e git+https://github.com/i2bc/ORFmine.git@v3.0.0#egg=orfmine
```

#### 3. Local Installation (for source code modification)

Clone the repository and install locally:

```bash
git clone https://github.com/i2bc/ORFmine.git
cd ORFmine
python3.9 -m pip install -e .
```

---

## Containerized Usage (Docker or Singularity)

To simplify setup, ORFmine is fully compatible with containerized environments.

- **Docker**:

```bash
$package_name $args --docker
```

- **Singularity**:

```bash
$package_name $args --singularity
```

---

## Documentation

For detailed installation instructions, usage examples, and pipeline configurations, visit the [full documentation](https://orfmine-docs-link.com).

---

## License and Citation

- **License**: ORFmine is distributed under the MIT License.
- **Citation**: If you use ORFmine in your research, please cite the following works:

> Papadopoulos, C., Chevrollier, N., Lopes, A. Exploring the peptide potential of genomes. Meth. Mol. Biol. (2022).  
> Papadopoulos, C., et al. The Ribosome Profiling landscape of yeast reveals a high diversity in pervasive translation. bioRxiv (2023).
