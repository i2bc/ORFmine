# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).



##  version [1.0.0](https://github.com/i2bc/orfmine/releases) - 2022-11-25
### Added
- softwares.ini file to inform on absolute local path to tango and/or iupred2a


### Changed
- use of docker image through a handler process - commands to use docker image are the same as without docker (breaking)
- output directory options for all orfmine callable scripts
- rename output of orfplot (replace "output.png" by "outdir/filename.png")
- place -o argument option in orfget with -outdir & -outname in orfget (not required by default)


##  version [0.8.7](https://github.com/i2bc/orfmine/releases) - 2022-11-12
### Added
- docker image of ORFmine (breaking)

### Changed
- whole new structure dependent on docker container filesystem  (breaking)


# Changelog

All notable changes to this project will be documented in this file.

## [3.0.0] - 2024-12-26

### Added
- Added new arguments to `argparse` and updated the CLI options.
- Integrated a reporting script into the package.
- Added detailed log management for STAR.
- Added an option to stop Snakemake execution if FASTQ files fail the selection threshold.
- Included the `--intron_length` option, applicable only for STAR.
- Created a new version of `orfdate.py` using `subprocess` instead of `Bio.Application`.
- Integrated `MultiQC` at the beginning and after trimming.
- Added scripts for checking read lengths (25-35) before processing with Ribowaltz.
- Added new parameter management: replaced `softwares.ini` with `--path_iupred` and `--path_tango`.

### Changed
- Reorganized Snakemake rules to avoid execution order-related errors.
- Updated Dockerfile: installed dependencies via `setup.py` and `requirements.txt` before activating the Conda environment.
- Updated `ORFmine_env.yml`: removed the `pip` section and distributed its contents into `requirements.txt` and `setup.py`.
- Suppressed error messages related to threshold selection across multiple FASTQ files.
- Updated STAR command used for rRNA removal.
- Revised scripts for histogram generation and periodicity plots (`plotting_periodicity.py`).
- Adjustments for `ORFold` to run properly in Docker:
  - Added Biopython 1.68.
- Renamed ambiguous outputs in Ribowaltz (e.g., `selected length`, `{sample}.bam`).
- Renamed files and functions:
  - From `Mapping/Filter_Unwanted_Sequence` to `Mapping_Unwanted_Sequence_Filtering`.

### Fixed
- Resolved an issue in `bam2reads.py` by correctly adding `.constant` and `.loader` modules.
- Fixed the inability to run `orfold -h` in Docker.
- Added DendroPy (version 4.5.0) to `requirements.txt` and `setup.py` for ORFdate.
- Addressed errors caused by empty values in `config.yml`.

### Deprecated
- Removed `softwares.ini` in favor of new parameters.

### Removed
- Deleted unnecessary GFF rules (`02_gff_name_editing.smk`).
- Eliminated unused sections from `ORFmine_env.yml`.

### Notes
- **Updated Documentation:** Explained the file structure, provided example commands for datasets, and added input/output details.
- **Tests Performed:** 
  - Verified mapping tools (Bowtie2, STAR, Hisat2) across standalone Snakemake and package-integrated versions.
  - Checked results after integrating Ribowaltz scripts.
  - Tested the Singularity image.

