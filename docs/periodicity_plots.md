# CDS Periodicity Plot Utility

## Overview 
This script is an additional tool for **ORFribo** that visualizes peridicity in ribosome profiling data. It generates bar plots showing in-frame and out-of-frame ribosome occupancy across coding sequences (CDS).

## Usage 

### **1. Individual Mode (Single File)**  
To generate a periodicity plot from a single `.tab` file:

```bash
python metagenome_plotting.py --individual --tab path/to/file.tab --cds_name CDS001
```

Parameters:

`--individual` : Enables individual mode.

`--tab` : Path to the **.tab** file.

`--cds_name` : Name of the CDS to analyze.
 


### 2. Processing multiple CDS from a File

If you have multiple *.tab* and multiple CDS names in a text file, use: 

```bash 
python metagenome_plotting.py --pooled --directory path/to/directory --cds_file cds_list.txt 
```

Parameters:

`--cds_file` : Path to a text file containaing one CDS name per line.
`--pooled` : Enables pooled mode to process multiple files.
`--directory`:  Path to the directory containing multiple **.tab** files.


## Output 

The script generates **bar plots** where:  

- **Orange bars**: In-frame reads (P0).  
- **Blue bars**: Out-of-frame reads (P1 & P2).  

The plot is saved as a `.png` file with the **CDS name**.
