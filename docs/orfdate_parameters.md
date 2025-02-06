## ORFdate Parameters

### **Mandatory**

- `--target TARGET, -T`  
  Taxonomic name of the focal species.

- `--mapping , -M`  
  CSV file matching FASTA (column 1) and tree (column 2) names.

- `--tree TREE, -W`  
  Newick file for the phylogeny tree.

---

### **Optional**

- `-h, --help`  
  Show this help message and exit.

- `--cpus , -P`  
  Total number of CPUs that can be used for the task.  
  *(Default: 1).*

- `--evalue , -E`  
  BLASTp e-value threshold.  
  *(Default: 1e-3).*

- `--min-coverage , -C`  
  Minimum query coverage threshold.  
  *(Default: 0.70).*

- `--has-underscores , -S`  
  Whether underscores are kept when reading the tree or considered as spaces.  
  *(Default: False).*  
  - If `True`, underscores in species names in the tree and CSV file are treated explicitly and not replaced by spaces.  
  - If `False`, underscores in tree labels are replaced by spaces.  
  Useful when the species name contains underscores (e.g., `Scer_annotationV2`).

- `--keep-files, -K`  
  Keep intermediary computed files such as BLASTdb and BLAST output.

- `--blast BLAST, -B`  
  Whether to perform BLASTp or not.  
  `True` to perform BLASTp, `False` otherwise.  
  *(Default: True).*

- `--out OUT, -O`  
  Output directory.  
  *(Default: `.`).*

