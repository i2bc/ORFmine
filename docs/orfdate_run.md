## Running ORFdate
ORFdate estimates the evolutionary ages of a set of ORFs (coding but also noncoding) based on phylostratigraphy. It takes as inputs:


* a phylogenetic distance tree in newick format, containing the species of interest (focal species) and its neighbors 
* a two column csv file (i.e. columns separated by commas) with the path to the sequence fasta files of each species as first column, and the associated species names as indicated in the phylogenetic tree as second column. It, therefore, involves the fasta files of the complete proteomes (amino acid sequences) of the focal and its neighbors to be stored in the corresponding directories
* the fasta files containing all the amino acid sequences of the focal and its neighboring species

All the inputs must be stored in the /database/ directory of the container, though a subdirectory can be created for more clarity (e.g. /database/orfdate_inputs/). The following command line will estimate the ages of all sequences present in the fasta file of the focal species (be they coding or noncoding) based on the phylogenetic tree provided by the user (an example of the different inputs can be found in ORFmine/examples/database/orfdate_inputs/). 



```
orfdate -focal focal_name -tree /database/orfdate_inputs/species_tree.nwk -names /database/orfdate_inputs/names.csv
```

with the -focal option corresponding to the name of the focal name as indicated in the tree; -tree option corresponding to the distance tree (newick format) of the focal and its neighbors and -names option corresponding to the two column csv file with the path/filenames of the fasta files and the corresponding species names. Please note that additional arguments including BLAST parameters can be provided. See [here](./orfdate_parameters.md) for more details.



