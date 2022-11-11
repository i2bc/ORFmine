## ORFdate input files

ORFdate takes as input files:


* a distance phylogenetic tree in newick format, containing the species of interest (focal species) and its neighbors (see ORFmine/examples/database/orfdate_inputs/ for an example). All branches must be associated with their relative distance on which the age estimation will be based.
* a two columns csv file (i.e. columns separated by commas) with the path to the sequence fasta files of each species as first column, and the associated species names as indicated in the phylogenetic tree as second column (see ORFmine/examples/database/orfdate_inputs/ for an example). It, therefore, involves the fasta files of the complete proteomes (amino acid sequences) of the focal and its neighbors to be stored in the corresponding directories.
* the fasta files containing all the amino acid sequences of the focal and its neighboring species.


All these files must be stored directly or in a subdirectory of the /database/ directory of the container (see [here](./orfmine_quickstart.md#prepare-your-folders) for more details on the preparation of the folders of the container).
