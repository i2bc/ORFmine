## ORFfold Parameters

### **Mandatory**

- `--faa`  
  FASTA file containing the amino acid sequences to treat.

- `--options`  
  Indicates which properties are to be calculated.  
  `H` for estimating the fold potential with HCA,  
  `I` for the estimation of the disorder propensity with IUPred,  
  and `T` for the aggregation propensity with Tango.  
  Combinations of letters are accepted if the user wants to calculate several properties at the same time.  
  Example: `--options HIT` will estimate all three properties.  
  *(Default: H).*

- `--path_tango` 
   Path to tango executable file

- `--path_iupred`
   Path to iupred executable file
---

### **Optional**

- `-h, --help`  
  Shows this help message and exits.

- `--gff`  
  GFF annotation file. The ID (i.e., annotation) of the sequences given in the input FASTA file must match the ID label in the GFF file (column #3).  
  ORFold generates one GFF file per studied property (fold potential, disorder, and/or aggregation propensities). Each GFF file contains, for the sequences in the input FASTA file, the corresponding property values (fold potential, disorder, or aggregation propensities).  
  The values are stored in column #9 of the output GFF files, which can subsequently be uploaded into a genome viewer.  
  *(See [here](./Run_orfold_advanced.md) for examples of this option).*

- `--keep`  
  ORFold uses IUPred and Tango for predicting disorder and aggregation propensities.  
  By default, ORFold does not save the output files of these methods to save storage.  
  However, with the `-keep` option, the user can save the Tango output files.  
  Example: `-keep T` will save Tango output files in the `TANGO` directory.

- `--sample, -N`  
  Working with large genomes can generate very large files.  
  To address this, ORFold allows generating a sample of the initial dataset and performing calculations on this specific subset.  
  This option enables the user to specify the number of randomly selected sequences to be treated by ORFold.

-  `--out `   
  Output directory (default: '.').
  

