## ORFget Parameters

### **Mandatory**

- `--fna`  
  Nucleotide FASTA file of the genome whose ORFs are to be extracted.

- `--gff`  
  GFF annotation file.

- `--outdir`  
  Root name of the output FASTA.

---

### **Optional**

- `-h, --help`  
  Shows this help message and exits.

- `--outname`  
  Output basename of the generated FASTA file(s).  
  *(default: Automatically built from the GFF filename and requested features).*

- `--features_include`  
  Genomic features or annotation patterns to be searched for ORF selection.  
  ORFs whose annotations match the entered patterns or genomic features are written in the output FASTA file.  
  *(default: All genomic features present in the input GFF, consequently all annotated ORFs will be written in the output GFF file).*  
  *(See the [ORFget section](./orfget_run.md) for more details).*

- `--features_exclude`  
  Genomic features or annotation patterns to be searched for ORF exclusion.  
  ORFs whose annotations match the entered patterns or genomic features are not written in the output FASTA file.  
  *(default: None, consequently all annotated ORFs will be written in the output GFF file).*  
  *(See the [ORFget section](./orfget_run.md) for more details).*

- `--chr_exclude`  
  To exclude a (or a subset of) seqID(s) from the ORF extraction.  
  ORFs belonging to these seqIDs (usually chromosomes or contigs) will not be written in the output FASTA file.  
  *(default: None).*

- `--N`  
  Number of ORFs randomly selected for writing, according to their genomic feature or pattern.

- `--type`  
  Sequences type of output (`prot` *(default)* / `nucl` / `both`).

- `--check`  
  DO NOT write sequences that DO NOT finish with a STOP codon.

- `--elongate`  
  Number of nucleotides to elongate towards 5' & 3' UTR.

