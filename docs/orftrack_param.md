## ORFtrack Parameters

### **Mandatory**

- `--fna`  
  Nucleotide FASTA file of the genome whose ORFs are to annotate.

- `--gff`  
  GFF annotation file.

---

### **Optional**

- `-h, --help`  
  Shows the help and exits.

- `--codon-table`  
  Codon table ID to be used for the considered chromosomes.  
  The ID must be consistent with the NCBI codon table IDs (see  
  [NCBI Codon Tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).  
  *(Default: The standard codon table, ID = 1).*

- `--chr`  
  List of seqID to be treated by ORFtrack (i.e., column #1 of the GFF file, generally chromosome or contig ID).  
  IDs must be separated by a space:  
  `--chr NC_001148.4 NC_001139.3`  
  In this case, ORFtrack will not treat the entire genome but only the specified seqIDs.  
  **Recommendation**: Use this option for large genomes to distribute calculations across multiple CPUs.

- `--chr-exclude`  
  List of seqID(s) to exclude from processing.  
  *(Default: None).*

- `--types_only`  
  Genomic feature(s) considered for the annotation step (`CDS` is included by default).  
  If multiple genomic features are considered, they must be separated by spaces.  
  Noncoding ORFs are annotated as intergenic (if they do not overlap any feature) or overlapping (if they overlap a given genomic feature).  
  *(See the [ORFget section](./orfget_run.md) for more details).*

- `--types_except`  
  Genomic feature(s) excluded from the annotation step (`gene` and `exon` are excluded by default).  
  If multiple genomic features are excluded, they must be separated by spaces.  
  Noncoding ORFs will be annotated as intergenic or overlapping another feature not listed here.  
  *(See the [ORFget section](./orfget_run.md) for more details).*

- `--orf_len`  
  Minimal number of nucleotides between two consecutive STOP codons to define an ORF.  
  *(Default: 60 nucleotides).*  
  *(See the [ORF definition section](./orftrack_orfdef.md) for more details).*

- `--co_ovp`  
  Minimal fraction of the ORF length that overlaps a genomic feature to annotate the ORF as overlapping it.  
  *(Default: 70%).*  
  *(See the [Overlap section](./orftrack_overlap.md) for more details).*

- `-out`  
  Output directory.

- `--show-types`  
  Prints all genomic features annotated in the input GFF file.

- `--show-chrs`  
  Prints all the seqIDs (usually corresponding to chromosome or contig IDs) present in the input GFF file.

- `--ofasta`  
  Writes amino acid and nucleotide FASTA files for ORFs.

