## Aims and general description

ORFMap scans a given genome in the 6 frames, and searches for 
all possible ORFs longer than a given size (default: 60 nucleotides
including the 3' STOP codon). It annotates them according to a set of genetic features (e.g. noncoding intergenic,
coding, noncoding and overlapping with a specific genetic feature - see
the ORF annotation section for more details ![ORF annotation](./orf_def_orfmap.md)). 
ORFmap takes as inputs a FASTA file containing the nucleotide
sequences of all chromosomes or contigs and their corresponding 
annotations in a GFF file. The program returns a new GFF file that contains all
identified ORFs. In addition, the amino acid sequences of 
all annotated ORFs or specific subsets of ORFs (i.e. only noncoding intergenic ORFs for example)
can be extracted and written in a FASTA file with ORFget, a tool 
provided with ORFmap. 
