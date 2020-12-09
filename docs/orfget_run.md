## Extraction and writing of ORF sequences with ORFget

ORFget is a tool provided with ORFmap that allows the user to extract
the protein or nulceotide sequences of specific subsets of ORFs 
according to their annotation categories. ORFget deals with regular 
expressions, thereby allowing different 
levels of annotation in a very easy fashion.

Here are presented some examples of selection of ORFs with ORFget.


### Extraction of the sequences of all the ORFs of a GFF file

The following command writes the protein sequences of all ORFs 
annotated in the input GFF file.


``` python
orfget   ALL SEQUENCES PROT

```
ORFget generates a FASTA file containing all the corresponding protein
sequences. 



### Extraction of the sequences of all noncoding ORFs identified with ORFmap

The following two commands, each enable the user to write the 
protein sequences of all noncoding 
ORFs no matter their status (i.e. intergenic or overlapping)
(see [here](./orfmap_annotation.md) for a description of all ORF categories).

``` python
orfget include intergenic  + nc_same_ovp + nc_opp_ovp

```
or 
``` python
orfget exclude c_CDS

```

### Extraction of the sequences of a specific subset of ORFs according to their annotation

The following instruction writes the protein sequences of the ORFs
which overlap with CDS on the same, or on the opposite strand.

``` python
orfget include  nc_same_ovp_CDS + nc_opp_ovp_CDS

```


Notice that using the argument "exclude" assumes that the selection 
operates on all genomic features except those that are excluded. 
Consequently, if the user wants to select all noncoding sequences
except those overlapping CDS, mRNAs, tRNAs, and rRNAs, he must 
exclude the coding ORFs (c_CDS) as well. Otherwise, they will be
kept.


``` python
orfget exclude C_CDS nc_same_ovp_tRNA + nc_same_ovp_rRNA +nc_opp_ovp_mRNA + C_CDS nc_opp_ovp_tRNA + nc_opp_ovp_rRNA 
+nc_opp_ovp_mRNA  

```

### Extraction of the sequences of a random subset of ORFs 

Sometimes, for computational time or storage reasons, the user does 
not want to deal with all the ORFs of a specific category. ORFget
can provide the user with a subset of N (to be defined by the user)
randomly selected ORFs from a specific ORF category. The last instruction
writes the sequences of 10000 randomly selected noncoding 
intergenic ORFs.


``` python
orfget 10000 intergenic ORFs

```
### Reconstruction of protein sequences
In addition, ORFget enables the reconstruction of all protein 
sequences of a genome (i.e. all isoforms) according to their 
definition in the original GFF file. The following instruction
writes all the resulting sequences in a FASTA file.


``` python
orfget build CDS
```
