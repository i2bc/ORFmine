## Extraction and writing of ORF sequences with ORFget

ORFget is a tool provided with ORFmap that allows the user to extract
the protein and/or nulceotide sequences of specific subsets of ORFs 
according to their annotation categories. ORFget deals with regular 
expressions, thereby allowing different levels of annotation in a 
very easy fashion.

ORFget has two principal options:

* **-features_include** : List of motifs to match at GFF features column and keep the sequence  
* **-features_exclude** : List of motifs to match at GFF features column and eliminate the sequence

The motifs can be explicit (for detailed selection) or more implicit (for general selection).<br><br>
For example, the motif **"nc"** appears in the features: **nc**_intergenic, **nc**_ovp_same-mRNA, **nc**_ovp_opp-mRNA and **nc**_ovp_same-tRNA.
As a result the option ```-feature_include nc``` will keep all the four feature categories. 
<br><br>
The option ```-feature_include nc_ovp``` will keep:
	
* **nc_ovp**_same-mRNA
* **nc_ovp**_opp-mRNA
* **nc_ovp**_same-tRNA

The option ```-feature_include nc_ovp_same``` will keep:

* **nc_ovp_same**-mRNA
* **nc_ovp_same**-tRNA

The option ```-feature_include mRNA``` will keep: 

* nc_ovp_same-**mRNA**
* nc_ovp_opp-**mRNA**

etc... 

Here are presented some examples of selection of ORFs with ORFget.


### Extraction of the sequences of all the ORFs of a GFF file

The following command writes the protein sequences of all ORFs 
annotated in the input GFF file.


``` python
orfget -fna genome.fasta -gff mapping_orf_genome.gff
```
ORFget generates a FASTA file containing all the corresponding protein
sequences. 



### Extraction of the sequences of all noncoding ORFs identified with ORFmap

The following commands, each enable the user to write the 
protein sequences of all noncoding 
ORFs no matter their status (i.e. intergenic or overlapping)
(see [here](./orfmap_annotation.md) for a description of all ORF categories).

``` bash
orfget -fna genome.fasta -gff mapping_orf_genome -features_include nc
```
or 
``` bash
orfget -fna genome.fasta -gff mapping_orf_genome -features_include nc_intergenic nc_ovp
```
or
``` bash
orfget -fna genome.fasta -gff mapping_orf_genome -features_exclude c_CDS
```

### Extraction of the sequences of a specific subset of ORFs according to their annotation

The following instruction writes the protein sequences of the ORFs
which overlap with CDS on the same, or on the opposite strand.

``` bash
orfget -fna genome.fasta -gff mapping_orf_genome -features_include nc_ovp_same_CDS nc_ovp_opp_CDS
```


Notice that using the argument "features_exclude" assumes that the selection 
operates on all genomic features except those that are excluded. 
Consequently, if the user wants to select all noncoding sequences
except those overlapping CDS, mRNAs, tRNAs, and rRNAs, he must 
exclude the coding ORFs (c_CDS) as well. Otherwise, they will be
kept.


``` bash
orfget -fna genome.fasta -gff mapping_orf_genome -features_exclude c_CDS nc_same_ovp_tRNA nc_same_ovp_rRNA nc_opp_ovp_mRNA nc_opp_ovp_tRNA nc_opp_ovp_rRNA nc_opp_ovp_mRNA  
```

### Extraction of the sequences of a random subset of ORFs 

Sometimes, for computational time or storage reasons, the user does 
not want to deal with all the ORFs of a specific category. ORFget
can provide the user with a subset of N (to be defined by the user)
randomly selected ORFs from a specific ORF category. The last instruction
writes the sequences of 10000 randomly selected noncoding 
intergenic ORFs.


``` python
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_include nc_intergenic -n 10000
```

### Reconstruction of protein sequences
In addition, ORFget enables the reconstruction of all protein 
sequences of a genome (i.e. all isoforms) according to their 
definition in the original GFF file. The following instruction
writes all the resulting sequences in a FASTA file.


``` python
orfget -fna genome.fasta -gff genome.gff -features_include CDS
```
