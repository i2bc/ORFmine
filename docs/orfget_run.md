## Extraction and writing of ORF sequences with ORFget

ORFget is a tool provided with ORFtrack that allows the user to extract
the protein and/or nulceotide sequences of specific subsets of ORFs 
according to their annotation categories
(see [here](./orftrack_annotation.md) 
for a description of all ORF categories). ORFget deals with annotation 
patterns, thereby allowing different levels of annotation in a 
very easy fashion.

ORFget has two principal options:

* ```-features_include```: list of motifs that will be used to define the 
  ORFs that will be included in the fasta 
  output. The sequences whose annotations include these patterns will 
  be retained in the output fasta file 
* ```-features_exclude```: list of motifs that will be used to define the 
  ORFs that will be excluded in the fasta 
  output. The sequences whose annotations include these patterns 
  will not be written in the output fasta file
  
The searched patterns can be specific (for a finer selection) or more general.<br><br>

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
For example, the motif <b>"nc"</b> which refers to all NonCoding ORFs appears in the features: <b>nc</b>_intergenic, <b>nc</b>_ovp_same-mRNA, <b>nc</b>_ovp_opp-mRNA and <b>nc</b>_ovp_same-tRNA.
As a result, the option <code>-feature_include nc</code> will keep all the four
feature categories. 
<br><br>
The option <code>-feature_include nc_ovp</code> will keep:
<ul>	
 <li><b>nc_ovp</b>_same-mRNA</li>
 <li><b>nc_ovp</b>_opp-mRNA</li>
 <li><b>nc_ovp</b>_same-tRNA</li>
</ul>

The option <code>-feature_include nc_ovp_same</code> will keep:
<ul>
 <li><b>nc_ovp_same</b>-mRNA</li>
 <li><b>nc_ovp_same</b>-tRNA</li>
</ul>

The option <code>-feature_include mRNA</code> will keep: 
<ul>
 <li>nc_ovp_same-<b>mRNA</b></li>
 <li>nc_ovp_opp-<b>mRNA</b></li>
</ul>

The option <code>-feature_exclude opp</code> will eliminate the nc_ovp_<b>opp</b>-mRNA and will keep:
<ul>
 <li>nc_intergenic</li>
 <li>nc_ovp_same-mRNA</li>
 <li>nc_ovp_same-tRNA</li>
etc... 
</p>
</div>
Here are presented some examples of selection of ORFs with ORFget.


### Extraction of the sequences of all the ORFs of a GFF file

The following command writes the amino acid sequences of all the ORFs 
 of the gff file annotated by ORFtrack. 

``` python
orfget -fna /database/genome.fasta -gff /database/mapping_orf_genome.gff
```
ORFget generates a fasta file containing all the corresponding amino acid
sequences. The output fasta file is written in the /database/ directory of the container.


<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
It can also handle gff files that were not generated by ORFtrack but in this case the user must be sure of the feature names to be indicated if using the -feature_include/exclude options. In this case, the features must correspond to those indicated in the 3rd column of the input gff file.
</p>
</div>


### Extraction of the sequences of all noncoding ORFs identified with ORFtrack

The following commands, each enable the user to write the 
amino acid sequences of all noncoding 
ORFs no matter their status (i.e. intergenic or overlapping)
(see [here](./orftrack_annotation.md) for a description of all ORF categories).

``` bash
orfget -fna /database/genome.fasta -gff /database/mapping_orf_genome.gff -features_include nc
```
or 
``` bash
orfget -fna /database/genome.fasta -gff /database/mapping_orf_genome.gff -features_include nc_intergenic nc_ovp
```
or
``` bash
orfget -fna /database/genome.fasta -gff /database/mapping_orf_genome.gff -features_exclude c_CDS
```

The output fasta file is written in the /database/ directory of the container and renamed based on the gff file rootname and the include features.

### Extraction of the sequences of a specific subset of ORFs according to their annotation

The following instruction writes the amino acid sequences of the ORFs
which overlap CDSs on the same, or on the opposite strand.

``` python
orfget -fna /database/genome.fasta -gff /database/mapping_orf_genome.gff -features_include nc_ovp_same-CDS nc_ovp_opp-CDS
```


Notice that using the argument "features_exclude" assumes that the selection 
operates on all genomic features except those that are excluded. 
Consequently, if the user wants to select all noncoding sequences
except those overlapping CDSs, mRNAs, tRNAs, and rRNAs, he/she must 
exclude the coding ORFs (c_CDS) as well. Otherwise, they will be
kept.


``` python
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_exclude c_CDS nc_same_ovp-tRNA nc_same_ovp-rRNA nc_opp_ovp-mRNA nc_opp_ovp-tRNA nc_opp_ovp-rRNA nc_opp_ovp-mRNA  
```
The output fasta file is written in the /database/ directory of the container and renamed based on the gff file rootname and the include features.

### Extraction of the sequences of a random subset of ORFs 

Sometimes, for computational time or storage reasons, the user does 
not want to deal with all the ORFs of a specific category. ORFget
can provide the user with a subset of N (to be defined by the user)
randomly selected ORFs from a specific ORF category. The last instruction
writes the sequences of 10000 randomly selected noncoding 
intergenic ORFs.


``` python
orfget -fna /database/genome.fasta -gff /database/mapping_orf_genome.gff -features_include nc_intergenic -n 10000
```

### Reconstruction of protein sequences
In addition, ORFget enables the reconstruction of all protein 
sequences of a genome (including all isoforms) according to their 
definition in the original gff file. The following instruction
writes all the resulting sequences in a fasta file. Please note that in this case, the input gff file is the initial one, not the one generated by ORFtrack which is ORF-centered and that contains C_CDS ORFs instead of the exact CDSs (i.e. ATG-STOP including isoforms) as indicated in the original gff file.


``` python
orfget -fna /database/genome.fasta -gff /database/genome.gff -features_include CDS
```

### Writing amino acid or nucleotide sequences
By default, ORFget will generate the amino acid sequences of the 
desired ORFs in a fasta file 
with the extension **.pfasta**. If the user wishes to generate the nucleotide
or even both nucleotide and amino acids sequences, he/she must use the 
option
```-type nucl``` and ```-type both```, respectively.


