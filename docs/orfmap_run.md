## Annotation of ORFs with ORFmap

### Running ORFmap on a complete genome

The following instruction annotates all the possible ORFs of the
input genome.


``` bash
orfmap -fna genome.fasta -gff genome.gff
```
Depending on the size of the genome, ORFmap takes a few minutes to
several hours to annotate all the ORFs of a genome. It returns 
a new GFF file containing the annotation of all the identified ORFs 
(including coding and noncoding ORFs). 


### Running ORFmap on a single chromosome or a subset of chromosomes

ORFmap can be launched on a single chromosome (e.g. chromosome ID XXX) 
with the following 
instruction:


``` bash
orfmap -fna genome.fasta -gff genome.gff  -chr chr_ID XXXX
```
This can be very useful if the user wants to run ORFmap on several 
CPUs. Also, it can be launched on a subset of chromosomes as follows:


``` bash
orfmap -fna genome.fasta -gff genome.gff  -chr CHR1_ID CHR2_ID2 CHR3_ID3
```



