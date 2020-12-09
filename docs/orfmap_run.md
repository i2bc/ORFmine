## Annotation of ORFs with ORFmap

### Running ORFmap on a complete genome

The following instruction annotates all the possible ORFs of the
input genome.


``` python
orfmap   orfmap -fna genome.fasta -gff genome.gff

```
Depending on the size of the genome, ORFmap takes a few minutes to
several hours to annotate all the ORFs of a genome. It returns 
a new GFF file containing the annotation of all the identified ORFs 
(including coding and noncoding ORFs). 


### Running ORFmap on a single chromosome or a subset of chromosomes

ORFmap can be launched on a single chromosome with the following 
instruction:


``` python
orfmap   orfmap -fna genome.fasta -gff genome.gff  CHR XXXXX

```
This can be very useful if the user wants to run ORFmap on several 
CPUs. Also, it can be launched on a subset of chromosomes as follows:


``` python
orfmap   orfmap -fna genome.fasta -gff genome.gff  CHR1 CHR2 CHR3

```



