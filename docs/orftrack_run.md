## Annotation of ORFs with ORFtrack

### Running ORFtrack on a complete genome

The following instruction annotates all the possible ORFs of the
input genome.


``` bash
orftrack -fna /database/genome.fasta -gff /database/genome.gff
```

Depending on the size of the genome, ORFtrack takes a few minutes to
several hours to annotate all the ORFs of a genome. It returns 
a new gff file containing the annotation of all the identified ORFs 
(including coding and noncoding ORFs). The output gff file is stored in the /database/ folder of the container and named mapping_orf_[genome].gff. It also provides a summary.log file that summarizes for each chromosome/contig all the ORF categories that have been detected.


### Running ORFtrack on a single chromosome or a subset of chromosomes

ORFtrack can be launched on a single seqID (usually chromosome or contig indicated in the first 
column of the input gff)(e.g. chromosome seqID: XXX) 
with the following instruction:


``` bash
orftrack -fna /database/genome.fasta -gff /database/genome.gff  -chr chr_ID_XXXX
```
This can be very useful if the user wants to run ORFtrack on several 
CPUs. Also, it can be launched on a subset of seqIDs as follows:


``` bash
orftrack -fna /database/genome.fasta -gff /database/genome.gff  -chr seqID1 seqID2 seqIDx
```



