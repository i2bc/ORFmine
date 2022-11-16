## ORFribo inputs

ORFribo takes as inputs several files that must be stored either in the /database/ or /fastq/ folders of the container (see [here](./orfmine_quickstart.md#prepare-your-folders) for more details on the preparation of the folders of the container).

#### The following files must be stored in the /database/ directory of the container

 * a fasta file containing the nucleotide sequences
of the chromosomes or contigs of a genome
 * the original annotation of the genome in a gff file (see 
[the GFF3 documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
for more details on the gff3 format).
* the gff file containing the ORFs to be analyzed. We strongly recommend using the gff file [generated by ORFtrack](./orfmine_quickstart.md#annotation-and-extraction-of-orfs-with-orftrack) that contains the ORF annotation of the input genome according to [ORFtrack definition](./orftrack_annotation.md) and parameters (ORF [length](./orftrack_orfdef.md) and [overlap](./orftrack_overlap.md)). That said, the user can provide his/her own ORF annotation file. In this case, the ORF category of each annotated ORF (even if there is only one category, it must be named explicitly) must be indicated in the 3rd column of the gff file (see the ORFtrack output as example in ORFmine/examples/database/mapping_orf_Scer.gff)
* A fasta file with the sequences you do not want to treat and that you want to remove from the mapping, usually rRNA sequences (Ex : Scer_rRNA.fa)


#### The following files must be stored in the fastq directory of the container
* Ribosome Profiling data (i.e. the fastq(.gz) files) *(the /fastq/ folder must only contain the fastq files without anything else inside)*




See [here](./orfmine_quickstart.md#prepare-your-folders) for more details on the preparation of the folders. An example of inputs can be found in the ORFmine/examples/ directory.