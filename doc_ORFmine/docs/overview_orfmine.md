## Introduction



ORFmine is an open-source package that aims at extracting, annotating,
and characterizing the fold potential and the structural properties of
all Open Reading Frames (ORF) encoded in a genome (including coding and
noncoding sequences). ORFmine consists of two independent programs, 
ORFmap and ORFold that can be used together or independently
(see here for an exmaple of application /usage_orfmine.md)).

ORFMap searches for all possible ORFs in the 6 frames of an input
genome, and annotate them according to a set of genetic features 
(e.g. noncoding intergenic ORFs, coding ORFs, noncoding ORFs 
that overlap with a specific genetic feature...). It provides
the user with a GFF file containing the annotations of all identified ORFs
that can be directly uploaded on a genome viewer for a visual inspection.
In addition, their amino acid OR nucleotide sequences can be extracted 
in a FASTA file (see complete documentation of ORFmap at XXXXX).

ORFold probes the fold potential and the disorder and aggregation 
propensities of a set of amino acid sequences.
The fold potential is estimated with the HCA method [REF], while the
disorder and aggregation propensities are calculated with IUPRED, and
TANGO respectively [REFiupred REFtango]. The specificity of ORFold lies
in the fact that the user can provide the amino acid sequences along 
with their corresponding annotations in a GFF file. In this
case, ORFold produces new GFF files, each containing for each annotated
sequence, their fold potential, disorder and aggregation propensities 
respectively, thereby enabling the manual inspection of these 
properties for all annotated ORFs/sequences in a genome viewer.



Used together, ORFmap and ORFold allow the characterization of the 
fold potential and structural properties of the potential peptides or proteins 
encoded in all ORFs of a genome regardless of their status (coding or
noncoding ORFs)



ORFold probes the fold potential and the disorder and aggregation 
propensities of a set of amino acid sequences provided in a FASTA file.
The fold potential is estimated with the HCA method [REF], while the
disorder and aggregation propensities are calculated with IUPRED, and
TANGO respectively [REFiupred REFtango]. The specificity of ORFold lies
in the fact that if the user provides the amino
acid sequences in a GFF file, ORFold it can handle a GFF file 
searches for 
all possible ORFs

Although HCA is very fast and can handle all ORFs of a small genome
in a few minutes, TANGO and IUPRED are much slower 
(X hours for all the ORFs of Escherichia coli str. K-12 substr. 
MG1655 on a personal computer). 
Consequently, the user can turn off the calculation of the disorder
and aggregation propensities. ORFold takes as input a FASTA file 
containing the amino acid sequences to treat. 
The corresponding GFF file can be given as well. In this case, the foldability potential will be stored/added/written in the GFF file. The latter can be loaded subsequently on a genome viewer such as IGV [REF], enabling the visualization of the distribution of the foldability potential on the whole genome/enabling a visual and manual analysis of the foldability potential of the whole genome. The output of ORFold is a table containing the foldability potential and/or the disorder and aggregation propensities of each input sequence and optionally a GFF file for mapping the foldability potential on the genome if the inputs include a GFF file/if the user provided/gave a GFF file as input/if a GFF file has been given as input/to visualize the distribution of the foldability potential of the different genome regions if the inputs include a GFF file/if the user provided/gave a GFF file as input/if a GFF file has been given as input. The program can handle several FASTA files at the same time and will generate as many outputs as given FASTA files (is it true?). Finally, ORFold can also provide the user with foldability plots which represent the distribution of the foldability of the input/entered sequences along with those of a dataset of globular proteins as reference. 




If the user is only interested
in the annotation and extraction of a genome's ORFs, whereas ORFold can estimate the foldability potential of any set of sequences no matter their annotation or genome localization/without using genomic information.


exploring
the potential of a complete genome  to give rise to novel with the 
extraction and annotation of all the possible ORF present in 
noncoding regions. The ORFmine package is available at xxx and 
