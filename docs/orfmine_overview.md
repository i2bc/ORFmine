## Overview of ORFmine



Recent studies attribute a new role to the noncoding genome in
the production of novel peptides. The widespread transcription
of noncoding regions and the pervasive translation of the resulting
RNAs offer a vast reservoir of novel peptides to the organisms.


ORFmine is an open-source package that aims at extracting, annotating,
and characterizing the fold potential and the structural properties of
all Open Reading Frames (ORF) encoded in a genome (including coding and
noncoding sequences). ORFmine consists of two independent programs, 
ORFmap and ORFold that can be used together or independently
(see [here](./orfmine_quickstart.md) for an example of 
application).

ORFMap searches for all possible ORFs longer than 60 nucleotides in the six frames of an input
genome, and annotate them according to a set of genomic features 
(e.g. noncoding intergenic ORFs, coding ORFs, noncoding ORFs 
that overlap with a specific genomic feature...). It provides
the user with a GFF file containing the annotations of all identified ORFs
that can be directly uploaded on a genome viewer for a visual inspection.
In addition, their amino acid and/or nucleotide sequences can be extracted 
in a FASTA file (for more details, see the complete 
documentation of ORFmap).

ORFold probes the fold potential and the disorder and aggregation 
propensities of a set of amino acid sequences.
The fold potential is estimated with the HCA method [REF], while the
disorder and aggregation propensities are calculated with IUPRred, and
TANGO respectively [REFiupred REFtango]. The specificity of ORFold lies
in the fact that the user can provide the amino acid sequences along 
with their corresponding annotations in a GFF file. In this
case, ORFold produces new GFF files, each containing for each annotated
sequence, their fold potential, disorder and aggregation propensities 
respectively, thereby enabling the manual inspection of these 
properties for all genomic features in a genome viewer
(for more details, see the complete 
documentation of ORFold).

