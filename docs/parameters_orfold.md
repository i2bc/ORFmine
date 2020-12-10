#Parameters of ORFold:
####-keep
ORFold uses IUPRED and TANGO for predicting the Disorder and Aggregation prone regions on each sequence.
By default, ORFold does NOT save the output files of these two methods for memory reasons. 
If the user wants to keep the IUPRED & TANGO output files he/she can activate this option through the label -keep. 
As for the -options label, TI or IT will save both files for IUPRED and TANGO while I or T will save only the IUPRED or TANGO (respectively) output files. The outputs are stored in the directories called IUPRED & TANGO, respectively.  

####-plot
By default, there are three reference datasets for evaluating the HCA score.
Disorder regions, Globular proteins and transmebrane helices. (Papadopoulos et al. 202?)
The user can project the foldability HCA scores of the sequences in the FASTA file(s) by activating the option -plot (type -plot True).

####-gff
This option is more advanced. 
The user can also give GFF file(s) which must correspond to the amino-acid sequences in the FASTA file(s). 
The name of each sequence must be identical with the "ID" label in the GFF file. 
This option gives a supplementary output which is a new GFF file of the sequences but this time colored based on their Score HCA (Foldability potential).
Blue for low Score and Red for high Score. 
This GFF can be visualized to adequate softwares (such as IGV) and map the foldability HCA score throughout the whole genome.  
This option can be perfectly used with the output FASTA and GFF files of the ORFmap tool for mapping the foldability potential of the non coding ORFs on the genome.


