# Running ORFold:

## Advanced run:

### Map the three methods along the genome of an organism

In the previous section we presented how to launch ORFold on one FASTA file of amino acid sequences. This is a simple example where the user has some sequences and wants to characterize their foldability potential. However, ORFold contains some more advanced options which are usefull in the case that these sequences are extracted by a genome wide analysis and the user posseses a GFF file of their genome localisation. In this case the user can get new GFF file(s) which will map the method(s) (HCA, IUPred or Tango) asked by the user (with the **-options** label) along the organism's genome.

The user can pass this GFF file (with the **-gff** label) and ORFold will generate new GFF file(s), identical with the initail one but which will contain the information of the method(s) asked with the **-options** label.  

To do so, the user has to type the following command:
```{}
orfold -fna sequences.fasta -options HIT -gff sequences.gff 
```

This command will generate the **sequences.tab** file already seen in the basic run section and additionally will generate three new files:

 1. sequences_HCA.gff
 2. sequences_IUPRED.gff
 3. sequences_TANGO.gff

These GFF files are identical as the one provided by the user with the difference that the sequences are colored based on their HCA score, disorder propensity and aggregation propensity respectivelly. Blue color indicates low values while red color indicate high values of each method. These GFF files can be visualized on genome visualisation tools (such as IGV) and can give a visual mapping of these three methods along the genome of an organism. 

### Multiple files handle:

ORFold can handle multiple input files (FASTA and GFF) and will try to associate them based on the following rules:

1. If the user provides the same number of FASTA and GFF files, ORFold will try to associate them based on their root name, nomatter the order of the files.

		orfold -fna sequences_Y.fasta sequences_X.fasta -options H -gff sequences_X.gff sequences_Y.gff
	
	In this case will associate:

	* sequence_Y.fasta with sequence_Y.gff
	* sequence_X.fasta with sequence_X.gff 
	
	&nbsp;
	

2. I If the user provides the same number of FASTA and GFF files and their root names are not identical, then ORFold will associate them based on the order the files are given.

		orfold -fna sequences_Y.fasta sequences_X.fasta -options H -gff sequences_A.gff sequences_B.gff

	In this case will associate:

	* sequence_Y.fasta with sequence_A.gff
	* sequence_X.fasta with sequence_B.gff

	&nbsp;

3. If the user provides multiple FASTA files and only one GFF file, then:

	* If the name of the GFF file matches with the name of one of the FASTA files then they are associated and the rest of the FASTA files stay with no association.
		
			orfold -fna sequences_Y.fasta sequences_X.fasta -options H -gff sequences_X.gff

		ORFold will associate:
	
		* sequences_X.fasta with sequnces_X.gff
		* sequences_Y.fasta with **nothing**
		
		&nbsp;

	* If the name of the GFF file does not match with any name of the FASTA files, then the GFF file will be associated with all the FASTA files considering that the FASTA files correspond to different subgroups of the same dataset.
			
			orfold -fna sequences_Y.fasta sequences_X.fasta -options H -gff sequences_B.gff

		ORFold will associate:

		* sequences_Y.fasta with sequences_B.gff
		* sequences_X.fasta with sequences_B.gff
		
		&nbsp;

4. If the user provides multiple FASTA and GFF files (but not the same number), then all the GFF files must have a FASTA file with the corresponding root name. Otherwise ORFold does not know how to associate the files and gives an ERROR message. 

		orfold -fna sequences_Y.fasta sequences_X.fasta sequences_Z.fasta -options H -gff sequences_Z.fasta sequences_Y.gff

	ORFold will associate:

	* sequences_Y.fasta with sequences_Y.gff
	* sequences_Z.fasta with sequences_Z.gff
	* sequences_X.fasta with nothing
	
	&nbsp;

		orfold -fna sequences_Y.fasta sequences_X.fasta sequences_Z.fasta -options H -gff sequences_B.fasta sequences_A.gff

	ORFold will give the following error message:
		
		Oups! You provided GFF file(s) which has no correspodance to FASTA

### Random sample of sequences 
Working with total genomes could eventually generate big amounts of sequences which will need an important excecution time for ORFold. If the user does not want to calculate the foldability of all the sequences but mostly prefers to have an indication based on a representative sample of sequences, he can create a random sample of the size of his preference. To do so, he must use the **-N** label and pass the size of sample he wishes to generate.

	orfold -fna sequences.fasta -options HIT -gff sequences.gff -N 3000

In this example, ORFold will calculate the three methods (HCA, IUPred and Tango) on a sample of 3000 sequences extracted randomly from the initial **sequences.fasta** file.    



