# Running ORFold:

## Advanced run:

### Mapping of the fold potential and the disorder and aggregation propensities along the genome of an organism

In the previous section we presented how to launch ORFold on a set
of amino acid sequences stored in a fasta file. However, 
the originality of ORFold relies on the fact that the user can manually 
inspect the distribution of the properties estimated with ORFold (fold potential,
and disorder and aggregation propensities) along a genome of interest. 
In this case, the user must provide the genome annotation file (gff) along with 
the input fasta file. ORFold will return new gff files (one per studied property)
that contain for the ORFs provided in the input fasta file, their corresponding 
property scores (fold potential, disorder or aggregation propensities). The 
values are stored in the column #9 of the output gff files. The gff files can be subsequently
uploaded on a genome viewer such as IGV [1].

The input gff file must be given with the **-gff** option as follows:

```{}
orfold -faa /database/sequences.fasta -options HIT -gff /database/sequences.gff 
```

ORFold generates a **sequences.tab** file containing the fold potential, and the 
disorder and aggregation propensities of each sequence present in the input fasta file.
Additionally, ORFold produces three new gff files:

 1. sequences_HCA.gff
 2. sequences_IUPRED.gff
 3. sequences_TANGO.gff

The output gff files are identical to the one provided by the user except that for the sequences present
at the same time in the fasta and the gff files that were given as inputs, the column #9
is now replaced by the fold potential, or the disorder or aggregation propensities calculated by ORFold.
That way, the corresponding sequences can be colored according to these values on
a genome viewer, thereby enabling the visual inspection of these properties along
the input genome. Notice that on IGV, blue indicates low values (for all mapped properties) 
while red indicates high values. 

 ![HCA Scale](./img/mapping/Scale.png)<br>
<em>Figure 1: Color scale for the HCA score values 
</em>

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">

Notice that the ID 
of the sequences given in the fasta file (i.e. annotation after the ">" in 
the fasta file) must be strictly identical to those of the corresponding sequences 
in the gff file (i.e. ID indicated in the column #3 of the gff file).  
         
    </p>
</div>

### Dealing with multiple files at the same time:

ORFold can handle multiple input files (fasta and gff) and will associate 
the fasta and gff files according to the following rules:

1. If the user provides the same number of fasta and gff files, ORFold associates them 
   based on their root name, no matter the order of the files.

		orfold -faa /database/sequences_Y.fasta /database/sequences_X.fasta -options H -gff /database/sequences_X.gff /database/sequences_Y.gff
	
	In this case, ORFold associates:

	* sequence_Y.fasta with sequence_Y.gff
	* sequence_X.fasta with sequence_X.gff 
	
	&nbsp;
	

2. I If the user provides the same number of fasta and gff files,
   but their root names are not identical, ORFold associates them 
   according to the order of the files in the command line.

		orfold -faa /database/sequences_Y.fasta /database/sequences_X.fasta -options H -gff /database/sequences_A.gff /database/sequences_B.gff

	In this case, ORFold associates:

	* sequence_Y.fasta with sequence_A.gff
	* sequence_X.fasta with sequence_B.gff

	&nbsp;

3. If the user provides multiple fasta files and only one gff file, then:

	* If the name of the gff file matches with the name of one of the fasta files,
	  the two files are associated, while the other fasta files are not associated
	  to the input gff file.
		
			orfold -faa /database/sequences_Y.fasta /database/sequences_X.fasta -options H -gff /database/sequences_X.gff

		ORFold associates:
	
		* sequences_X.fasta with sequences_X.gff
		* sequences_Y.fasta with **nothing**
		
		&nbsp;

	* If the name of the gff file does not match with any of the names of 
	  the fasta files, then the gff file will be associated with all the fasta
	  files, considering that the fasta files correspond to different 
	  subgroups of the same dataset.
			
			orfold -faa /database/sequences_Y.fasta /database/sequences_X.fasta -options H -gff /database/sequences_B.gff

		ORFold associates:

		* sequences_Y.fasta with sequences_B.gff
		* sequences_X.fasta with sequences_B.gff
		
		&nbsp;

4. If the user provides multiple fasta and gff files (but not the same number), all gff files must have a corresponding fasta file with the same root name. 
   Otherwise, ORFold will give an ERROR message. 

		orfold -faa /database/sequences_Y.fasta /database/sequences_X.fasta /database/sequences_Z.fasta -options H -gff /database/sequences_Z.fasta /database/sequences_Y.gff

	ORFold associates:

	* sequences_Y.fasta with sequences_Y.gff
	* sequences_Z.fasta with sequences_Z.gff
	* sequences_X.fasta with nothing
	
	&nbsp;

		orfold -faa /database/sequences_Y.fasta /database/sequences_X.fasta /database/sequences_Z.fasta -options H -gff /database/sequences_B.fasta /database/sequences_A.gff

	ORFold will give the following error message:
		
		Oups! You provided GFF file(s) which has/have no correspondance to the input FASTA files

All these examples assume that the input files are stored in the /database/ directory of the container. All the outputs are written in the /workdir/orfold/ directory of the container.

###  Running ORFold on subsets of randomly selected sequences 
Working with complete genomes could generate big amounts of sequences 
which can dramatically increase the computational time of ORFold when dealing with
large genomes (especially if the estimation of the disorder and aggregation propensities are
activated). If the user does not want to treat all 
the sequences, he/she can create a random sample of the input sequences, large enough
 to have an estimation of the distribution of the studied properties of its 
dataset from a representative sample of the input sequences. To do so, 
the user must indicate the number of sequences that are to be randomly selected
with the **-N** option. For a representative dataset, we recommend selecting at least
10000 sequences.

	orfold -faa /database/sequences.fasta -options HIT -gff /database/sequences.gff -N 10000

In this example, ORFold will estimate the fold potential, and the disorder and 
aggregation propensities on a sample of 10000 sequences extracted randomly 
from the initial **sequences.fasta** file.    

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
	
	If the user works with more than one fasta file and wishes to create 
	random samples for all the input sequence files, he/she has to indicate in 
	the -N option the size for each input file explicitly (in the same order 
	as the inputs passed in the -fna option).
	<br>
```{}
orfold -faa /database/sequences_X.pfasta /database/sequences_Y.pfasta -options H -N 1500 3000
```	
Also, if the user wants to sample two subsets of same sizes, he/she has to indicate the subset sizes explicitly for each input
```{}
orfold -faa /database/sequences_X.pfasta /database/sequences_Y.pfasta -options H -N 1500 1500
```
	If the user whishes to calculate the fold potential of all the sequences 
	of one of the given inputs, he has to indicate it with the "all" flag (again with respect to
    the order of input files)
```bash
orfold -faa /database/sequences_X.pfasta /database/sequences_Y.pfasta -options H -N all 3000
```
In this case, ORFold will calculate the fold potential for <b>all</b> the sequences 
in the sequences_X file while will generate a random sample of 3000 sequences for the 
sequences_Y file.   
	
    </p>
</div>

References

1. Robinson JT, Thorvaldsdóttir H, Winckler W, et al (2011)
   Integrative genomics viewer. Nature biotechnology 29:24–26