# ORFMap
ORFMap - A tool aimed at scanning a genome for stop-codons delimited sequences (ORFs) and annotating them.

## Summary
* <p><a href="#descr">Description</a></p>
* <p><a href="#installation">Installation</a></p>
* <p><a href="#usage_descr">Usage description</a></p>
* <p><a href="#usage_ex">Some usage examples</a></p>

<h2><a name="descr">Description</a></h2>

From a genomic fasta file and its associated GFF, the program first scans the genome to retrieve all sequences
delimited by stop codons. Only sequences of at least 60 nucleotides long are kept by default.

Those so-called ORF sequences are then annotated depending upon GFF element type(s) used as a reference.
The CDS element type is always used as a reference but others can be added.

By default an ORF sequence has 5 possible annotations:

| ORF annotation | Condition |
| --- | --- |
| c_CDS |if the ORF overlap with a CDS in the same phase |
| nc_5-CDS | if the 5' extremity of the c_CDS is at least 60 nucleotides long |
| nc_3-CDS | if the 3' extremity of the c_CDS is at least 60 nucleotides long |
| nc_ovp-CDS | if the ORF overlap with a CDS in a different phase |
| nc_intergenic | if the ORF do not overlap with anything |
 
**Note:** 
If an ORF sequence is tagged as 'c_CDS', this sequence is further processed to be cut at its 5' and 3' extremities that do not overlap with the CDS. If their length is above or equal to 60 nucleotides, then those subsequences can be assigned as nc_5-CDS and/or nc_3-CDS.
 <br></br>
 <br></br>
 
The user can also specify what GFF element type(s) can be used as reference(s) to annotate ORF sequences in addition to the CDS type. For instance, if the user adds the tRNA element type, ORF sequences could now be assigned as nc_ovp-tRNA if they overlap with a tRNA. Thus 6 assignments would now be possible for an ORF sequence:

| ORF annotation | Condition |
| --- | --- |
| c_CDS |if the ORF overlap with a CDS in the same phase |
| nc_5-CDS | if the 5' extremity of the c_CDS is at least 60 nucleotides long |
| nc_3-CDS | if the 3' extremity of the c_CDS is at least 60 nucleotides long |
| nc_ovp-CDS | if the ORF overlap with a CDS in a different phase |
| nc_ovp-tRNA | if the ORF overlap with a tRNA |
| nc_intergenic | if the ORF do not overlap with anything |

**Note on default parameters**:
* CDS is the only element type used as a reference to annotate ORF sequences.
* the minimum nucleotide number required to consider an ORF sequence is set at 60 nucleotides
* an ORF sequence is considered as overlapping with an element (e.g. CDS) if at least 70 % of its sequence overlap with the element or if this element is totally included within the ORF sequence


<h2><a name="installation">Installation</a></h2>

### 1. Download and uncompress the latest release archive

#### Download the latest release
Latest release: 
[ ![](./documentation/images/download-flat/16x16.png "Click to download the latest release")](https://github.com/nchenche/orfmap/releases/latest/)

#### Uncompress the archive
If you downloaded:
* the *.zip* file: ```unzip orfmap-x.x.x.zip```
* the *.tar.gz* file: ```tar xzvf orfmap-x.x.x.tar.gz```


### 2. Create an isolated environment
Although not strictly necessary, this step is highly recommended (it will allow you to work on different projects without having
any conflicting library versions).
 
#### Install virtualenv
``` python
python3 -m pip install virtualenv
```

#### Create a virtual python3 environment
```bash
virtualenv -p python3 my_env
```

#### Activate the created environment
```bash
source my_env/bin/activate
```

Once activated, any python library you'll install using pip will be installed solely in this isolated environment.
Every time you'll need to work with libraries installed in this environment (i.e. work on your project), you'll have
to activate it. 

Once you're done working on your project, simply type `deactivate` to exit the environment.


### 3. Install ORFMap in your isolated environment

Be sure you're virtual environment is activated, and then follow the procedure described below.

#### Go to the ORFMap directory
 
```bash
cd orfmap-x.x.x/
```

#### Install 

```python
python setup.py install
```

or 
```python
pip install .
```


<h2><a name="usage_descr">Usage description</a></h2>

To see all options available:

```
run_orfmap -h
```

This command will show:

<pre>usage: run_orfmap [-h] -fna [FNA] -gff [GFF] [-chr [CHR]] [-types_only TYPES_ONLY [TYPES_ONLY ...]]
                  [-types_except TYPES_EXCEPT [TYPES_EXCEPT ...]] [-o_include O_INCLUDE [O_INCLUDE ...]] [-o_exclude O_EXCLUDE [O_EXCLUDE ...]]
                  [-orf_len [ORF_LEN]] [-co_ovp [CO_OVP]] [-out [OUT]] [--show-types] [--show-chrs]

Genomic mapping of pseudo-ORF

optional arguments:
  -h, --help            show this help message and exit
  -fna [FNA]            Genomic fasta file (.fna)
  -gff [GFF]            GFF annotation file (.gff)
  -chr [CHR]            Chromosome name
  -types_only TYPES_ONLY [TYPES_ONLY ...]
                        Type feature(s) to use as reference(s) (&apos;CDS&apos; in included by default).
  -types_except TYPES_EXCEPT [TYPES_EXCEPT ...]
                        Type feature(s) to not consider as reference(s) (None by default).
  -o_include O_INCLUDE [O_INCLUDE ...]
                        Type feature(s) and/or Status attribute(s) desired to be written in the output (all by default).
  -o_exclude O_EXCLUDE [O_EXCLUDE ...]
                        Type feature(s) and/or Status attribute(s) desired to be excluded (None by default).
  -orf_len [ORF_LEN]    Minimum number of nucleotides required to define a sequence between two consecutive stop codons as an ORF sequence (60
                        nucleotides by default).
  -co_ovp [CO_OVP]      Cutoff defining the minimum CDS overlapping ORF fraction required to label on ORF as a CDS. By default, an ORF sequence
                        will be tagged as a CDS if at least 70 per cent of its sequence overlap with the CDS sequence.
  -out [OUT]            Output directory
  --show-types          Print all type features
  --show-chrs           Print all chromosome names
</pre>

Except -fna and -gff arguments that are mandatory, all others are optional.


### Basic run

ORFMap requires two input files: 
* a genomic fasta file (-fna)
* its associated GFF file (-gff).


The most basic run can be executed by typing:

```
run_orfmap -fna mygenome.fna -gff mygenome.gff
```

All of the ORF sequences are annotated relative to the CDS element type only. Thus 5 possible annotations are possible:

| ORF annotation | Condition |
| --- | --- |
| c_CDS |if the ORF overlap with a CDS in the same phase |
| nc_5-CDS | if the 5' extremity of the c_CDS is at least 60 nucleotides long |
| nc_3-CDS | if the 3' extremity of the c_CDS is at least 60 nucleotides long |
| nc_ovp-CDS | if the ORF overlap with a CDS in a different phase |
| nc_intergenic | if the ORF do not overlap with anything |


The output will be two separated files with the prefix "mapping_orf_":
* mapping_orf_mygenome.fa: 	a proteic fasta file of all the ORFs sequences found
* mapping_orf_mygenome.gff:	A GFF file describing all the ORFs sequences found
  
By default, the two output files will contain all possible 5 annotations mentionned above.


<h2><a name="usage_ex">Some usage examples</a></h2>

By default, all element types (except 'region' and 'chromosome') found in the GFF file are used as reference
to annotate ORF sequences. If an ORF sequence overlaps with more than 2 elements, then the ORF sequence will be assigned
according to the element with which it overlaps the most. For instance, let's say an ORF sequence overlaps at 85% with
a tRNA and at 90% with a sRNA, then the ORF will be assigned as nc-ovp_sRNA.
Note that the CDS element type always has the priority relative to any other element types. Therefore, if an ORF 
sequence overlaps at 72% with a CDS and at 95% with an other element that is not a CDS, then the ORF will be assigned as
c_CDS. When an ORF sequence entirely overlaps with multiple elements, then the choice for its  assignment is quite
arbitrary : the ORF will be assigned depending on the first element met in the GFF. That case could appear for 
intrinsically related elements such as gene, exon and mRNA. For example, let's say an ORF sequence equally overlaps with
an exon and a gene region (but there's no overlap with the CDS part). Since the gene normally appears firt in the GFF 
file, the ORF will be assigned as nc-ovp_gene. In order to avoid those special cases, an option allows the user specify 
element types that should not be considered as reference for the ORF assignment. 




In the case where an ORF sequence overlaps 

##### Use tRNA and snRNA element as a reference to annotate ORF sequences:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -types_only tRNA snRNA -out myResults
```



##### Use tRNA and snRNA element as a reference to annotate ORF sequences:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -types_only tRNA snRNA -out myResults
```

##### Write in output files only ORF sequences mapped as nc_ovp-tRNA and nc_ovp-snRNA:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -types_only tRNA snRNA -o_include nc_ovp-tRNA nc_ovp-snRNA -out myResults
```

##### Write in output files all ORF sequences except those mapped as c_CDS:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -type tRNA snRNA -o_exclude c_CDS -out myResults
```

##### or:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -type tRNA snRNA -o_exclude coding -out myResults
```

<em>Note</em>:
<p>
-o_include and -o_exclude take either feature types or a status attribute as arguments.
Feature types have to be amongst the possible annotations for ORF sequences (e.g. c_CDS, nc_5-CDS, nc_intergenic...)
 while status attribute is either 'coding' or 'non-coding' ('coding' refers to c_CDS and 'non-coding' refers to the other ones).
 </p>


##### Assign ORF seqences if stop-to-stop length is at least 50 nucleotides:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -orf_len 50
```

##### Consider an ORF sequence as overlapping with any element if at least 60 % of its sequence overlap with the element:
```
run_orfmap -fna mygenome.fna -gff mygenome.gff -co_ovp 0.6
```


