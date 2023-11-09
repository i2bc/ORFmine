# Running example



In this section, we present several examples of the use of the different programs of ORFmine:

* [get started](#get-started): downloading the examples and launching the container
* [annotation and extraction](#ORF_annot) of all the ORFs (coding and noncoding) of a genome with ORFtrack and its companion program ORFget
* [prediction of the fold potential](#ORF_pred) of all noncoding ORFs of a genome with ORFold
* [characterization of the translation activity](#ORF_ribo) of all intergenic ORFs of a genome with ORFribo
* [estimation of the ages](#ORF_ages) of all the CDSs of a genome with ORFdate



## Get started!

To guide the user, all command line examples will rely on the examples directory that can be found [here](http://bim.i2bc.paris-saclay.fr/anne-lopes/ORFmine_examples). If you want to follow the guidelines, please follow the procedure below:

```bash
# create a directory that will host the examples archive and go in it
mkdir ~/orfmine_tutorial
cd orfmine_tutorial

# download the example archive if not already done
wget http://bim.i2bc.paris-saclay.fr/anne-lopes/ORFmine_examples/ORFmine_examples.zip

# otherwise, if you have already downloaded the example archive, place it in the orfmine_tutorial/ directory

# untar the archive
unzip ORFmine_examples.zip

# go into the example directory 
cd ORFmine/examples


# launch the container (docker usage case).
docker run --rm -it -v $(pwd)/workdir:/workdir/ -v $(pwd)/database:/database/ -v $(pwd)/fastq:/fastq/ lopesi2bc/orfmine:latest /bin/bash

# or with singularity (considering your .sif image is there: ~/orfmine/orfmine_latest.sif)
singularity shell --pwd /workdir/ -B workdir:/workdir/ -B database:/database/ -B fastq:/fastq/ ~/orfmine/orfmine_latest.sif
```


After following this procedure, you could juste copy/paste all of the commands below to test them directly.

<a name="ORF_annot"></a>

## Annotation and extraction of ORFs with ORFtrack

### Annotation of ORFs

ORFtrack takes as input (i) a fasta file containing the nucleotide sequences of all chromosomes or contigs, and (ii) the corresponding annotations in a gff file (see [the GFF3 documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) for more details on the GFF3 format). Please notice that we strongly recommend using GFF3 formats.
The following command runs ORFtrack on the complete genome of *Saccharomyces cerevisiae* (*Scer*) available in the ORFmine/examples folder:

``` bash
orftrack -fna /database/Scer.fna -gff /database/Scer.gff
```
ORFtrack generates a novel gff file mapping_orf_Scer.gff that contains the annotations of all the identified ORFs (see [here](./orftrack_annotation.md) for more details on the annotation process). The output file is located in /database/. Indeed, all the ORFtrack and ORFget outputs are stored in /database/ since they might constitute the inputs of the other ORFmine programs. A copy of the output can be found in the
ORFmine/examples/database/ folder). Please note that ORFtrack also generates a summary.log file that lists all the ORF categories it has identified.


### Extraction of the sequences of noncoding ORFs in a fasta format
The amino acid sequences of all annotated ORFs or specific subsets of ORFs (i.e. only noncoding intergenic ORFs for example) can be extracted and written in a fasta file with ORFget, a tool provided with ORFtrack. The following instruction writes the amino acid sequences of all yeast noncoding ORFs (including noncoding intergenic ORFs and those that overlap with a genomic feature). ORFget therefore needs the genome sequence (nucleotides) from which the ORF sequences will be extracted and translated into amino acids and the ORFtrack output gff file with the annotations and coordinates of all the identified ORFs. The -features_include option indicates the features to be matched in the names of the ORFs that we want to extract. In this case, we are interested in all noncoding ORFs (i.e. all ORFs containing the flag "nc" - see [here](./orftrack_annotation.md) for more details on the annotation rules of the ORFs), so we indicate it to ORFget through the "-features_include nc" option (see [here](./orfget_run.md) for more examples).


``` bash
orfget -fna /database/Scer.fna -gff /database/mapping_orf_Scer.gff -features_include nc
```

ORFget produces a fasta file (mapping_orf_Scer_nc.pfasta) containing the amino acid sequences
of all noncoding ORFs annotated by ORFtrack (see [here](./orfget_run.md) for more examples
on the use of ORFget).The output fasta is stored in /database/ directory of the container and renamed based on the gff file rootname and the include features. The user can also extract the nucleotide sequences through the "-type" option (-type nucl) or both amino acid and nucleotide sequences (-type both).

### Extraction of the sequences of the CDSs of *yeast* in a fasta format
In addition, ORFget enables the reconstruction of all coding
sequences (CDSs) of the input genome (in this examlpe, Scer) according to their annotation in the **original**
gff file. When different isoforms are described, it will generate all of them according to the gff file. It, therefore, needs as input the original gff file (not the one created by ORFtrack which only contains STOP-to-STOP ORFs) and the genome sequence in FASTA format.


``` bash
orfget -fna /database/Scer.fna -gff /database/Scer.gff -features_include CDS
```

The output is generated in the /database/ directory of the container. An example of this output is available in the ORFmine/examples/database/ directory.

<br><br>

<a name="ORF_pred"></a>
## Characterization of the fold potential and disorder and aggregation propensities of the noncoding ORFs with ORFold


### Running ORFold
ORFold only needs a fasta file containing the amino acid
sequences to characterize. In the following example, ORFold will estimate the
foldability and the disorder and aggregation propensities of all the amino acid sequences potentially encoded in the noncoding ORFs of yeast.

``` bash
orfold -faa /database/mapping_orf_Scer_nc.pfasta  -options H
```

with -options H to indicate that we want to estimate the foldability (H for HCA). 


ORFold generates one table that is stored in the /workdir/orfold directory of the container
(a copy of the output can be found in ORFmine/examples/workdir/orfold directory). The table contains for each input
sequence, the fold potential
as calculated by HCA [1]
(see [here](./How_it_works_orfold.md) for more details on the calculation of this score).


In addition, a gff
file containing the annotation of the input sequences can be provided
as well. In this case, ORFold writes the fold potential
and disorder and aggregation propensities in the 9th column of three gff files respectively.
The latter can be subsequently uploaded on a genome viewer like IGV for a
visual inspection of the distribution of these properties along
the genome.

The following instruction enables the user to probe the fold
potential and disorder and aggregation propensities of all the yeast noncoding
ORFs annotated by ORFtrack and takes as inputs the corresponding fasta file and the gff file generated by ORFtrack (mapping_orf_Scer.gff located in /database/).

``` bash
orfold -faa /database/mapping_orf_Scer_nc.pfasta -options H -gff /database/mapping_orf_Scer.gff
```


In addition to the score table, ORFold has generated a new gff file containing for each ORF its fold potential. The output gff file is stored in the /workdir/orfold directory of the container and can be subsequently uploaded on a genome viewer.


### Visualization of the distribution of the fold potential with ORFplot

ORFplot enables the visualization of the distribution of the fold
potential of the amino acid sequences potentially encoded in noncoding ORFs along with the one of a reference dataset
of globular proteins taken from Mészáros et al. [2]. It takes as input the table generatd by ORFold that must be located in the /workdir/orfold/ directory of the container.

```bash
orfplot -tab /workdir/orfold/mapping_orf_Scer_nc.tab -names  "Yeast noncoding ORFs"
```

The output graphic is generated in PNG and PDF formats and stored in the /workdir/orfold directory of the container (see
[here](./Plot_orfold.md) for other examples of use of ORFplot). Each
distribution is compared with the one of the globular protein dataset
with a Kolmogorov Smirnov test. Asterisks on the plot denote level
of significance: * < 0.05, \**\** < 0.01, *** < 0.001.




<br><br>



<a name="ORF_ribo"></a>
## Characterization of the translation activity of all intergenic ORFs of yeast


ORFribo, based on ribosome-profiling data, probes the translation activity of all the ORFs of interest of a genome (be they coding or noncoding) that are annotated in a gff file. We strongly recommend using ORFtrack to [generate the gff file](#annotation-and-extraction-of-orfs-with-orftrack) since ORFtrack and ORFribo were jointly developed so that the output of ORFtrack fits the prerequisites of ORFribo. In this example, we investigate the translation of all intergenic ORFs of yeast annotated as "nc_intergenic" by ORFtrack in its output mapping_orf_Scer.gff. An example of this file that has been [generated by ORFtrack](#annotation-and-extraction-of-orfs-with-orftrack) and that will be used as input by ORFribo is available in the ORFmine/examples/database/ folder.


Here is presented a simple example that enables one to run ORFribo from inputs and a pre-configurated config.yaml file which are all available in the ORFmine/examples/ folder. However, before running ORFribo on this example dataset, we recommend the user read in details the [phylosophy](./orfribo_objectives.md) and [main steps](./How_it_works_orfribo.md) of ORFribo.

### Input files
ORFribo relies on the fastq file(s) (Ribo-Seq dataset(s)) that must be stored in the /fastq/ folder of the container - an example of such file is available in the ORFmine/examples/fastq folder

* /fastq/SRR1520313_17031088.fastq.gz


and 4 additional files that must be stored in the /database/ folder of the container. Examples of such files are available in the ORFmine/examples/database/ folder.

* Scer.fna: fasta file containing the nucleotide sequences
of the chromosomes of the yeast genome
* Scer.gff: the original annotation of yeast in a gff file
* mapping_orf_Scer.gff: the gff file containing all the ORFs detected by ORFtrack and including the nc_intergenic ORFs for which we want to estimate the translation activity. It therefore involves that you have previously generated this file [as presented here](#annotation-and-extraction-of-orfs-with-orftrack)
* Scer_rRNA.fa: A fasta file with the sequences that we want to remove from the mapping step, here rRNA sequences


### Preparing the config.yaml file
ORFribo is very easy to handle and only needs a configuration file to be edited before running it. The latter named *config.yaml* contains the parameters that can be adjusted by the user. The full description of this file is available [here](./orfribo_configuration.md). For this example, a pre-filled configuration file is present is the ORFmine/examples/workdir/ directory.


Please place the config.yaml file in the */workdir/* directory of your container and edit it if needed following the instructions presented [here](./orfribo_configuration.md):


### Running ORFribo

The only parameters to adjust in the command line are:

* The number of virtual CPUs/threads to reserve for the analysis

* The approximate amount of memory to reserve for the analysis (in GB)

ORFribo is launched like this:
``` bash
orfribo CPU MEMORY
```
For example, for an analysis launched with the use of 6 CPUs and around 10GB of memory, it would be:
``` bash
orfribo 6 10
```

Be careful to have enough memory on your computer/cluster for you analysis as a lack of RAM would lead the pipeline to crash without giving any error message. So if the program stops without any obvious reason this might be the reason why, and one should restart the pipeline with more resources.



### Main outputs


ORFribo generates many intermediate files that can be useful for further analysis. Their description is available [here](./orfribo_outputs.md). The main output is the table named all_samples_genome.25-35.mean70_median70_reads_concatenated.tab that summarizes the results for each ORF of interest (e.g. nb and fractions of F0, F1 and F2 reads). An example of this output table can be found in the examples.zip archive in the ORFmine/examples/workdir/orfribo/RESULTS/Bam2Reads_genome_output/ directory. A full description of the summary table can be found [here](./orfribo_outputs.md#main-outable-table).



<a name="ORF_ages"></a>
## Running ORFdate
ORFdate estimates the evolutionary ages of a set of ORFs (coding but also noncoding) based on phylostratigraphy. It takes as inputs:


* a distance phylogenetic tree in newick format, containing the species of interest (focal species) and its neighbors (see ORFmine/examples/database/orfdate_inputs/ for an example). All branches must be associated with their relative distance on which the age estimation will be based.
* a two columns csv file (i.e. columns separated by commas) with the path to the sequence fasta files of each species as first column, and the associated species names as indicated in the phylogenetic tree as second column (see ORFmine/examples/database/orfdate_inputs/ for an example). It, therefore, involves the fasta files of the complete proteomes (amino acid sequences) of the focal and its neighbors to be stored in the corresponding directories.
* the fasta files containing all the amino acid sequences of the focal and its neighboring species.

All the inputs must be stored in the /database/ directory, though a subdirectory can be created for more clarity (e.g. /database/orfdate_inputs/). The following command line will estimate the ages of 42 yeast CDSs (subsample of yeast CDSs that can be found in the file ORFmine/examples/database/orfdate_inputs/sample_42_Scer_CDS.pfasta) based on the phylogenetic tree of the Saccharomyces sensu stricto (i.e. <i>S. cerevisiae, S. paradoxus,  S. mikatae, S. kudrievezii, S. arboricola and S. bayanus </i>)(the corresponding tree can be found in ORFmine/examples/database/orfdate_inputs/Saccharomycew_species.nwk). Please note that the last node of the tree will be used as the upper limit for the age estimation. Consequently, all sequences with a match in the farest species, with respect to the focal, are associated with the same upper bounded estimated age regardless of whether they may have other matches in more distant species outside the input tree (see Papadopoulos et al, [3] for more details).



```
orfdate -focal Scer -tree /database/orfdate_inputs/Saccharomyces_species.nwk -names /database/orfdate_inputs/names.csv
```

with the -focal option corresponding to the name of the focal species as indicated in the tree. Please note that additional arguments including BLAST parameters can be provided. See [here](./orfdate_parameters.md) for more details.

An example of the names.csv file is available in the ORFmine/examples/database/orfdate_inputs/ directory as well as the 5 fasta files containing all the sequences of the 5 neighboring species of the focal.

ORFdate generates two outputs in csv format located in the /workdir/orfdate/ directory. Examples of these outputs can be found in the ORFmine/exmaples/workdir/orfdate/ folder.

* sample_42_Scer_CDS_hits.csv which reports the number of significant matches for all sequences of the focal species (rows) across the neighbors species of the tree (columns)
* sample_42_Scer_CDS_dated.csv which reports for each sequence (col 1), the farest neighbor species for which a significant math has been found (col2), and the time of divergence between this neighbor and the focal species in Mya (col3)

ORFdate generates also intermediate files (blast databases and blast outputs) which are stored in the /workdir/orfdate/blastdb/ and /workdir/orfdate/blastout/ directories. Examples of such outputs can be found in the ORFmine/examples/workdir/orfdate/ directory.

<br><br><br>





#### References


1. Bitard-Feildel, T. & Callebaut, I. HCAtk and pyHCA: A Toolkit and Python API for the Hydrophobic Cluster Analysis of Protein Sequences. bioRxiv 249995 (2018).
2. 	Mészáros B, Erdős G, Dosztányi Z (2018) IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46:W329–W3372. 
3. Papadopoulos, C., Arbes, H., Chevrollier, N., Blanchet, S., Cornu, D., Roginski, P., Rabier, C., Atia, S., Lespinet, O., Namy, O., Lopes, A. (submitted).
