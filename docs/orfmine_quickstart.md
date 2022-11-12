
# Quick start

Here we first explain how to [prepare the folders](#prepare_folders) necessary for ORFmine and [launch the container](#start_container).


In addition, are presented several examples of the use of the different programs of ORFmine:


* [annotation and extraction](#ORF_annot) of all the ORFs (coding and noncoding) of a genome with ORFtrack and its companion program ORFget
* [prediction of the fold potential](#ORF_pred) of all noncoding ORFs of a genome with ORFold
* [characterization of the translation activity](#ORF_ribo) of all intergenic ORFs of a genome with ORFribo
* [estimation of the ages](#ORF_ages) of all the CDSs of a genome with ORFdate



All the inputs used in these examples are
available in the ORFmine/examples/ folder that can be dowloaded [here]().

<a name="prepare_folders"></a>
## Prepare your directories

First, you need to prepare your directories (that will depend on the usage) with the appropriate data. These directories will then be bound to your container when launching the container, so that your data will be accessible to the ORFmine programs (see [here](#start_container) for launching the container).

Your future container will have the following architecture:

<pre>
/
├── database
├── fastq
└── workdir
    ├── orfdate
    ├── orfold
    └── orfribo
</pre>


### database directory (any usage)
Several files must be provided in order to run the different ORFmine programs. They need to be
placed in the same directory whatever the programs to be used. This directory will then be bound to the /database directory of your future container.

* For ORFtrack:
    * A fasta file with the reference genome sequence (e.g. genome.fa)
    * A gff file with the reference genome annotation (e.g. genome.gff3)
* For ORFribo:
    * A fasta file with the reference genome sequence (e.g. genome.fa)
    * A gff file with the reference genome annotation (e.g. genome.gff3)
    * A fasta file with the sequences you do not want to treat and that you want to remove from the mapping, usually rRNA sequences (e.g. outRNA.fa)
* For ORFold:
    * A fasta file with the amino acid sequences to treat
    * A gff file (optional) if you want to map the HCA, Tango or IUPred scores along the chromosome and to see them with a genome browser like IGV.
* For ORFdate:
    * Fasta files of the complete proteomes (amino acid sequences) of the focal species and its neighbors
    * The corresponding tree (nwk format)
    * Table (2 columns) with the paths to the fasta files of proteomes and the names of the corresponding species

More details are given in the corresponding sections.

All these files whatever the programs you will run must be in the same local directory (though subdirectories can be created for more clarity but they must be in the same global directory) that will constitute the /database directory of the future container. You are free to choose the name and the path of the folder in your local machine.

### working directory (any usage)
When creating the container, please be aware that you will choose a directory that will be your working directory and where some of the outputs will be written.

* Particular case of ORFribo:
    * The working directory must contain the configuration file called *config.yaml*, completed by the user (see [here](#ORF_ribo) for more details on how to complete the *config.yaml*).

This directory will then be bound on the /workdir directory of your future
container. You are free to choose the name and the path of the folder in your
local machine.

### fastq directory (ORFribo)
If you want to use ORFribo, you must place the Ribosome Profiling data (i.e. the fastq(.gz) files) in a specific folder **without anything else inside it**.
You can have this folder wherever you want on your computer and name
it as you like. This directory will then be bound to the /fastq directory of your future container.


<a name="start_container"></a>

## Launch the container

To launch the container and get an interactive shell inside it, you can use either
the "singularity " or the "docker" commands. Every folder that you have prepared [before](#prepare_folders) must be bound to their corresponding folders inside the container through the -v or -B options (see the command lines below). In the following command lines, every italics part has to be adapted from your computer's organization (names and paths) and the other parts are part of the container structure and **must remain as they are!**

**Using singularity**
<pre>
<code class="language-bash">singularity shell --pwd /workdir/ \
-B <i>PATH_TO_YOUR_WORKING_DIR</i>:/workdir/ \
-B <i>PATH_TO_YOUR_DATA_DIR</i>:/database/ \
-B <i>PATH_TO_YOUR_FASTQ_DIR</i>:/fastq/ \
<i>PATH_TO_YOUR_ORFMINE_IMAGE</i>.sif</code>
</pre>


**Using docker**

The same logic applies with the docker command except that **absolute paths** must be given:

<pre>
<code class="language-bash">docker run --rm -it \
-v <i>ABSOLUTE_PATH_TO_YOUR_WORKING_DIR</i>:/workdir/ \
-v <i>ABSOLUTE_PATH_TO_YOUR_DATA_DIR</i>:/database/ \
-v <i>ABSOLUTE_PATH_TO_YOUR_FASTQ_DIR</i>:/fastq/ \
annelopes94/orfmine:v0.8.7</code>
</pre>


Please note that if you do not want to use ORFribo, you do not need to link the fastq folder to the container and can remove the binding part of the command line: <code><i>PATH_TO_YOUR_FASTQ_DIR</i>:/fastq/</code>. See [here](./orfmine_binding.md) for more details on how to bind your local folders to those of the container.


#### Binding Tango and/or IUPred

If you want to use Tango and/or IUPred with ORFfold, you will need to have those softwares installed on your local computer (choose the Unix compatible source code for Tango). 

Once installed, please place the parent directory content of those softwares inside a directory of your choice. Below is a description of the workflow we followed:

```bash
# let's consider Tango is installed in ~\tango and IUPred in ~/iupred2a

# create a directory that will contain Tango and IUPred source codes content
mkdir ~/softwares

# place Tango and IUPred source codes content in softwares
cp ~/tango ~/softwares
cp -r ~/iupred2a/* ~/softwares
```

To use Tango and/or IUPred on the container, you must launch your container by binding the directory created above to a specific directory on the container

For singularity:
<pre>
<code class="language-bash">singularity shell --pwd /workdir/ \
-B <i>PATH_TO_YOUR_WORKING_DIR</i>:/workdir/ \
-B <i>PATH_TO_YOUR_DATA_DIR</i>:/database/ \
-B <i>PATH_TO_YOUR_FASTQ_DIR</i>:/fastq/ \
-B <i>PATH_TO_YOUR_DIR_CONTAINING_TANGO_OR_IUPRED</i>:/ORFmine/orfold_v1/orfold/softwares/ \
<i>PATH_TO_YOUR_ORFMINE_IMAGE</i>.sif</code>
</pre>
For docker:
<pre>
<code class="language-bash">docker run --rm -it \
-v <i>ABSOLUTE_PATH_TO_YOUR_WORKING_DIR</i>:/workdir/ \
-v <i>ABSOLUTE_PATH_TO_YOUR_DATA_DIR</i>:/database/ \
-v <i>ABSOLUTE_PATH_TO_YOUR_FASTQ_DIR</i>:/fastq/ \
-v <i>ABSOLUTE_PATH_TO_YOUR_DIR_CONTAINING_TANGO_OR_IUPRED</i>:/ORFmine/orfold_v1/orfold/softwares/ \
annelopes94/orfmine:v0.8.7</code>
</pre>


<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">

    <ul>
    <li> When launching the container with docker, the <i>annelopes94/orfmine:v0.8.7</i> part can be changed if you want to launch a specific version of ORFmine. For example : <i>annelopes94/orfmine:v0.8.6</i>.
     </li>

        <li>  If you do not know which versions you have on your computer,
        you can print them by typing <code>docker image ls</code>
         </li>

         <li>
         WARNING: please keep in mind that, once launched, you work <i>inside</i> the container. Filenames in the
command lines must, therefore, be indicated either alone (i.e. without any path) or prefixed with the path of the folder as it is
inside the container (i.e. path relative to the container not to your local machine - e.g. /database, /fastq, /workdir).
For example, when running ORFold, the fasta filenames must be indicated in the command line as follows:
<code>orfold -faa sequences.fasta</code>
OR
<code>orfold -faa /database/sequences.fasta</code>
         Again, do not indicate the path of the sequences.fasta file according to your local machine architecture - the container only refers to its relative paths.

         </li>

          <li>
          If you want autocompletion when writing the filenames, you need to
        precise the path of the folder containing the corresponding files (relative path according to the container).
         </li>

    <li>
    You can verify the content of your directories as follows:
<code>ls /workdir/</code>
This should display all the files and folders of the working directory.
Or :
<code>ls /database/</code>
This should display all the files and folders of the data folder containing all the data files (fasta, gff, ...).
    </li>
     </ul>
    </p>
</div>


<br>


# Examples of ORFmine usage

To guide the user, all command line examples will rely on the examples directory that can be found [here](). If you want to follow the guidelines, please follow the procedure below:
```bash
# place examples.zip on a directory and unzip it
mkdir ~/orfmine_tutorial
mv ~/Downloads/examples.zip ~/orfmine_tutorial
cd ~/orfmine_tutorial/
unzip examples.zip
cd examples/

# launch the container (docker usage case).
docker run --rm -it -v $(pwd)/workdir:/workdir/ -v $(pwd)/database:/database/ -v $(pwd)/fastq:/fastq/ annelopes94/orfmine:v0.8.7

# or with singularity (considering your .sif image is there: ~/orfmine/orfmine_v0.8.7.sif)
singularity shell --pwd /workdir/ -B workdir:/workdir/ -B database</i>:/database/ -B fastq</i>:/fastq/ ~/orfmine/orfmine_v0.8.7.sif
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

* /fastq/XYZ.fastq


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


ORFribo generates many intermediate files that can be useful for further analysis. Their description is available [here](./orfribo_outputs.md). The main output is the table named genome_mean70_median70_reads_concatenated.tab that summarizes the results for each ORF of interest (e.g. nb and fractions of F0, F1 and F2 reads). An example of this output table can be found in the /workdir/orfribo/RESULTS/Bam2Reads_genome_output/SRR1520313_17031088/ directory. A full description of the summary table can be found [here](./orfribo_outputs.md#main_table).



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
