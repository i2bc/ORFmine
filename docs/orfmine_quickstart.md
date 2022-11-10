Here we first explain how to [prepare the folders](#prepare_folders) and [launch the container](#start_container).


In addition, are presented several examples of the use of the different programs of ORFmine:


* [annotation and extraction](#ORF_annot) of all the ORFs (coding and noncoding) of a genome with ORFtrack and its companion program ORFget
* [prediction of the fold potential and disorder and aggregation propensities](#ORF_pred) of all noncoding ORFs of a genome with ORFold
* [characterization of the translation activity](#ORF_ribo) of all intergenic ORFs of a genome with ORFribo
* [estimation of the ages](#ORF_ages) of all the CDSs of a genome with ORFdate



All the inputs used in these examples are
available in the ORFmine/examples/ folder.

<a name="prepare_folders"></a>
## Prepare your folders

First, you need to prepare your directories (that will depend on the usage) with the appropriate data. These directories will then be bound to your container when launching the container, so that the associated data will be accessible to the ORFmine programs (see [here](#start_container) for launching the container).

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


### database folder (any usage)
Several files must be provided in order to run the different ORFmine programs. They need to be placed in the same directory whatever the programs to be used. This directory will then be bound to the /database directory of your future container.

* For ORFtrack:
    * A fasta file with the reference genome sequence (Ex : genome.fa)
    * A gff file with the reference genome annotation (Ex : genome.gff3)
* For ORFribo:
    * A fasta file with the reference genome sequence (Ex : genome.fa)
    * A gff file with the reference genome annotation (Ex : genome.gff3)
    * A fasta file with the sequences you do not want to treat and that you want to remove from the mapping, usually rRNA sequences (Ex : outRNA.fa)
* For ORFold:
    * A fasta file with the amino acid sequences to treat
    * A gff file (optional) if you want to map the HCA, Tango or IUPred scores along the chromosome and to see them with a genome browser like IGV.
* For ORFdate:
    * Fasta files of the complete proteomes (amino acid sequences) of the focal and its neighbors
    * The corresponding tree (nwk format)
    * Table (2 columns) with the path to the proteome fasta files and the names of the species

More details are given in the corresponding sections

All these files whatever the programs you will run must be in the same local directory (though subdirectories can be created for more clarity but they must be in the same global directory) that will constitute the /database directory of the future container. You are free to choose the name and the path of the folder in your local machine.

### working directory (any usage)
When creating the container, please be aware that you will choose a directory that will be your working directory and where some of the outputs will be written.

* Particular case of ORFribo:
    * The working directory must contain the configuration file called *config.yaml*, completed by the user (see [here](#ORF_ribo) for more details on how to complete the config.yaml).

This directory will then be bound on the /workdir directory of your future container. You are free to choose the name and the path of the folder in your local machine.

### fastq folder (ORFribo)
If you want to use ORFribo, you must place the Ribosome Profiling data (i.e. the fastq(.gz) files) in a specific folder **without anything else inside it**.
You can have this folder wherever you want on your computer and name it as you like. This directory will then be bound to the /fastq directory of your future container.


<a name="start_container"></a>
## Launch the container

To launch the container and get an interactive shell inside it, you can use either the "singularity " or the "docker" commands. Every folder that you have prepared [before](#prepare_folders) must be bound to their corresponding folders inside the container through the -v or -B options (see the command line below). In the following command lines, every italics part has to be adapted from your computer's organization (names and paths) and the other parts are part of the container structure and **must remain as they are!**

**Using singularity**

<code>
singularity shell --pwd /workdir/ -B */path/to/your/local/working/folder/*:/workdir/ -B */path/to/your/local/data/folder/*:/database/ -B */path/to/your/local/fastq/folder/*:/fastq/ */path/to/ORFmine/image*.sif
</code>


**Using docker**

The same logic applies with the docker command but the paths **must be absolute paths**:

<code>
docker run --rm -it -v */absolute/path/to/your/local/working/folder/*:/workdir/ -v */absolute/path/to/your/local/data/folder/*:/database/ -v */absolute/path/to/your/local/fastq/folder/*:/fastq/ *orfmine:latest* /bin/bash
</code>

Please note that if you do not want to use ORFribo, you do not need to link the fastq folder to the container and can remove the -v or -B /absolute/path/to/your/local/fastq/folder/:/fastq/ part of the command line. See [here](./orfmine_binding.md) for more details on how to bind your local folders to those of the container.

You are now ready to use all the ORFmine programs and start your analyses!




<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">

    <ul>
    <li> When launching the container, the <i>orfmine:latest</i> part can be changed if you want to launch a specific version of ORFmine. For example : <i>orfmine:v0.7.11</i>.
     </li>

        <li>  If you do not know which versions you have on your computer,
        you can print them by typing :
        `docker image ls`
         </li>

         <li>
         WARNING: please keep in mind that, once launched, you work **inside** the container. Filenames in the command lines must, therefore, be indicated either alone (i.e. without any path) or prefixed with the path of the folder as it is inside the container (i.e. path relative to the container not to your local machine - e.g. /database, /fastq, /workdir). For example, when running ORFold, the fasta filenames must be indicated in the command line as follows:
<code>orfold -fna sequences.fasta</code>
OR
<code>orfold -fna /database/sequences.fasta</code>
         </li>

          <li>
          Please do not indicate the path of the sequences.fasta file according to your local machine architecture - the container only refers to its relative paths.

         </li>

          <li>
          If you want autocompletion when writing the filenames, you need to precise the path of the folder containing the corresponding files (relative path according to the container).
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




<br><br>

<a name="ORF_annot"></a>
## Annotation and extraction of ORFs with ORFtrack

### Annotation of ORFs

Before starting, be sure your container [is activated](#start_container) and your sequence inputs (fasta and gff files) are stored in the container's /database directory (see [here](#prepare_folders) for the organization of data). All the following commands rely on paths relative to the docker architecture, and therefore, assume the gff and fasta inputs are stored in */database/*.


ORFtrack takes as input (i) a fasta file containing the nucleotide sequences of all chromosomes or contigs, and (ii) the corresponding annotations in a gff file (see [the GFF3 documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) for more details on the GFF3 format). Please notice that we strongly recommend using GFF3 formats.
The following command runs ORFtrack on the complete genome of Saccharomyces cerevisiae (*Scer*) available in the "examples" folder:

``` bash
orftrack -fna /database/Scer.fna -gff /database/Scer.gff
```
ORFtrack generates a novel gff file mapping_orf_Scer.gff that contains the annotations of all the identified ORFs (see [here](./orftrack_annotation.md) for more details on the annotation process). The output file is located in /database. Indeed, all the ORFtrack and ORFget outputs are stored in /database since they might constitute the inputs of the other ORFmine programs. A copy of the output can be found in the ORFmine/examples/database folder). Please note that ORFtrack also generates a summary.log file that lists all the ORF categories it identified.


### Extraction of noncoding ORFs in a fasta file
The amino acid sequences of all annotated ORFs or specific subsets of ORFs (i.e. only noncoding intergenic ORFs for example) can be extracted and written in a fasta file with ORFget, a tool provided with ORFtrack. The following instruction writes the amino acid sequences of all yeast noncoding ORFs (including noncoding intergenic ORFs and those that overlap with a genomic feature). ORFget therefore needs the genome sequence (nucleotides) from which the ORF sequences will be extracted and translated into amino acids and the ORFtrack output gff file with the annotations and coordinates of all the identified ORFs. The -features_include option indicates the features to be matched in the names of the ORFs that we want to extract. In this case, we are interested in all noncoding ORFs (i.e. all ORFs containing the flag "nc" - see [here](./orftrack_annotation.md) for more details on the annotation rules of the ORFs), so we indicate it to ORFget through the "-features_include nc" option (see [here](./orfget_run.md) for more examples).


``` bash
orfget -fna /database/Scer.fna -gff /database/mapping_orf_Scer.gff -features_include nc
```

ORFget produces a fasta file (Scer_noncoding.pfasta) containing the amino acid sequences of all noncoding ORFs annotated by ORFtrack (see [here](./orfget_run.md) for more examples on the use of ORFget).The output fasta is stored in /database directory of the container. The user can also extract the nucleotide sequences through the "-type" option (-type nucl) or both amino acid and nucleotide sequences (-type both).

### Extraction of protein sequences of *yeast* in a fasta file
In addition, ORFget enables the reconstruction of all coding sequences (CDSs) of the input genome (in this examlpe, Scer) according to their annotation in the **original** gff file. When different isoforms are described, it will generate all of it according to the gff file. It, therefore, needs as input the original gff file (not the one created by ORFtrack which only contains STOP-to-STOP ORFs) and the genome sequence in FASTA format.


``` bash
orfget -fna /database/Scer.fna -gff /database/Scer.gff -features_include CDS
```

The output is stored in the /database directory of the container.

<br><br>

<a name="ORF_pred"></a>
## Characterization of the fold potential and disorder and aggregation propensities of the noncoding ORFs with ORFold


### Running ORFold
ORFold only needs a fasta file containing the amino acid sequences to characterize. In the following example, ORFold will estimate the foldability and the disorder and aggregation propensities of all the noncoding ORFs of yeast.

``` bash
orfold -fna /database/Scer_noncoding.pfasta  -options HIT
```
with -options HIT to indicate that we want to estimate the foldability (H for HCA), the disorder (I for IUPred) and the aggregation propensity (T for Tango). "-options H" will only estimate the foldablity.

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        The option -options HIT implies that IUPred and Tango have been downloaded and that the full installation (including IUPred and Tango) has been undertaken.
    </p>
</div>


ORFold generates one table that is available in the /workdir/orfold directory of the container (a copy of the output can be found in ORFmine/examples/workdir/orfold directory). The table contains for each input sequence, the fold potential, disorder and aggregation propensities based on HCA [3], IUPred [4][5][6] and TANGO [7][8][9], respectively (see [here](./How_it_works_orfold.md) for more details on the calculation of these scores).


In addition, a gff file containing the annotation of the input sequences can be provided as well. In the latter case, ORFold writes the fold potential and disorder and aggregation propensities in three gff files respectively. The latter can be subsequently uploaded on a genome viewer like IGV for a visual inspection of the distribution of these properties along the genome.

The following instruction enables the user to probe the fold potential and disorder and aggregation propensities of all the yeast noncoding ORFs annotated by ORFtrack and takes as inputs the corresponding fasta file and the gff file generated by ORFtrack (mapping_orf_Scer.gff located in /database).

``` bash
orfold -fna /database/Scer_noncoding.pfasta -options HIT -gff /database/mapping_orf_Scer.gff
```


In addition to the score table, ORFold also generates three new gff files containing for each ORF its fold potential and disorder and aggregation propensities respectively. The output gff files are stored in the /workdir/orfold directory of the container and can be subsequently uploaded on a genome viewer.


### Visualization of the distribution of the fold potential with ORFplot

ORFplot enables the visualization of the distribution of the fold potential of the noncoding ORFs and protein sequences along with the one of a reference dataset of globular proteins taken from Mészáros et al. [8].

``` bash
orfplot -tab /workdir/orfold/Scer_noncoding.tab -names  "Yeast noncoding ORFs"
```

The output graphic is generated in PNG and PDF formats and stored in the /workdir/orfold directory of the container (see [here](./Plot_orfold.md) for other examples of use of ORFplot). Each distribution is compared with the one of the globular protein dataset with a Kolmogorov Smirnov test. Asterisks on the plot denote level of significance: * < 0.05, \**\** < 0.01, *** < 0.001.




<br><br>



<a name="ORF_ribo"></a>
## Characterization of the translation activity of all XXX of a genome


### ORFribo purpose and main steps
ORFribo is designed to analyse, based on ribosome-profiling data, the translation activity of all the ORFs of a genome (coding and noncoding) that are annotated in a gff file. We strongly recommend using ORFtrack to [generate the gff file](#annotation-and-extraction-of-orfs-with-orftrack) since ORFtrack and ORFribo were jointly developed so that the output of ORFtrack fits the prerequisite of ORFribo. However, the user might provide its own gff. In this case, the name of the ORF category to be analyze (e.g. intergenic, xyz...) must be given in the 3rd column of the gff file and in the "final_counts" attribute of the config.yaml file (more details in the [dedicated section](#ORF_categories) - an example of the gff output generated by ORFtrack can be found in the ORFmine/examples/database directory). We recall that ORFribo does not predict the translation status of ORFs. In contrast, it calculates the number of reads that map on the genomic coordinates of each ORF of the ORF category of interest and detects their corresponding frame. As a result, for each ORF, it provides the number of reads that map in-frame (i.e. in the reading frame of the ORF or in its F0 frame) and the number of reads that map in its +1 and +2 frames. The user can then classify and/or infer the translation status of each ORF according to these numbers. In other words, ORFribo has no a priori and aims at calculating these numbers for you, you are then free to take decision based on these numbers and your own criteria.


The main steps of ORFribo are:

* Determination of the P-site offset: For each read size or kmer, ORFribo tries to find the distance in nucleotides from the 5'-end of the read (begining of the Ribosome Protected Fragment - RPF) to the first base of the ribosome's P-site, thanks to riboWaltz [12]. The detection of the first base of the ribosome's P-site enables the identification of the translated codon and, thus, of the translated frame. This is necessary to get the quality controls about phasing in coding regions.


* Alignment on coding sequences (CDSs) and detection of kmers of good quality: Read phasing is a very important step to make sure we are able to detect the frame that is translated among the three frames of the RNA. In coding regions, the frame that is expected to be translated (i.e. the coding one) is already known and will be used to estimate the quality of the experiment and more precisely, to identify read kmers of good quality. To do so, ORFribo will retain for the next alignment step, only kmers whose P-site position as predicted by riboWaltz indicate, for more than 70% of the corresponding kmers, codons that belong to the phase of the CDSs (P0, i.e. reading frame of the CDS). 70% is the default value but it can be modified by the user through the config.yaml file.


* Alignment on all ORFs of interest: Read kmers which have passed the phasing filter on CDSs are aligned on the set of [ORFs indicated by the user](#ORF_categories) (e.g. intergenic ORFs, alternative ORFs...). The ORF categories to be analyzed by ORFribo must correspond to those identified by ORFtrack and provided in the 3rd column of its gff output (see an example ORFmine/examples/database/mapping_orf_Scer.gff) (see [here](./orftrack_annotation.md) for more details on the ORF categories and annotation process). If another gff file is provided, please have a look on the format of the ORFtrack gff output. The P-site offset detected for each kmer in the previous step allows to determine the proportion of reads in each phase of every tested ORF, thereby estimating its fractions of in-frame reads (F0 or P0), and in its +1 and +2 frames. This can help the user identify ORFs with high levels of translation specificity (high fractions of F0 reads) or classify ORFs according to their levels of translation specificity (i.e. different levels of F0 reads).

### Before launching ORFribo: prepare the config.yaml
ORFribo is very easy to handle and only needs a configuration file to be edited before running it. The latter named *config.yaml* contains the parameters that can be adjusted by the user. Here we describe how to fill this file (a pre-filled configuration file is present is the ORFmine/examples/workdir/ directory, corresponding to the example dataset).
Please note that every file's name must be written without path and they must be placed in your local folder [linked](orfmine_binding.md) to the /database/ directory of the container. The fastq files must be stored in your local fastq/ directory that [has been linked](orfmine_binding.md) to the /fastq directory of the container when having launched the container.

Please place the config.yaml file in the */workdir/* directory of your container and edit it as follows:

### Prepare example configuration
An already filled configuration file is available : *ORFmine/examples/workdir/config.yaml*.
For more information on how to fill the configuration file, you can go [here](./orfribo_configuration.md).

### Running ORFribo

The only parameters to adjust in the command line are:

* The number of virtual CPUs / threads to reserve for the analysis

* The approximate amount of memory to reserve for the analysis (in GB)

ORFribo is launched like this :
``` bash
orfribo CPU MEMORY
```
For example, for an analysis launched with the use of 6 CPUs and around 10GB of memory, it would be :
``` bash
orfribo 6 10
```

Be careful to have enough memory on your computer/cluster for you analysis as a lack of RAM would lead the pipeline to crash without giving any error message. So if the program stops without any obvious reason this might be the reason why, and one should restart the pipeline with more resources.

<br><br>

### Main outputs















### Running ORFdate
ORFdate estimates the evolutionary ages of a set of ORFs (coding but also noncoding) based on phylostratigraphy. It takes as inputs:

* the name of the focal species as indicated in the tree (e.g. "Saccharomyces cerevisiae")
* a phylogenetic tree in newick format, containing the species of interest (focal species) and its neighbors (see ORFmine/examples/database/orfdate_inputs/ for an example)
* a two columns csv file (i.e. columns separated by commas) with the path to the sequence fasta files of each species as first column, and the associated species names as indicated in the phylogenetic tree as second column (see ORFmine/examples/database/orfdate_inputs/ for an example). It, therefore, involves the fasta files of the complete proteomes (amino acid sequences) of the focal and its neighbors to be stored in the corresponding directories
* optional arguments including BLAST parameters. See [here](./orfdate_run.md) for more details.

All the inputs must be stored in the /database directory, though a subdirectory can be created for more clarity (e.g. /database/orfdate_inputs/). The following command line will estimate the ages of 42 yeast CDSs (subsample of yeast CDSs that can be found in the file ORFmine/examples/database/orfdate_inputs/sample_42_Scer_CDS.pfasta) based on the phylogenetic tree of the Saccharomyces sensu stricto (i.e. <i>S. cerevisiae, S. paradoxus,  S. mikatae, S. kudrievezi, S. arboricola and S. bayanus </i>)(the corresponding tree can be found in ORFmine/examples/database/orfdate_inputs/saccharomycetaceae.nwk). Please note that last node of the tree will be used as the upper limit for the age estimation. Consequently, all sequences with a match in the farest species, with respect to the focal, are associated with the same upper bounded estimated age regardless of whether they may have other matches in more distant species outside the input tree (see Papadopoulos et al, [11] for more details).

```
orfdate -focal "Saccharomyces cerevisiae" -tree /database/orfdate_inputs/saccharomycetaceae.nwk -names /database/orfdate_inputs/name_correspondences.csv
```

The name_correspondences.csv file is available in the ORFmine/examples/database/orfdate_inputs/ directory as well as the 5 fasta files containing all the sequences of the 5 neighboring species of the focal.

ORFdates generates two outputs in csv format located in the /workdir/orfdate/ directory:

* *_hits.csv which reports the number of significant matches for all sequences of the focal species (rows) across the neighbors species of the tree (columns)
* *_dated.csv which reports for each sequence (col 1), the farest neighbor species for which a significant math has been found (col2), and the time of divergence between this neighbor and the focal species in Mya (col3)

with "*" being the basename of the focal fasta file. ORFdate generates tmp files (blast databases and blast outputs) which are stored in the /workdir/orfdate/blast_database/ and /workdir/orfdate/blast_out/ directories. Examples of such outputs can be found in the ORFmine/examples/workdir/orfdate/ directory.

<br><br><br>





#### References


1. Merkel, Dirk. "Docker: lightweight linux containers for consistent development and deployment." Linux j 239.2 (2014): 2.
2. Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459. Published 2017 May 11. doi:10.1371/journal.pone.0177459
3. Bitard-Feildel, T. & Callebaut, I. HCAtk and pyHCA: A Toolkit and Python API for the Hydrophobic Cluster Analysis of Protein Sequences. bioRxiv 249995 (2018).
4. Dosztanyi, Z., Csizmok, V., Tompa, P. & Simon, I. The pairwise energy content estimated from amino acid composition discriminates between folded and intrinsically unstructured proteins. Journal of molecular biology 347, 827–839 (2005).
5. Dosztányi, Z. Prediction of protein disorder based on IUPred. Protein Science 27, 331– 340 (2018).
6. Mészáros, B., Erdős, G. & Dosztányi, Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46, W329–W337 (2018).
7. Fernandez-Escamilla, A.-M., Rousseau, F., Schymkowitz, J. & Serrano, L. Prediction of sequence-dependent and mutational effects on the aggregation of peptides and proteins. Nature biotechnology 22, 1302–1306 (2004).
8. Linding, R., Schymkowitz, J., Rousseau, F., Diella, F. & Serrano, L. A comparative study of the relationship between protein structure and β-aggregation in globular and intrinsically disordered proteins. Journal of molecular biology 342, 345–353 (2004).
9. Rousseau, F., Schymkowitz, J. & Serrano, L. Protein aggregation and amyloidosis: confusion of the kinds? Current opinion in structural biology 16, 118–126 (2006).
10. Mészáros B, Erdős G, Dosztányi Z (2018) IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46:W329–W337
11. Papadopoulos, C., Arbes, H., Chevrollier, N., Blanchet, S., Cornu, D., Roginski, P., Rabier, C., Atia, S., Lespinet, O., Namy, O., Lopes, A. (submitted).
12. Lauria F, Tebaldi T, Bernabò P, Groen EJN, Gillingwater TH, Viero G. riboWaltz: Optimization of ribosome P-site positioning in ribosome profiling data. PLoS Comput Biol. 2018 Aug 13;14(8):e1006169.
