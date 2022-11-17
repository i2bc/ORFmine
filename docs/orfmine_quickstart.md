
# Quick start

Here we first explain how to [prepare the folders](#prepare_folders) necessary for ORFmine and [launch the container](#start_container).



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
lopesi2bc/orfmine:latest /bin/bash</code>
</pre>


Please note that if you do not want to use ORFribo, you do not need to link the fastq folder to the container and can remove the binding part of the command line: <code><i>PATH_TO_YOUR_FASTQ_DIR</i>:/fastq/</code>. See [here](./orfmine_binding.md) for more details on how to bind your local folders to those of the container.


#### Binding Tango and/or IUPred

If you want to use Tango and/or IUPred with ORFfold, you will need to have those softwares installed on your local computer (choose the Unix compatible source code for Tango). 

Once installed, please place the parent directory content of those softwares inside a directory of your choice. Below is a description of the workflow we followed:

```bash
# let's consider Tango is installed in ~\tango and IUPred in ~/iupred2a

# create a directory that will contain Tango and IUPred source codes content
mkdir ~/softwares

# place Tango UNIX source code content in softwares/
cp ~/tango/tango_x86_64_release ~/softwares

# place IUPred necessary codes content in softwares/
cp ~/iupred2a/iupred2a_lib.py ~/softwares
cp ~/iupred2a/iupred2a.py ~/softwares
cp -r ~/iupred2a/data ~/softwares
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
lopesi2bc/orfmine:latest /bin/bash</code>
</pre>




You are ready to run your analyses! That said, we strongly recommend starting with a [small example](./orfmine_example.md) we prepared for you.

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">

    <ul>
    <li> When launching the container with docker, the <i>lopesi2bc/orfmine:latest</i> part can be changed if you want to launch a specific version of ORFmine. For example : <i>lopesi2bc/orfmine:v0.8.7</i>.
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

