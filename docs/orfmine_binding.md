

### Binding the local folders to the container

When creating your container with the singularity or docker commands, the folders you have prepared on your local machine are linked to the folders already created inside
the container through the *-B* or *-v* arguments (for singularity or docker
respectively). The paths of the local and container folders are separated with a colon as follows:
`/path/on/your/local/machine/:/corresponding/path/inside/the/container`

Pleae note that you must indicate the complete path of your local folder, not its parent folder.

* workdir (mandatory):
    */path/of/your/local/working/folder/* to be linked to **/workdir/**  
    This is needed everytime you
    start the container since it links your working directory with that of the container (named "workdir").
    `/path/of/your/local/working/folder/:/workdir/`
* database (mandatory):
    */path/of/your/data/folder/* to be linked to **/database/**  
    This is needed everytime you
    start the container since it links the folder containing your data (gff, fasta files...) to the
    "database" folder of the container.
    `/path/of/your/local/data/folder/:/database/`
* fastq (only if you use ORFribo):
    */path/of/your/local/fastq/folder/* to be linked to **/fastq/**  
    Needed when using ORFribo in order to indicate where the Ribosome Profilng data (i.e. fastq.gz files) are stored.
    This folder will be linked to the "fastq" folder of the container.
     `/path/of/your/local/fastq/folder/:/fastq/`



**Examples:**

Let's say all the genome sequence and annotation files are stored in a directory
named "Genomes" and located in the "/home/ProjectXYZ/" directory. This will be
indicated in the singularity or docker command lines through the
-B/-v argument (depending on the use of docker (-v) or singularity (-B)
respectively) as follows:

`-v /home/ProjectXYZ/Genomes/:/database/` for docker
`-B /home/ProjectXYZ/Genomes/:/database/` for singularity

If the Ribosome Profiling fastq files are stored in the folder
"WT_vs_Mutant" located in the "/home/Ribosome-profiling/October2022/"
directory, it will be indicated in the singularity or docker command lines as follows:

`-v /home/Ribosome-profiling/October2022/WT_vs_Mutant/:/fastq/` for docker
`-B /home/Ribosome-profiling/October2022/WT_vs_Mutant/:/fastq/` for singularity

Examples of the comlete command line can be found [here](./orfmine_quickstart.md#launch-the-container).
