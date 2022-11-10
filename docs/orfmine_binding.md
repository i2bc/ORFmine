### Binding the local folders to the container

When you will create your container with the singularity or docker commands, you will bind the folders you have prepared on your local machine to the folders already created inside the container with the *-B* or *-v* arguments (for singularity or docker respectively) by separting them with a colon in the command like this:
`/path/on/your/local/machine/:/corresponding/path/inside/the/container`

You must precise the exact local folder you want to bind, not its parent folder.

* workdir (mandatory):
    */path/to/working/folder/* linked to **/workdir/**  
    This is needed everytime you start the container since it links your working directory with that of the container (named "workdir").
    `/path/on/your/working/folder/:/workdir/`
* database (mandatory):
    */path/to/data/folder/* linked to **/database/**  
    This is needed everytime you start the container since it links the folder containing your data (gff, fasta files...) to the "database" folder inside the container.
    `/path/on/your/data/folder/:/database/`
* fastq (only if you use ORFribo):
    */path/to/fastq/folder/* linked to **/fastq/**  
    Needed when using ORFribo in order to indicate where the Ribosome Profilng (i.e. fastq.gz files) are stored.
    This folder will be linked to the "fastq" folder inside the container.
     `/path/on/your/fastq/folder/:/fastq/`



**Examples:**

Let's say all the genome sequence and annotation files are stored in a directory named "Genomes" and located in the "/home/ProjectXYZ/" directory. This will be indicated in the singularity or docker command lines through the -B/-v argument (depending on the use of docker (-v) or singularity (-B) respectively) as follows (see below for the complete command lines):

`-v /home/ProjectXYZ/Genomes/:/database/` for docker
`-B /home/ProjectXYZ/Genomes/:/database/` for singularity

If the Ribosome Profiling fastq files are stored in the folder "WT_vs_Mutant" located in the "/home/Ribosome-profiling/October2022/" directory, it will be indicated in the singularity or docker command lines as follows (see below for the complete command lines):

`-v /home/Ribosome-profiling/October2022/WT_vs_Mutant/:/fastq/` for docker
`-B /home/Ribosome-profiling/October2022/WT_vs_Mutant/:/fastq/` for singularity


To exit the container just type :
``` bash
exit
```
