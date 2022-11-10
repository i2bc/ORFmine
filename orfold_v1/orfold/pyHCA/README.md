![Python 3.5](https://img.shields.io/badge/Python-3.5-blue.svg) ![Python 3.6](https://img.shields.io/badge/Python-3.6-blue.svg) ![travis-ci](https://travis-ci.org/T-B-F/pyHCA.svg?branch=master)

pyHCA and hcatk (HCA toolkit) are a python library and executable for the Hydrophobic Cluster Analysis of protein sequences.
pyHCA implements various class and function for handling sequences and analyses them
hcatk provides off the hands set of program to perform hydrophobic cluster analyses.

Requires
========

- python3      >= 3.6
- Biopython    >= 1.65
- scipy        >= 1.0
- numpy        >= 1.14
- scikit-learn >= 0.19
- requests     >= 2.18
- ete3         >= 3.1.1


Installation
============

A quick install can be perform using:

    pip3 install .
 

The 'domOnTree' functionality of **hcatk** requires that ete3, and therefore PyQt4, are installed.
However, the ete3 can be difficult to install as some features requires PyQt4 and sip.
Please refer to the official ete3 insallation guidelines ![http://etetoolkit.org/download/](http://etetoolkit.org/download/) for any support.
On Mac OS X, you will also need to install XQuartz to use ete3, please refer to ![XQuartz documentation](http://www.xquartz.org/).

If you do not plan to use this feature you can proceed to the installation without PyQt4 and ete3.

We recommend you to work on a conda virtual environment to properly build the non Python extention of the ete3 package and afterward install pyHCA in this new environment.


A Dockerfile is also available to install the executable inside a container:

    docker build -t hcatk .
    docker run hcatk -h
    docker run --rm -v /path/to/data:/data hcatk draw -i /data/orc1.fasta -o /data/orc1_output.svg


Currently ete3/PyQt4 don't seem to be properly installed using this method, therefore the 'domOnTree' is not available if hcatk is built using docker.


Example
*******

download and install conda from miniconda website using the correct installer (64 bits / 32 bits, MacOSC / Linux):

on MacOSX

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

on Linux (64 bits installer)

    https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

create a virtual environment and switch to the environment

    conda create -n test_pyHCA python=3.6 pip
    source activate test_pyHCA
    
and install pyqt4 before running pyHCA installer

    conda install -n test_pyHCA pyqt=4.11.4
    cd <path to pyHCA directory>
    pip install .
    


Usage
=====


dispred
--------

Predict protein disorder at the amino acid level using HCA and physico-chemical features from AAIndex

    $ hcatk dispred -h

    usage: hcatk dispred [-h] -i FASTAFILE -o OUTPUTFILE -m MODEL [--verbose]


Arguments:
**********

    -h, --help            show help message and exit
    
required arguments:

    -i FASTAFILE   the fasta file
    -o OUTPUTFILE  output file in svg format
    -m MODEL       model to use

optional arguments:

    -v --verbose      print information


Example:
********

    $ hcatk dispred -i data/orc1.fasta -o data/orc1.allowDSoverlap.txt -m pyHCA/model_weights/allowDSoverlap.h5


Models `allowDSoverlap.h5` and `noDSoverlap.h5` can be found in the `pyHCA/model_weights/` directory of the github repository.


Output format:
**************

The output is formated in a fasta like file, with an header storing the protein name

    >protein_name

followed by four columns (index, amino acid, score, prediction: 0=ordered, 1=disordered)

	1 M 0.202 0
	2 V 0.223 0
	3 N 0.270 1
	4 K 0.292 1
	5 E 0.353 1
	6 N 0.450 1
	7 A 0.312 1


segment
-------

Use the composition in hydrophobic cluster of a sequence to detect domains.

    $ hcatk segment -h

    usage: hcatk segment [-h] -i INPUTF -o OUTPUTF [-v]
                          [-m {cluster,domain}]
                          [-t {aminoacid,nucleotide}]

Arguments:
**********

    -h, --help            show help message and exit
    
required arguments:

    -i INPUTF             an amino-acid sequence files in fasta format
    -o OUTPUTF            the output file with annotation
    
optional arguments:

    -v                    keep temporary results
    -m {cluster,domain}   method to use, cluster: will report *the hydrophobic
                          clusters found in the sequence, domain: will delineate
                          domains based on the hydrophobic cluster profile of
                          the sequence
    -t {aminoacid,nucleotide}
                          the type of the biological sequences passed in the
                          input file

Example:
********

    $ hcatk segment -i data/orc1.fasta -o data/orc1.hca -m domain


Output Format:
**************

The output is formated in a fasta like style, with an header storing the protein name, size and HCA score over the full protein sequence:

    >protein_name protein_size protein_hca_score
    
followed by four columns:

    domain  524     527     nan 9.5
    domain  552     923     0.0032921164246364487 2.3978494623655915
    cluster 1       2       11
    cluster 10      17      10001011
    cluster 23      23      1
    
The first column correspond to the hca element identified, either a domain or a cluster.
The second and third columns correspond to the start and stop (indexed from 1 to the sequence length) of the hca element.
The fourth column either corresponds to a p-value if the element is a hca domain or to the hydrophobic cluster in binary mode if the element is a cluster.
The p-value of the domain is computed against a reference distribution made of intrinsically disordered proteins and describe the "degree of foldability" associated to the hca domain element.
A nan value is returned if the domain has less than 30 residues.
The fourth column corresponds to the HCA score of the domain sequence.

/!\ Warning /!\

HCA scores are provided for the whole protein sequence and domain with less than 30 amino acids for information only.
The scores and p-values associated where designed on globular protein domains and intrinsically disordered regions of comparable lengths and with more than 30 amino acids.



    
draw
----

Draw the HCA diagram of a sequence.
Optionnaly, can display the domain annotation of a sequence if provied.

    $ hcatk draw -h

    
    usage: hcatk draw [-h] -i FASTAFILE [-w WINDOW] [-d DOMAIN] [-f {pfam,seghca}]
                  [--color-msa {rainbow,identity}] -o OUTPUTFILE

Arguments:
**********

    -h, --help            show help message and exit
    
required arguments:
    
    -i FASTAFILE          the fasta file
    -o OUTPUTFILE         output file in svg format

optional arguments:

    -w WINDOW             sequence len before breaking the sequence to the next
                          plot (-1 the whole sequence are used, minimum size is
                          80)
    -d DOMAIN             [optionnal] provide domain annoation
    -f {pfam,seghca}      the domain file format
    --cons-msa {aa,hca}   method to use to compare sequences (aa, hca)

Example:
********
    
    $ hcatk draw -i data/PF00533_sub.txt -o data/PF00533_sub.svg --cons-msa aa

Output Format:
**************    

The output is a svg text file and can be vizualised using any svg viewer (inkscape, modern web browser ...).
    

tremolo
-------

Use Tremolo-HCA to find remote homologous proteins with domain context.

    $ hcatk tremolo -h

    usage: hcatk [-h] -i INPUTFASTA [-d DOMAINS [DOMAINS ...]] 
                 [-w WORKDIR] [--pp2ipr P2IPR] [--cpnfig CONFIGFILE] 
                 [-o OUTPUT] --target-db


Arguments:
**********

    -h, --help            show help message and exit

required arguments:

    -i INPUTFASTA         input fasta file
    -o OUTPUT             output file
    -w WORKDIR            working directory

optional arguments:

    -d DOMAINS [DOMAINS ...] list of domain positions (start and stop inclusive
                             and separated by comma : -d 1,10 20,30 60,100. 
                             If not provided the search will be performed on 
                             each domain found after segmentation of the input 
                             sequence. To use the whole protein use -d all.
    --p2ipr P2IPR            path to the Interpro annotation of UniproKB proteins,
                             gzip format supported.
    -config CONFIGFILE       path to the configuration file for optional arguments
    --target-db HTARGETDB    path to the target sequences database

Example:
********
    
    $ echo "using hhblits"
    $ hcatk tremolo -i data/orc1.fasta -d 10,143 --p2ipr data/protein2ipr.dat.gz --config tremolo_config.ini --target-db hhsuite/uniprot20_2016_02/uniprot20_2016_02 -o data/orc1_tremolo.txt -w tremolo_tmp


Output Format:
**************    

Tremolo output format is made of multiple parts.

A header part storing the protein sequence information used as query:

    Qname   protein_sequence_name
    Qdesc   protein_description
    Qseq    protein_sequence

    Qdom    domain_index       domain_start      domain_stop

If multiple domains were given as query (through the -d option), multiple Qdom fields will be present:

    Qdom    1       10      50
    Qdom    2       60      80
    
The header is followed by a two summaries of the results.
The first summary give for each hit of each query domain, the name of the protein target, the id of the hit (in case of multiple hits per target) , the domain arrangement of the protein target (with interpro ID and domain specific ID spearated by a "/"), the e-value and the bit-score of the hit.
The second summary give some information on the protein domain arrangement diversity associated with the domain query.
The second summary fields start with "INFO" and look like:

    INFO    1       IPR001025/PF01426;IPR001025/SM00439;IPR001025/PS51038   11
    
After the header and the summary, a detailled description of the results are provided per domain query sorted by domain arrangement and protein target:

    ## <- start of a domain arrangement results
    >protein_target protein_size
    Qdom    1       IPR001025/PS51038       167     286      <- Qdom <query domain index> <interpro domain/db specific name> <start> <stop>
    Qdom    1       IPR001025/SM00439       167     286      <- Qdom <query domain index> <interpro domain/db specific name> <start> <stop>
    Qdom    1       IPR001025/PF01426       168     261      <- ...
    Qdom    1       IPR008395/PF05641       386     459     
    Qdom    1       IPR014002/SM00743       466     524     
    Hit     1       0       7e-06   97.16   79.55   17.0    0.195 <- Hit <query domain index> <hit index, starts from 0 and increase if multiple hits> <HHblits e-value> <HHblits proba score> <HHblits bit-score> <identity> <similarity>
    HitQali 1       0       22      134     HKNVYFYQKCIYGPLTLSVGDFILVSNADAAE ... <- HitQali <query domain index> <hit index> <hit start> <hit stop> <sequence>
    HitQcon 1       0       22      134     ~k~~~fy~kc~~~~~~i~vGdfVLIen~D~ae ...
    HitTcon 1       0       154     248     gKqLkHYpsFcRNGtTI~VqSFVfVMake--- ...
    HitTali 1       0       154     248     GKQLKHYPSFCRNGTTISVQSFVFVMSKE--- ...
    
    //  <---- End of a hit
    ## <---- End / beginning of domain arrangement
    
    >tr|V5HB40|V5HB40_IXORI 420
    Qdom    1       IPR001025/PS51038       1       103
    Hit     1       0       7.9e-06 97.21   74.04   22.0    0.412
    HitQali 1       0       77      134     REPCRAIVQWYSWPKAIPHNKYDDDEVAIDF ...
    HitQcon 1       0       77      134     ~~~krA~VQWfsR~~eiP~~kr~ll~r~~~~ ...
    HitTcon 1       0       7       63      kdhrfvtvqwylrvtelpptqqgrlghcdyf ...
    HitTali 1       0       7       63      KDHRFVTVQWYLRVTELPPTQQGRLGHCDYF ...
    
    // <---- End of a hit but next target hit is a different protein but with the same domain arrangement so no double # symbol
    
    >tr|H2ZNF4|H2ZNF4_CIOSA 222
    Qdom    1       IPR001025/PS51038       53      176
    Hit     1       0       3.7e-05 96.89   65.49   21.0    0.288
    HitQali 1       0       36      133     LTLSVGDFILVSNADAAEPDTVSGCDVARIL ...
    HitQcon 1       0       36      133     ~~i~vGdfVLIen~D~aepd~~d~~~VAki~ ...
    HitTcon 1       0       54      139     nlisigdgvviacges-----kqdfylaqvs ...
    HitTali 1       0       54      139     NLISIGDGVVIACGES-----KQDFYLAQVS ...
    //
    ...


domOnTree
---------

Vizualise Tremolo-HCA results on a taxonomic tree with protein domain arrangement information


    $ hcatk domOnTree -h

    usage: hcatk [-h] [-i TREMOLORES] [-t TREEFILE] [-s PROT2SPECIES]
                 [-n NCBITAXID [NCBITAXID ...]] [-o OUTPUT]


Arguments:
**********

    -h, --help            show help message and exit
    
required arguments:

    -i TREMOLORES         tremolo results with domain matchs
    -o OUTPUT             phylogenetic tree with tremolo hits
    
optional arguments:


    -t TREEFILE           phylogenetic tree with node as ncbi taxonomic ids
    -s PROT2SPECIES       file with prot to species informations
    -n NCBITAXID [NCBITAXID ...]
                          list of node for which leaves will be merged (internal
                          node need to be in tree)


Example:
********

    $ hcatk domOnTree -i data/tremolo_result.txt -o data/tremolo_result.pdf
                          
Output Format:
**************    

The output is a pdf file representing the domain arrangement of each sequence associated with the protein domain query used in tremolo. 
Each sequence is positioned on a taxonomic tree according to the species to which the sequence belongs to.


Additional ressources
---------------------

The interpo domain annoation can be downloaded at:
    wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/protein2ipr.dat.gz

HHblits of the HH-suite package can be downlad at (v3 or higher):
    git clone git@github.com:soedinglab/hh-suite.git

And the uniprot hhblits compatible database at:
    http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
