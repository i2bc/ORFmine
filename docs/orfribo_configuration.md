## Usage


### Prepare your parameters in the configuration file

In case you plan on using ORFribo, you must adjust all the parameters in the configuration file *config.yaml* which you can find in the *orfribo/* folder on GitHub [here] (https://github.com/i2bc/ORFmine/orfribo/config.yaml).

In this file, every file name must be written without path and these files must be placed in the folder linked to the */database/* directory in the container.

#### Project name  
**project_name**: the name you want for your project (please avoid blank spaces or special characters)

#### Name of database subfolder file  
You must enter the full name (**with extensions**) without the path of files added in the database subfolder previously created.   

**fasta**: reference_genome_fasta_file.fa  

**gff**: corresponding_gff_annotation_file.gff3  

**gff_intergenic**: the name of ORFtrack output fasta file

**fasta_outRNA**: unwanted_sequences_fasta_file.fa  

#### Pipeline option selection  
During the ORFribo process, data is trimmed and selected depending on the read lengths.

**already_trimmed**: If your data contains reads already trimmed of their adapter, you can set this option on “yes”. Else, set it on "no".

**adapt_sequence**: If they are not trimmed, you should specify the sequence of the adapter in quotes on the line here like "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA". If you do not put anything between the quotes, RiboDoc will try to find the adapter itself but this can sometimes lead to a wrong adapter sequence.

#### Reads length to analyse
You also have to define the range for read lengths selection. Default values select reads from 25 to 35 bases long.  

**readsLength_min**: minimum read length.   

**readsLength_max**: maximum read length.   

#### Names of coding features in GFF
You might also need to specify features keywords in the GFF file to fit your GFF file format :   

**gff_cds_feature**: Feature corresponding to CDS in the annotation file. "CDS" is the default value (can sometimes be "ORF").     

**gff_name_attribut**: Name of the genes features in the GFF. Default is "Name" but it can sometimes be "gene_name" or else.     

#### Statistical settings  
**orfstats_mean_threshold**: Minimum mean of in-frame reads in coding region to select specific read lengths (default is 70)

**orfstats_median_threshold**: Minimum median of in-frame reads in coding region to select specific read lengths (default is 70)

#### Names of features or ORF categories to be analyzed (noncoding but also coding)
**final_counts**: List of ORF categories for which you want to investigate the translation activity. The ORF categories listed here must correspond to those provided by ORFtrack in the 3rd column of the output gff file or indicated in ORFtrack's *summary.log*. Examples of these two files can be found in the ORFmine/examples/database/ directory as mapping_orf_Scer.gff or summary.log files. By default, ORFribo will probe the translation activity of all intergenic ORFs which are referred to as "nc_intergenic" (see [here](./orftrack_annotation.md) for more details on the ORF categories and annotation process). If you want to also probe that of noncoding ORFs lying in the alternative frames of CDSs on the same strand, one should write their names in a single quoted block, each feature separated by a single blank space.
For example, if you want to also probe that of noncoding ORFs lying in the alternative frames of CDSs on the same strand, you must add the "nc_ovp_same-CDS", like this : "nc_intergenic nc_ovp_same-CDS".


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
