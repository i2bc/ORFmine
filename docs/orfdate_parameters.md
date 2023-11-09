## ORFdate parameters


<b>Mandatory</b>

```-focal FOCAL  ```        taxonomic name of the focal species <br>
```-names NAMES  ```        csv file matching fasta (col1) and tree (col2) names<br>
```-tree TREE   ```         newick file for the phylogeny tree<br>
 


<b>Optional</b>


  ```-h, --help```            show this help message and exit<br>
 ``` -ncpus NCPUS ```         total number of cpus that can be used fo the task (default: 1)<br>
```  -evalue EVALUE   ```     BLASTp evalue threshold (default: 1e-3)<br>
 ``` -query_cov QUERY_COV ``` minimum query coverage threshold (default: 0.70)<br>
 ``` -preserve_underscores ``` PRESERVE_UNDERSCORES (default: False). If True, 
                        the underscores that are present in the species names in the tree and csv file are considered explicitly and not replaced by spaces. Usually underscoeres in the tree labels are replaced by spaces. Usefull when the species name contains underscores (e.g. Scer_annotationV2). 