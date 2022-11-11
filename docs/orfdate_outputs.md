## ORFdate outputs:

ORFdate generates two outputs in csv format located in the /workdir/orfdate/ directory. Examples of these outputs can be found in the ORFmine/exmaples/workdir/orfdate/ folder.

* *_hits.csv which reports the number of significant matches for all sequences of the focal species (rows) across the neighbors species of the tree (columns)
* *_dated.csv which reports for each sequence (col 1), the farest neighbor species for which a significant math has been found (col2), and the time of divergence between this neighbor and the focal species in Mya (col3)

with "*" being the basename of the focal fasta file. ORFdate generates intermediate files (blast databases and blast outputs) which are stored in the /workdir/orfdate/blast_database/ and /workdir/orfdate/blast_out/ directories. Examples of such outputs can be found in the ORFmine/examples/workdir/orfdate/ directory.

