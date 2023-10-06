import argparse




def get_args():
    parser = argparse.ArgumentParser(description="Estimate the evolutionary age of ORFs in a genome based on phylostratigraphy using ORFdate.")

    parser.add_argument("--target", required=True, help="Taxonomic name of the focal species")
    parser.add_argument("--mapping", required=True, help="CSV file matching fasta (col1) and tree (col2) names", default=False)
    parser.add_argument("--tree", required=True, help="Newick file for the phylogeny tree")
    parser.add_argument("--out", required=False, type=str, default=".", help="Output directory (default='.')")
    parser.add_argument("--cpus", help="Total number of cpus that can be used fo the task (default=1)", default=1, type=int)
    parser.add_argument("--blast", help="Wether to perform BLASTp or not. True to perform BLASTp, False otherwise (default=True)", type=bool, default = True)
    parser.add_argument("--evalue", help="BLASTp evalue threshold (default=0.001)", default=0.001, type=float)
    parser.add_argument("--min-coverage", help="Minimum query coverage threshold (default=0.7)", default=0.7, type=int)
    parser.add_argument("--has-underscores", help="Wether underscores are kept when reading the tree, or considered as spaces. True to keep underscores, False otherwise (default=False)", default = False, type=bool)
    parser.add_argument('--keep-files', action='store_true', default=False, help="Add '--keep' to keep intermediary computed files such as blastdb and blast outputs")

    container_group = parser.add_mutually_exclusive_group()
    container_group.add_argument("--docker", action='store_true', default=False, help="Flag used to run computations on a docker container")
    container_group.add_argument("--singularity", action='store_true', default=False, help="Flag used to run computations on a singularity container")

    return parser.parse_args()
