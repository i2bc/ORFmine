#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 14:40:40 2022

@author: paul.roginski
"""
# TODO : determine the number of cpu to use per blast
# TODO : gérer les fichiers de correspondance incomplets

import argparse
from ast import Dict, List
import pandas as pd
import numpy as np
import multiprocessing
import dendropy
from functools import partial
import operator
from pathlib import Path
import os
import shutil
import sys

from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastpCommandline



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-focal", required=True, help="Taxonomic name of the focal species")
    parser.add_argument("-names", required=True, help="CSV file matching fasta (col1) and tree (col2) names", default=False)
    parser.add_argument("-tree", required=True, help="Newick file for the phylogeny tree")
    parser.add_argument("-out", required=False, type=str, default="./", help="Output directory (default='./')")
    parser.add_argument("-ncpus", help="Total number of cpus that can be used fo the task (default=1)", default=1, type=int)
    parser.add_argument("-blast", help="Wether to perform BLASTp or not. True to perform BLASTp, False otherwise (default=True)", type=bool, default = True)
    parser.add_argument("-evalue", help="BLASTp evalue threshold (default=0.001)", default=0.001, type=float)
    parser.add_argument("-query_cov", help="Minimum query coverage threshold (default=0.7)", default=0.7, type=int)
    parser.add_argument("-preserve_underscores", help="Wether underscores are kept when reading the tree, or considered as spaces. True to keep underscores, False otherwise (default=False)", default = False, type=bool)
    parser.add_argument('--keep', action='store_true', default=False, dest='bool_keep', help="Add '--keep' to keep intermediary computed files such as blastdb and blast outputs")

    return parser.parse_args()


def generate_tree(tree_file: str, names_df: pd.DataFrame, preserve_underscores: bool) -> dendropy.Tree:
    """Generate and print a tree object that will be used to infer hits distance.

    Args:
        tree_file (str): Path to a phylogenetic tree filename (.newick format) 
        names_df (pd.DataFrame): pandas 2D dataframe

    Returns:
        dendropy.Tree: {}
    """.format(dendropy.Tree)

    tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=preserve_underscores)
    print("Original tree")
    print(tree.as_ascii_plot())

    # remove from tree taxa absent in names_df 
    taxa = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]
    extra_taxa = [ label for label in taxa if label not in names_df['tree'].tolist()]

    if len(extra_taxa) > 0 :
        print("Extra taxa are : {}".format(extra_taxa))
        tree.prune_taxa_with_labels(extra_taxa)

        # Also need to correct the "taxon_namespace" attribute
        tree.taxon_namespace = [ tree.taxon_namespace[i] for i in range(0, len(tree.taxon_namespace)) if tree.taxon_namespace[i].label in names_df['tree'].tolist()]
        print("Corrected tree")
        print(tree.as_ascii_plot())

    return tree


def map_pyhlo_distance(tree: dendropy.Tree, focal_name: str) -> Dict:
    """ Map phylogenetic distance relative to a species corresponding to a given focal name

    Args:
        tree (dendropy.Tree): instance of a dendropy.Tree arborescence 
        focal_name (str): Name of the focal species (e.g. "Scer SGD ORF") 

    Returns:
        Dict: dictionnary with species as keys and distance to the focal species as values
    """

    # Retrieve in the taxon_namespace attribute of the focal species.
    try :
        focal = tree.taxon_namespace.get_taxon(label=focal_name)
    except :
        labels = [tree.taxon_namespace[i].label for i in range(len(tree.taxon_namespace))]
        print("Cannot find the -focal {} in the names of the tree :\n{}".format(focal_name, labels))
        sys.exit(1)

    pdc = tree.phylogenetic_distance_matrix()  # distance matrix object

    # Build a dictionnary with species as keys and distance to the focal species as values
    distance_to_focal = { t.label:pdc(focal, t) for t in tree.taxon_namespace }

    if sum(distance_to_focal.values()) == 0 :
        print("No distance is specified in the newick file. Exiting.")

    return dict(sorted(distance_to_focal.items(), key=operator.itemgetter(1)))  # Sorted version


def blastp(fasta, focal_fasta, names_df, out_path, num_threads=4, min_cov=0.7, evalue=0.001, is_blast=True):

        database_folder = Path(out_path) / "blastdb"
        blastout_folder = Path(out_path) / "blastout"

        # Taxon name
        tree_name = names_df['tree'][np.where(names_df['fasta'] == fasta )[0][0]]
        print("{}...".format(tree_name))

        # BLASTp output file's name
        output_file = blastout_folder / str("{}_vs_{}.out".format(Path(focal_fasta).stem, Path(fasta).stem))

        if is_blast :
            db = database_folder / Path(fasta).stem

            # Build a protein BLAST database from the subject fasta
            cline = NcbimakeblastdbCommandline(dbtype="prot", input_file=fasta, out=db)
            cline()

            # Perfom a local BLASTp of the focal fasta against the current fasta's database
            cline = NcbiblastpCommandline(
                query=focal_fasta,
                db=db,
                out=output_file,
                outfmt="6 qseqid sseqid evalue qlen qstart qend slen sstart send length bitscore score",
                evalue=10,
                num_threads=num_threads
            )
            cline()


        # Parse the BLASTp output file and count the number of hits for each query ORF
        hits = {}
        with open(output_file, 'r') as f:
            for line in f.readlines() :
                if line[0] != "#":
                    full_line = line.split("\t")

                    # Coverage of the query = (q.stop - q.start) / query length
                    query_cov = (int(full_line[5]) - int(full_line[4])) / int(full_line[3])

                    # evalue is at position 3 (index 2)
                    if float(full_line[2]) < evalue and query_cov >= min_cov :
                        query = full_line[0]
                        if query in hits.keys() :
                            hits[query] += 1
                        else :
                            hits[query] = 1

        # The returned object will be merged with others into a dataframe
        return [tree_name, hits]


def parse_blast_hits(hits: List) -> pd.DataFrame:
    """Build a dataframe from a set of blastp hits.
    This dataframe, once transposed and renamed
    will have focal ORFs as rows and tree's taxa as columns.

    Args:
        hits (List): List of Dict

    Returns:
        pd.DataFrame: pandas 2D dataframe with ORFs as rows and tree's taxa as columns.
    """
    print("Merging all hits in one dataframe...")
    hits_df = pd.DataFrame([ele[1] for ele in hits])
    hits_df = hits_df.transpose()
    hits_df.columns = [ele[0] for ele in hits]
    hits_df = hits_df.fillna(0)
    hits_df = hits_df.astype(int)

    return hits_df


def write_outputs(hits: Dict, distances: Dict, out_basename: str) -> None:

    # re-order the columns of hits_df according to their distance to the focal species
    hits = hits[list(distances.keys())]

    # Write hits to a csv file for the user
    print("Writing the hits tableau...")
    hits.to_csv("{}_hits.csv".format(out_basename), index=True, index_label='seq')


    # Generate an empty dataframe that will be filled with :
    # first column = the list of the farthest species with a hit
    # second column = the distance of these species to the focal species.
    # print("Building the output csv...")
    farest_df = pd.DataFrame(index=hits.index, columns = ["farest_hit","distance"])

    # For each ORF
    for orf in hits.index :        
        row = hits.loc[orf,:]  # Its corresponding row in hits        
        hit_names = [name for name in hits.columns if row[name] > 0 ]  # Species-with-hits's names        
        hit_distances = [ distances[hit_name] for hit_name in hit_names ]  # Their distance to the focal species        
        max_distance = max(hit_distances)  # The maximum distance to the focal species        
        farest_names = [name for name in hit_names if distances[name] == max_distance ]  # Names of the species with hits that are at this maximum distance
        farest_df.loc[orf,"farest_hit"] = "|".join(farest_names)
        farest_df.loc[orf,"distance"] = max_distance

    # Write hits_df to a csv file for the user
    print("Writing the output csv...")
    farest_df.to_csv("{}_dated.csv".format(out_basename), index=True, index_label='seq')


def main():
    args = get_args()
    
    focal_name = args.focal  # taxonomic name of the focal species
    tree_file = args.tree  # newick file for the phylogeny tree
    ncpus = args.ncpus  # Number of cpus to use
    is_blast = args.blast  # boolean flag to inform on wether to perform blastp or not 
    evalue = args.evalue  # BLASTp evalue threshold    
    min_cov = args.query_cov  # minimum query coverage threshold
    out_path = Path(args.out)  # main output directory
    keep_files = args.bool_keep  # boolean flag to inform on wether to keep intermediary files or not
    preserve_underscores = args.preserve_underscores

    # Dataframe with fasta files' names as first column and tree taxons' names as second column
    names_df = pd.read_csv(args.names, names = ["fasta","tree"])

    # fasta name of the focal species
    try :
        focal_fasta = names_df['fasta'][np.where(names_df['tree'] == focal_name )[0][0]]
    except :
        print("The focal name {} is not found.\nNames provided in -names :\n{}".format(focal_name, names_df['tree']))
        sys.exit(1)


    # Generate and print a tree object that will be used to infer hits distance.
    tree = generate_tree(tree_file=tree_file, names_df=names_df, preserve_underscores=preserve_underscores)

    # get phylogenetic distances
    distance_to_focal = map_pyhlo_distance(tree=tree, focal_name=focal_name)

    if is_blast :
        print("Performing BLASTp on provided fastas...")
        database_folder = Path(out_path) / "blastdb"
        database_folder.mkdir(parents=True, exist_ok=True)

        blastout_folder = Path(out_path) / "blastout"
        blastout_folder.mkdir(parents=True, exist_ok=True)

    else :
        print("Skipping BLASTp.")
        
    # Perform the blastp function on each fasta
    fastas = [fasta for fasta in names_df['fasta']]
    with multiprocessing.Pool(processes=ncpus) as pool:
        all_hits = pool.map(partial(blastp, focal_fasta=focal_fasta, names_df=names_df, out_path=str(out_path), min_cov=min_cov, evalue=evalue, is_blast=is_blast), fastas)
        
    # all_hits = [blastp(fasta, focal_fasta, names_df, out_path, ncpus=ncpus, min_cov=min_cov, evalue=evalue, is_blast=is_blast) for fasta in names_df['fasta']]

    if not keep_files:
        shutil.rmtree(blastout_folder)
        shutil.rmtree(database_folder)

    # Build a dataframe from all of blastp hits
    hits_df = parse_blast_hits(hits=all_hits)

    # write outputs
    write_outputs(hits=hits_df, distances=distance_to_focal, out_basename=str(out_path / Path(focal_fasta).stem))
