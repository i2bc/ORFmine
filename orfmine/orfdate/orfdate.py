#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 14:40:40 2022

@author: paul.roginski
"""
# TODO : determine the number of cpu to use per blast
# TODO : gérer les fichiers de correspondance incomplets

# Standard library imports
from datetime import datetime
from functools import partial
import multiprocessing
from pathlib import Path
import operator
import shutil
import sys
from typing import Dict, List, Union

# Third-party imports
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
import dendropy
import numpy as np
import pandas as pd
from orfmine import DOCKER_IMAGE

# Self-party imports
from orfmine.orfdate.lib import arguments
from orfmine.utilities.container import ContainerCLI


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
    extra_taxa = [ label for label in taxa if label not in names_df['taxon'].tolist()]

    if len(extra_taxa) > 0 :
        print("Extra taxa are : {}".format(extra_taxa))
        tree.prune_taxa_with_labels(extra_taxa)

        # Also need to correct the "taxon_namespace" attribute
        tree.taxon_namespace = [ tree.taxon_namespace[i] for i in range(0, len(tree.taxon_namespace)) if tree.taxon_namespace[i].label in names_df['taxon'].tolist()]
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
        print("Cannot find the target taxon {} in the names of the tree :\n{}".format(focal_name, labels))
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
        tree_name = names_df['taxon'][np.where(names_df['fasta'] == fasta )[0][0]]
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


def remap_to_dest(df: pd.DataFrame, colname: str, dest: str):
    df[colname] = df[colname].apply(lambda x: str(Path(dest) / Path(x).name))

    return df

def get_mapping_df(csv_file: Union[str, Path], colnames: List[str]):
    # Dataframe with fasta files' names as first column and tree taxons' names as second column
    mapping_df = None
    try:
        mapping_df = pd.read_csv(csv_file, names=colnames)
    except FileNotFoundError:
        print(f"Error: unable to find the csv mapping file: {csv_file}")
        print("\nPlease ensure the path exists.")
    except pd.errors.EmptyDataError:
        print(f"Error: the csv mapping file is empty: {csv_file}")
    except pd.errors.ParserError as e:
        print(f"Error while parsing the csv mapping file: {csv_file}")
        print("Details:", str(e))

    # Check if mapping_df was assigned and is not empty
    if mapping_df is None or mapping_df.empty:
        print("Error: The DataFrame was not properly initialized.")
        exit(1)

    return mapping_df


def check_files(files: Union[str, Path]):
    files_not_found = False
    for _file in files:
        if not Path(_file).resolve().exists():
            print(f"Error: unable to find {_file}")
            files_not_found = True

    if files_not_found:
        print(f"\nPlease check your mapping csv file and ensure all the paths are correct.")
        exit(1)


def run_orfdate(
        target: str, newick_file: Union[str, Path], out_path: Union[str, Path], mapping_file: Union[str, Path],
        ncpus: int=1, is_blast: bool=True, evalue: float=0.001, min_cov: float=0.7,
        keep_files: bool=False, preserve_underscores: bool=False, is_container: bool=False):
    
    # read csv file to get a df mapping fasta files and their user-defined reference name
    mapping_df = get_mapping_df(csv_file=mapping_file, colnames=["fasta", "taxon"])

    # if run on the container, path in the csv files must match their expected location on the container;
    # their location has been set in /input on the container filesystem
    if is_container:
        remap_to_dest(df=mapping_df, colname="fasta", dest="/input")
    
    # check if files exist
    check_files(files=mapping_df["fasta"].to_list())

    # fasta name of the focal species
    try :
        focal_fasta = mapping_df['fasta'][np.where(mapping_df['taxon'] == target )[0][0]]
    except :
        print("The focal name {} is not found.\nNames provided in --mapping :\n{}".format(target, mapping_df['taxon']))
        sys.exit(1)

    # Generate and print a tree object that will be used to infer hits distance.
    tree = generate_tree(tree_file=newick_file, names_df=mapping_df, preserve_underscores=preserve_underscores)

    # get phylogenetic distances
    distance_to_focal = map_pyhlo_distance(tree=tree, focal_name=target)

    if is_blast :
        print("Performing BLASTp on provided fastas...")
        database_folder = Path(out_path) / "blastdb"
        database_folder.mkdir(parents=True, exist_ok=True)

        blastout_folder = Path(out_path) / "blastout"
        blastout_folder.mkdir(parents=True, exist_ok=True)
    else :
        print("Skipping BLASTp.")
        
    # Perform the blastp function on each fasta
    fastas = [fasta for fasta in mapping_df['fasta']]
    with multiprocessing.Pool(processes=ncpus) as pool:
        all_hits = pool.map(partial(blastp, focal_fasta=focal_fasta, names_df=mapping_df, out_path=str(out_path), min_cov=min_cov, evalue=evalue, is_blast=is_blast), fastas)
        
    # all_hits = [blastp(fasta, focal_fasta, names_df, out_path, ncpus=ncpus, min_cov=min_cov, evalue=evalue, is_blast=is_blast) for fasta in names_df['fasta']]

    if not keep_files:
        shutil.rmtree(blastout_folder)
        shutil.rmtree(database_folder)

    # Build a dataframe from all of blastp hits
    hits_df = parse_blast_hits(hits=all_hits)

    # write outputs
    write_outputs(hits=hits_df, distances=distance_to_focal, out_basename=str(out_path / Path(focal_fasta).stem))


def run_orfdate_containerized(args: arguments.argparse.Namespace):
    # get the paths in the csv file and check if they exist
    mapping_df = get_mapping_df(csv_file=args.mapping, colnames=["fasta", "taxon"])
    check_files(files=mapping_df["fasta"].to_list())

    # remap paths on the container filesystem  
    extra_bindings = {x:"/input" for x in mapping_df["fasta"].to_list()}

    # add special arguments to inform orfdate that csv file must be modified to match the expected mapped paths on the container
    args.is_container = True
    cmd_args = sys.argv + ["--is-container"]  # required since ContainerCLI must rewrite the original argv

    # instantiate containerCLI handler
    cli = ContainerCLI(
            input_args=["--mapping", "--tree"],
            output_arg="--out",
            args=args,
            cmd_args=cmd_args,
            image_base=DOCKER_IMAGE,
            prog="orfdate",
            container_type="docker" if args.docker else "singularity",
            extra_bindings=extra_bindings,
            dev_mode=args.dev,
            package_binding={"orfmine": "/home/orfuser/orfmine/orfmine"}
        )

    cli.show()
    if not args.dry_run:
        cli.run()


def main():
    start_time = datetime.now()

    args = arguments.get_args()

    if args.docker or args.singularity:
        run_orfdate_containerized(args=args)
    else:
        run_orfdate(
            target=args.target,
            newick_file=args.tree,
            out_path=Path(args.out),
            mapping_file=args.mapping,
            ncpus=args.cpus,
            is_blast=args.blast,
            evalue=args.evalue,
            min_cov=args.min_coverage,
            keep_files=args.keep_files,
            preserve_underscores=args.has_underscores,
            is_container=args.is_container
        )

    end_time = datetime.now()
    print('\nDuration: {}'.format(end_time - start_time))
