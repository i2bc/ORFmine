#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 14:40:40 2022

@author: paul.roginski
"""

# Standard library imports
from datetime import datetime
from functools import partial
import multiprocessing
from pathlib import Path
import operator
import shutil
import sys
from typing import Dict, List, Union
import subprocess  # Ajout pour remplacer Bio.Application

# Third-party imports
import dendropy
import numpy as np
import pandas as pd
from orfmine import DOCKER_IMAGE

# Self-party imports
from orfmine.orfdate.lib import arguments
from orfmine.utilities.container import ContainerCLI


def generate_tree(tree_file: str, names_df: pd.DataFrame, preserve_underscores: bool) -> dendropy.Tree:
    tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=preserve_underscores)
    print("Original tree")
    print(tree.as_ascii_plot())

    # Remove from tree taxa absent in names_df
    taxa = [tree.taxon_namespace[i].label for i in range(len(tree.taxon_namespace))]
    extra_taxa = [label for label in taxa if label not in names_df['taxon'].tolist()]

    if extra_taxa:
        print(f"Extra taxa are: {extra_taxa}")
        tree.prune_taxa_with_labels(extra_taxa)
        tree.taxon_namespace = [tree.taxon_namespace[i] for i in range(len(tree.taxon_namespace))
                                if tree.taxon_namespace[i].label in names_df['taxon'].tolist()]
        print("Corrected tree")
        print(tree.as_ascii_plot())

    return tree


def map_pyhlo_distance(tree: dendropy.Tree, focal_name: str) -> Dict:
    try:
        focal = tree.taxon_namespace.get_taxon(label=focal_name)
    except:
        labels = [tree.taxon_namespace[i].label for i in range(len(tree.taxon_namespace))]
        print(f"Cannot find the target taxon {focal_name} in the names of the tree:\n{labels}")
        sys.exit(1)

    pdc = tree.phylogenetic_distance_matrix()
    distance_to_focal = {t.label: pdc(focal, t) for t in tree.taxon_namespace}

    if sum(distance_to_focal.values()) == 0:
        print("No distance is specified in the newick file. Exiting.")
        sys.exit(1)

    return dict(sorted(distance_to_focal.items(), key=operator.itemgetter(1)))


def blastp(fasta, focal_fasta, names_df, out_path, num_threads=4, min_cov=0.7, evalue=0.001, is_blast=True):
    database_folder = Path(out_path) / "blastdb"
    blastout_folder = Path(out_path) / "blastout"

    # Taxon name
    tree_name = names_df['taxon'][np.where(names_df['fasta'] == fasta)[0][0]]
    print(f"{tree_name}...")

    # BLASTp output file's name
    output_file = blastout_folder / f"{Path(focal_fasta).stem}_vs_{Path(fasta).stem}.out"

    if is_blast:
        db = database_folder / Path(fasta).stem

        # Build a protein BLAST database
        subprocess.run([
            "makeblastdb",
            "-dbtype", "prot",
            "-in", str(fasta),
            "-out", str(db)
        ], check=True)

        # Perform BLASTp
        subprocess.run([
            "blastp",
            "-query", str(focal_fasta),
            "-db", str(db),
            "-out", str(output_file),
            "-outfmt", "6 qseqid sseqid evalue qlen qstart qend slen sstart send length bitscore score",
            "-evalue", str(evalue),
            "-num_threads", str(num_threads)
        ], check=True)

    # Parse BLASTp output
    hits = {}
    with open(output_file, 'r') as f:
        for line in f:
            if line[0] != "#":
                full_line = line.strip().split("\t")
                query_cov = (int(full_line[5]) - int(full_line[4])) / int(full_line[3])
                if float(full_line[2]) < evalue and query_cov >= min_cov:
                    query = full_line[0]
                    hits[query] = hits.get(query, 0) + 1

    return [tree_name, hits]


def parse_blast_hits(hits: List) -> pd.DataFrame:
    print("Merging all hits in one dataframe...")
    hits_df = pd.DataFrame([ele[1] for ele in hits]).transpose()
    hits_df.columns = [ele[0] for ele in hits]
    hits_df = hits_df.fillna(0).astype(int)
    return hits_df


def write_outputs(hits: pd.DataFrame, distances: Dict, out_basename: str) -> None:
    hits = hits[list(distances.keys())]
    hits.to_csv(f"{out_basename}_hits.csv", index=True, index_label='seq')

    farest_df = pd.DataFrame(index=hits.index, columns=["farest_hit", "distance"])
    for orf in hits.index:
        row = hits.loc[orf, :]
        hit_names = [name for name in hits.columns if row[name] > 0]
        hit_distances = [distances[name] for name in hit_names]
        max_distance = max(hit_distances, default=0)
        farest_names = [name for name in hit_names if distances[name] == max_distance]
        farest_df.loc[orf, "farest_hit"] = "|".join(farest_names)
        farest_df.loc[orf, "distance"] = max_distance

    farest_df.to_csv(f"{out_basename}_dated.csv", index=True, index_label='seq')


def remap_to_dest(df: pd.DataFrame, colname: str, dest: str):
    df[colname] = df[colname].apply(lambda x: str(Path(dest) / Path(x).name))
    return df


def get_mapping_df(csv_file: Union[str, Path], colnames: List[str]) -> pd.DataFrame:
    try:
        mapping_df = pd.read_csv(csv_file, names=colnames)
        if mapping_df.empty:
            raise ValueError("The mapping file is empty.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    return mapping_df


def check_files(files: Union[str, Path]):
    for _file in files:
        if not Path(_file).exists():
            print(f"Error: unable to find {_file}")
            sys.exit(1)


def run_orfdate(target: str, newick_file: Union[str, Path], out_path: Union[str, Path], mapping_file: Union[str, Path],
                ncpus: int = 1, is_blast: bool = True, evalue: float = 0.001, min_cov: float = 0.7,
                keep_files: bool = False, preserve_underscores: bool = False, is_container: bool = False):
    mapping_df = get_mapping_df(csv_file=mapping_file, colnames=["fasta", "taxon"])
    check_files(mapping_df["fasta"].tolist())

    focal_fasta = mapping_df['fasta'][np.where(mapping_df['taxon'] == target)[0][0]]
    tree = generate_tree(tree_file=newick_file, names_df=mapping_df, preserve_underscores=preserve_underscores)
    distance_to_focal = map_pyhlo_distance(tree=tree, focal_name=target)

    if is_blast:
        database_folder = Path(out_path) / "blastdb"
        database_folder.mkdir(parents=True, exist_ok=True)
        blastout_folder = Path(out_path) / "blastout"
        blastout_folder.mkdir(parents=True, exist_ok=True)

    fastas = mapping_df["fasta"].tolist()
    with multiprocessing.Pool(processes=ncpus) as pool:
        all_hits = pool.map(partial(blastp, focal_fasta=focal_fasta, names_df=mapping_df, out_path=out_path,
                                    min_cov=min_cov, evalue=evalue, is_blast=is_blast), fastas)

    if not keep_files:
        shutil.rmtree(database_folder)
        shutil.rmtree(blastout_folder)

    hits_df = parse_blast_hits(hits=all_hits)
    write_outputs(hits=hits_df, distances=distance_to_focal, out_basename=str(Path(out_path) / Path(focal_fasta).stem))


def main():
    args = arguments.get_args()

    if args.docker or args.singularity:
        print("Running in containerized mode.")
        run_orfdate_containerized(args=args)
    else:
        run_orfdate(
            target=args.target,
            newick_file=args.tree,
            out_path=args.out,
            mapping_file=args.mapping,
            ncpus=args.cpus,
            is_blast=args.blast,
            evalue=args.evalue,
            min_cov=args.min_coverage,
            keep_files=args.keep_files,
            preserve_underscores=args.has_underscores
        )


if __name__ == "__main__":
    main()
