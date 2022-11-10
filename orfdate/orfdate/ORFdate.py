#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 14:40:40 2022

@author: paul.roginski
"""
# TODO : determine the number of cpu to use per blast
# TODO : gérer les fichiers de correspondance incomplets

import argparse
import pandas as pd
import numpy as np
# import multiprocessing
import dendropy
import operator
import os
import sys

def run():
    # Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-focal", required=True, help="taxonomic name of the focal species")
    parser.add_argument("-names", required=True, help="csv file matching fasta (col1) and tree (col2) names", default=False)
    parser.add_argument("-tree", required=True, help="newick file for the phylogeny tree")
    parser.add_argument("-ncpus", help="total number of cpus that can be used fo the task", default = 1, type=int)
    parser.add_argument("-blast", help="(logical) wether or not to perform BLASTp", default = True)
    parser.add_argument("-evalue", help="BLASTp evalue threshold", default = 0.001, type=float)
    parser.add_argument("-query_cov", help="minimum query coverage threshold", default = 0.7, type=int)
    parser.add_argument("-preserve_underscores", help="should underscores be kept when reading the tree, or considered as spaces", default = False, type=bool)
    parser.add_argument("-db", help="directory for blastdb files", default = "/workdir/orfdate/blastdb", type=str)
    parser.add_argument("-blastout", help="directory for blast output files", default = "/workdir/orfdate/blastout", type=str)

    args = parser.parse_args()


    # taxonomic name of the focal species
    focal_name = args.focal

    # Dataframe with fasta files' names as first column and tree taxons' names as second column
    names_df = pd.read_csv(args.names, names = ["fasta","tree"])

    # fasta name of the focal species
    try :
        focal_fasta = names_df['fasta'][np.where(names_df['tree'] == focal_name )[0][0]]
    except :
        print("The focal name {} is not found.\nNames provided in -names :\n{}".format(focal_name,names_df['tree']))
        sys.exit(1)

    # newick file for the phylogeny tree
    tree_file = args.tree

    # Number of cpus to use
    ncpus = args.ncpus

    # blast (logical) wether or not to perform BLAST
    ua = str(args.blast).upper()
    if 'TRUE'.startswith(ua):
       blast = True
    elif 'FALSE'.startswith(ua):
       blast = False
    else:
       print("Error : -blast argument must be True of False.")

    # BLASTp evalue threshold
    evalue = args.evalue

    # minimum query coverage threshold
    min_cov = args.query_cov


    def path_plus_base(path) :
        return(os.path.splitext(path)[0])
    def basename(path):
        return(os.path.basename(path).split(".")[0])





    # Generate and print a tree object that will be used to infer hits distance.
    tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores= args.preserve_underscores)
    print("Original tree")
    print(tree.as_ascii_plot())


    # Identify among the taxa of the tree, those that are absent in name_df and
    # remove them from the tree.
    taxa = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]
    extra_taxa = [ label for label in taxa if label not in names_df['tree'].tolist()]
    if len(extra_taxa) > 0 :
        print("Extra taxa are : {}".format(extra_taxa))
        tree.prune_taxa_with_labels(extra_taxa)
        print("Corrected tree")
        print(tree.as_ascii_plot())
        # Also need to correct the "taxon_namespace" attribute
        tree.taxon_namespace = [ tree.taxon_namespace[i] for i in range(0,len(tree.taxon_namespace)) if tree.taxon_namespace[i].label in names_df['tree'].tolist()]


    # Retrieve in the taxon_namespace attribute of the tree object, the element
    # corresponding with the focal species.
    try :
        focal = tree.taxon_namespace.get_taxon(label = focal_name)
    except :
        labels = [tree.taxon_namespace[i].label for i in range(len(tree.taxon_namespace))]
        print("Cannot find the -focal {} in the names of the tree :\n{}".format(focal_name,labels))
        sys.exit(1)


    # Build a distance matrix object
    pdc = tree.phylogenetic_distance_matrix()


    # Build a dictionnary with species as keys and distance to the focal species as values
    distance_to_focal = { t.label:pdc(focal, t) for t in tree.taxon_namespace}
    # Sorted version
    distance_to_focal = dict(sorted(distance_to_focal.items(), key=operator.itemgetter(1)))
    # if no distance is specified :
    if sum(distance_to_focal.values()) == 0 :
        print("No distance is specified in the newick file. Exitting.")
        sys.exit(1)


    if blast :
        print("Performing BLASTp on provided fastas...")
        from Bio.Blast.Applications import NcbimakeblastdbCommandline
        from Bio.Blast.Applications import NcbiblastpCommandline

        database_folder = path_plus_base(focal_fasta)
        if args.db : database_folder = args.db
        if database_folder[-1] == "/" : database_folder = database_folder[0:-1]

        blastout_folder = path_plus_base(focal_fasta)
        if  args.blastout : blastout_folder = args.blastout
        if blastout_folder[-1] == "/" : blastout_folder = blastout_folder[0:-1]
        if not os.path.exists(blastout_folder):
            os.makedirs(blastout_folder)

    else :
        print("Skipping BLASTp.")


    def blastp(fasta):

            # Taxon name
            tree_name = names_df['tree'][np.where(names_df['fasta'] == fasta )[0][0]]
            print("{}...".format(tree_name))

            # BLASTp output file's name
            output_file = "{}/{}_vs_{}.out".format(blastout_folder,
                                                   path_plus_base(focal_fasta),
                                                   basename(fasta))

            if blast :

                db = database_folder+"/"+basename(fasta)

                # Build a protein BLAST database from the subject fasta
                cline = NcbimakeblastdbCommandline(dbtype="prot",
                                                   input_file = fasta,
                                                   out = db)
                cline()

                # Perfom a local BLASTp of the focal fasta against the current fasta's database
                cline = NcbiblastpCommandline(query=focal_fasta,
                                              db= db,
                                              out=output_file,
                                              outfmt = "6 qseqid sseqid evalue qlen qstart qend slen sstart send length bitscore score",
                                              evalue = 10,
                                              num_threads = args.ncpus)
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
            return([tree_name, hits])




    # Perform the blastp function on each fasta (parallelized on ncpus).
    # The list of returned objects is assigned to "all_hits"
    # pool = multiprocessing.Pool(ncpus)
    # all_hits = pool.map(blastp, names_df['fasta'])
    all_hits = [blastp(fasta) for fasta in names_df['fasta']]

    # Build a dataframe from "all_hits". This dataframe, once transposed and
    # renamed (columns), will have focal ORFs as rows and tree's taxa as columns.
    print("Merging all hits...")
    hits_df = pd.DataFrame([hits[1] for hits in all_hits])
    hits_df = hits_df.transpose()
    hits_df.columns = [hits[0] for hits in all_hits]
    hits_df = hits_df.fillna(0)
    hits_df = hits_df.astype(int)

    # re-order the columns of hits_df according to their distance to the focal species
    hits_df = hits_df[list(distance_to_focal.keys())]

    # Write hits_df to a csv file for the user
    print("Writing the hits tableau...")
    hits_df.to_csv("{}_hits.csv".format("/workdir/orfdate/"+basename(focal_fasta)),index=True, index_label='seq')


    # Generate an empty dataframe that will be filled with :
    # first column = the list of the farsest species with a hit
    # second column = the distance of these species to the focal species.
    print("Building the output csv...")
    farest_df = pd.DataFrame(index=hits_df.index,
                              columns = ["farest_hit","distance"])

    # For each ORF
    for orf in hits_df.index :

        # Its corresponding row in hits_df
        row = hits_df.loc[orf,:]

        # Species-with-hits's names
        hit_names = [name for name in hits_df.columns if row[name] > 0 ]

        # Their distance to the focal species
        hit_distances = [ distance_to_focal[hit_name] for hit_name in hit_names ]

        # The maximum distance to the focal species
        max_distance = max(hit_distances)

        # Names of the species with hits that are at this maximum distance
        farest_names = [name for name in hit_names if distance_to_focal[name] == max_distance ]

        farest_df.loc[orf,"farest_hit"] = "|".join(farest_names)
        farest_df.loc[orf,"distance"] = max_distance / 2 # because patristic distance is the sum of the lengths of the branches that link two nodes in a tree

    # Write hits_df to a csv file for the user
    print("Writing the output csv...")
    farest_df.to_csv("{}_dated.csv".format("/workdir/orfdate/"+basename(focal_fasta)),index=True, index_label='seq')
