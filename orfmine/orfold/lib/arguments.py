import argparse


def get_args():
    """
    Returns:
        Parameters
    """
    parser = argparse.ArgumentParser(description='ORF Foldability Calculation', allow_abbrev=False)

    parser.add_argument("-faa", type=str, action='store', required=True, nargs="*", help="FASTA file containing the amino acid sequences to treat")
    parser.add_argument("-gff", required=False, type=str, nargs="*", default=[], help="GFF annotation file")
    parser.add_argument("-options", type=list, required=True, nargs="?", help="Which properties are to be calculated. H for HCA (Default), I for IUPred, T for Tango")
    parser.add_argument("-out", required=False, nargs="?", default='.', type=str, help="Output directory ('.' by default).")
    parser.add_argument("--barcodes", required=False, action='store_true', default=False, help=argparse.SUPPRESS)
    parser.add_argument("--keep", required=False, action='store_true', default=False, help="Option for keeping the Tango output files")
    parser.add_argument("-N", required=False, type=str, nargs="*", default=["all"], help="Size of sample(s) per FASTA file")
    parser.add_argument("-D", "--dry-run", required=False, action='store_true', default=False, help="Flag used to show the docker command line. Must be used in conjonction with '--docker' or '--singularity'")

    container_group = parser.add_mutually_exclusive_group()
    container_group.add_argument("--docker", action='store_true', default=False, help="Flag used to run computations on a docker container")
    container_group.add_argument("--singularity", action='store_true', default=False, help="Flag used to run computations on a singularity container")

    args = parser.parse_args()

    if args.options is None:
        parser.error("-options requires at least one argument amongst: H, I, T (e.g. `-options HI` or `-options -HT` or `-options HIT` ...)")

    
    return args


