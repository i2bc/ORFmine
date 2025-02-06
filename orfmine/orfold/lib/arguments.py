import argparse
from orfmine.utilities.container import add_container_args


def get_args():
    """
    Returns:
        Parameters
    """
    parser = argparse.ArgumentParser(description='ORF Foldability Calculation', allow_abbrev=False)

    parser.add_argument("--faa", "-F", type=str, required=True, help="FASTA file containing the amino acid sequences to treat")
    parser.add_argument("--gff", "-G", type=str, required=False, help="GFF annotation file")
    parser.add_argument("--options", "-P", type=str, required=True, help="Which properties are to be calculated. H for HCA (Default), I for IUPred, T for Tango (e.g. `-options HI` or `-options -HT` or `-options HIT` ...)")
    parser.add_argument("--out", "-O", type=str, required=False, default='.', help="Output directory ('.' by default).")
    parser.add_argument("--barcodes", required=False, action='store_true', default=False, help=argparse.SUPPRESS)
    parser.add_argument("--keep", "-K", required=False, action='store_true', default=False, help="Option for keeping the Tango output files")
    parser.add_argument("--sample", "-N", type=int, required=False, default=-1, help="Size of the sample to use for the fasta sequences. Defaults to -1 (i.e. all sequences).")
    parser.add_argument("--path_tango", type=str, required=False, help="Path to the Tango executable")
    parser.add_argument("--path_iupred", type=str, required=False, help="Path to the IUPred directory")



    parser = add_container_args(parser=parser)

    args = parser.parse_args()

    if args.options is None or args.options not in "HIT":
        parser.error("-options requires at least one argument amongst: H, I, T (e.g. `-options HI` or `-options -HT` or `-options HIT` ...)")
    #if "T" in args.options and not args.path_tango:
    #    parser.error("--path_tango is required when 'T' option is specified in --options")
    #if "I" in args.options and not args.path_iupred:
    #    parser.error("--path_iupred is required when 'I' option is specified in --options")
    
    return args


