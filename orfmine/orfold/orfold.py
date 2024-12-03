#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ORFold - ORF Foldability Calculation

@author: christospapadopoulos
"""

import sys
from datetime import datetime
from pathlib import Path
from functools import partial
import multiprocessing
from typing import Union, Dict, List
import argparse

from orfmine.utilities.container import add_container_args
from orfmine.orfold.lib import utils as orfold_utils


def get_args():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description='ORF Foldability Calculation', allow_abbrev=False)

    parser.add_argument("--faa", "-F", type=str, required=True, help="FASTA file containing the amino acid sequences to treat")
    parser.add_argument("--gff", "-G", type=str, required=False, help="GFF annotation file")
    parser.add_argument("--options", "-P", type=str, required=False, default="H", 
                        help="Which properties are to be calculated. H for HCA (default), I for IUPred, T for Tango")
    parser.add_argument("--out", "-O", type=str, required=False, default='.', help="Output directory (default: '.').")
    parser.add_argument("--keep", "-K", required=False, action='store_true', default=False, help="Option for keeping the Tango output files")
    parser.add_argument("--sample", "-N", type=int, required=False, default=-1, help="Size of the sample to use for the fasta sequences. Defaults to -1 (all sequences).")

    parser = add_container_args(parser=parser)

    args = parser.parse_args()

    # Validate options for --options
    if args.options is None or not set(args.options).issubset("HIT"):
        parser.error("--options requires at least one argument among: H, I, T (e.g., --options HI or --options HIT)")

    return args


def process_fasta_file(fasta_file: Union[str, Path], out_path: Union[str, Path], options: str="H", sample_size: int=-1, to_keep: bool=False):
    """
    Processes the provided FASTA file based on the selected options.

    Args:
        fasta_file (Union[str, Path]): Input FASTA file.
        out_path (Union[str, Path]): Output directory.
        options (str): Selected options for analysis ("H", "I", "T").
        sample_size (int): Number of sequences to process. Defaults to all (-1).
        to_keep (bool): Whether to keep intermediate files.

    Returns:
        Dict: Scores for each sequence.
    """
    print(f"Processing FASTA file: {fasta_file}")
    out_path = Path(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    # Simulated processing logic (replace with actual logic as needed)
    scores = {}
    fasta_sequences = orfold_utils.fasta_generator(fasta_file)
    for header, sequence in fasta_sequences:
        scores[header] = {opt: "SimulatedValue" for opt in options}
    return scores


def run_orfold(fasta_file: Union[str, Path], out_path: Union[str, Path], options: str="H", sample_size: int=-1, to_keep: bool=False):
    """
    Main logic for processing ORF foldability.

    Args:
        fasta_file (Union[str, Path]): Input FASTA file.
        out_path (Union[str, Path]): Output directory.
        options (str): Selected options for analysis ("H", "I", "T").
        sample_size (int): Number of sequences to process. Defaults to all (-1).
        to_keep (bool): Whether to keep intermediate files.
    """
    scores = process_fasta_file(fasta_file, out_path, options, sample_size, to_keep)
    print(f"Processed scores: {scores}")


def main():
    """
    Main entry point for the script.
    """
    # If no arguments are provided, display help
    if len(sys.argv) == 1:
        get_args().print_help()
        sys.exit(0)

    # Parse arguments
    args = get_args()

    # Run the main logic
    start_time = datetime.now()
    run_orfold(
        fasta_file=args.faa,
        out_path=args.out,
        options=args.options,
        sample_size=args.sample,
        to_keep=args.keep
    )
    end_time = datetime.now()
    print(f"Duration: {end_time - start_time}")


if __name__ == "__main__":
    main()

