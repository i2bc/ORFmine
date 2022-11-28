import argparse
import os
from pathlib import Path
import re
import sys

import psutil

from packages.orftrack.lib import gff_parser
from packages.orftrack.lib import fasta_parser


def memory_usage():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)

    return mem


def get_args():
    """
    Returns:
        Parameters
    """

    parser = argparse.ArgumentParser(description='Returns amino acid and nucleic fasta files from genomic data')
    group = parser.add_mutually_exclusive_group()

    parser.add_argument("-fna", required=True, nargs="?", help="Genomic fasta file (.fna) ")
    parser.add_argument("-gff", required=True, nargs="?", help="GFF annotation file (.gff)")
    parser.add_argument("-chrs", required=False, nargs="*", type=str, default=[], help="Chromosome names to processed. By default, all chromosomes are processed (default: []).")
    parser.add_argument("-outdir", required=False, nargs="?", default='./', type=str, help="Output directory")
    parser.add_argument("-outname", required=False, nargs="?", default="", type=str, help="Basename for output files (default: '_proteins' suffix added to fasta file name)")
    
    parser.add_argument(
            "-features_include",
            type=str,
            required=False, 
            nargs="*",
            default=[],
            help="Annotation features to be considered (By definition is all)"
    )

    group.add_argument(
        "-P", "--proteic",
        required=False,
        action='store_true',
        default=False,
        help="For only amino acids output format (.chromosome_to_process)."
    )

    group.add_argument(
        "-N", "--nucleic",
        required=False,
        action='store_true',
        default=False,
        help="For nucleic output format (.nfasta)."
    )

    parser.add_argument(
            "-elongate",
            type=int,
            required=False, 
            nargs="?",
            default=None,
            help="Number of nucleotides the sequence of matching GFF elements must be elongated from both extremities (default: None)"
    )

    args = parser.parse_args()

    return args


def main():
    args = get_args()
    genomic_fna = args.fna
    genomic_gff = args.gff
    chromosomes_to_process = args.chrs
    features = args.features_include
    elongation = args.elongate

    out_formats = [".pfasta", ".nfasta"]
    if args.proteic:
        out_formats.remove(".nfasta")
    elif args.nucleic:
        out_formats.remove(".pfasta")

    outpath = Path(args.outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    features_in_name = "_" + "_".join(features)
    if args.outname == "":
        basename_out = str(outpath / f"{Path(genomic_gff).stem}{features_in_name}")
    else:
        basename_out = str(outpath / f"{args.outname}{features_in_name}")

    print(f"Chromosomes to be processed: {', '.join(chromosomes_to_process)}") 
    print(f"Asked features: {', '.join(features)}\n")

    print("Parsing fasta file...")
    fasta_hash = fasta_parser.parse(fasta_filename=genomic_fna)
    
    print("Processing GFF file...")
    chromosomes = process_gff(
        gff_filename=genomic_gff,
        fasta_hash=fasta_hash,
        asked_chrs=chromosomes_to_process,
        features=features,
        basename_out=basename_out,
        out_formats=out_formats,
        elongation=elongation
    )

    if chromosomes:
        print("\nWriting CDS sequences...")
        with open(basename_out + '.pfasta', 'w') as faa_file:
            with open(basename_out + '.nfasta', 'w') as fna_file:
                for chr_name in sorted(chromosomes):
                    proteins = chromosomes[chr_name].group_cds()

                    for name, cds_elements in proteins.items():
                        header = '>' + name + '\n'
                        fasta_aa = header + ''.join([cds.translate() for cds in cds_elements]) + '\n'
                        fasta_nuc = header + ''.join([cds.sequence() for cds in cds_elements]) + '\n'
                        faa_file.write(fasta_aa)
                        fna_file.write(fasta_nuc)



def process_gff(gff_filename: str = None, fasta_hash: dict = {}, asked_chrs: list = [], features: list=[], basename_out: str = "", out_formats: list=[".pfasta", ".nfasta"], elongation: int = None) -> dict:
    """Return a dictionary with chromosome name as key, and gff_parser.Chromosome instance as value
    if asked GFF elements that are CDS only
    and / or 
    write on the fly fasta / gff file for asked GFF elements that are not CDS.

    Args:
        gff_filename (str, optional): GFF filename. Defaults to None.
        fasta_hash (str, optional): Fasta hash infos. Defaults to None.
        asked_chrs (list, optional): Chromosome names to processed. By default, all chromosomes are processed (default: []).
        features (list, optional): List of GFF features to process. Empty list means all features will be treated. Defaults to [].
        basename_out (str, optional): Output file basename. Defaults to "".
        out_formats (list, optional): List of output formats to be generated. Outputs can be .fna, .pfasta or both.

    Returns:
        dict: Dictionary with chromosome name as key, and gff_parser.Chromosome instance as value
              with CDS GFF elements only.  
    """
    chromosomes = {}
    features_to_keep = "|".join(features)

    if elongation:
        basename_out += "_elongated"
        out_formats = [".nfasta", ".gff"]
    
    try:
        out_files = { _ext:open(basename_out+_ext, "w") for _ext in out_formats }

        with open(gff_filename, 'r') as gff_file:
            for line in gff_file:
                if not line.startswith('#'):
                    chr_name = line.split('\t')[0]
                    element_type = line.split('\t')[2]

                    # process only asked chromosome name
                    if asked_chrs and chr_name not in asked_chrs:
                        print(f" - chromosome {chr_name} is not {asked_chrs}. Continue...", end='\r')
                        continue

                    # process only gff elements in features_to_keep:
                    if not re.findall(features_to_keep, element_type):
                        continue

                    gff_element = gff_parser.GffElement(gff_line=line, fasta_chr=fasta_hash[chr_name])

                    if elongation:
                        gff_element.start = gff_element.start - elongation
                        gff_element.end = gff_element.end + elongation

                    if not elongation and gff_element.type == "CDS":           
                        # create gff_parser.Chromosome instance 
                        if chr_name not in chromosomes:
                            print(f" - adding chromosome {chr_name}", end='\r')
                            chromosome = gff_parser.Chromosome(id_=chr_name, fasta_chr=fasta_hash[chr_name])
                            chromosomes[chr_name] = chromosome
                            chromosome.source = line.split('\t')[1]
                        chromosome.add(gff_element=gff_element)

                    else:
                        if ".pfasta" in out_files:
                            out_files[".pfasta"].write(gff_element.get_fastaline())
                        if ".nfasta" in out_files:
                            out_files[".nfasta"].write(gff_element.get_fastanuc_line())
                        if ".gff" in out_files:
                            out_files[".gff"].write(gff_element.get_gffline())

                    # print(f"## GFF data memory usage: {memory_usage()}", end="\r") 147594  151006, 151097  151166

    finally:
        for _, _file in out_files.items():
            _file.close()
            
    return chromosomes


if __name__ == '__main__':
    sys.exit(main())

    inpath = Path("/home/nchenche/projects/orfmine_workdir/examples/database")
    genomic_fna = inpath / "Scer.fna"
    genomic_gff = inpath / "Scer.gff"
    chromosome = ["I"]
    features = []

    fasta_hash = fasta_parser.parse(fasta_filename=genomic_fna)
    chromosomes = process_gff(
        gff_filename=genomic_gff,
        fasta_hash=fasta_hash,
        asked_chrs=chromosome,
        features=features,
    )

    chr = chromosomes[chromosome]

    proteins = chr.group_cds()

    bugprot = proteins["CDS:YAL001C"]
    cds1, cds2 = bugprot