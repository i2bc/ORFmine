import argparse
import os
from pathlib import Path
import re
import sys
from threading import Thread
from typing import List
import time

import psutil

from packages.orftrack.lib import gff_parser
from packages.orftrack.lib import fasta_parser

import multiprocessing


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
            required=True, 
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

    parser.add_argument(
        "-chr_exclude",
        type=str,
        required=False, 
        nargs="*",
        default = [],
        help="Chromosomes to be excluded (By definition is None)"
    )

    args = parser.parse_args()

    return args


GFF_DESCR = {}


def memory_usage():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)

    return mem


def set_gff_descr(gff_fname):
    global GFF_DESCR
    GFF_DESCR = {}
        
    with open(gff_fname, 'rb') as gff_file:
        line = gff_file.readline().decode(encoding='utf-8')
        while line:
            if not line.startswith('#'):
                name = line.split('\t')[0]
                pos_chr = gff_file.tell()-len(line)
                if name not in GFF_DESCR:
                    GFF_DESCR[name] = pos_chr
                
            line = gff_file.readline().decode(encoding='utf-8')


class CDSQueue:
    def __init__(self) -> None:
        self.cds_list: List[gff_parser.GffElement] = []
        self.stored_strand: str = ""
        self.stored_parent: str = ""
        self.completed: bool = False
        self.cds_completed: List[gff_parser.GffElement] = []

    def add(self, cds: gff_parser.GffElement=None):
        if cds.parent != self.stored_parent:
            self.update(cds)
        else:
            self.cds_completed.clear()

        if cds.strand == "-":
            if not self.cds_list:
                self.cds_list.append(cds)
            else:
                if cds.start < self.cds_list[0].start:
                    self.cds_list.append(cds)
                else:
                    self.cds_list.insert(0, cds)
        else:
            self.cds_list.append(cds)


    def update(self, cds: gff_parser.GffElement=None):
        self.stored_parent = cds.parent

        self.cds_completed = [x for x in self.cds_list]
        self.cds_list.clear()
    
    def get_fasta(self, _type="nucleic"):
        try:
            header = '>' + self.cds_completed[0].name + '\n'
        except:
            print(f"Exception, list error: {self.cds_completed}")
        if _type == "nucleic":
            return header + ''.join([cds.sequence() for cds in self.cds_completed]) + '\n'
        else:
            return header + ''.join([cds.translate() for cds in self.cds_completed]) + '\n'   

    def elongate(self, value):
        # elongation of CDS in borders only
        if self.cds_completed[0].strand == "-":
            cds_first = self.cds_completed[-1]
            cds_last = self.cds_completed[0]
        else:
            cds_first = self.cds_completed[0]
            cds_last = self.cds_completed[-1]

        cds_first.start = cds_first.start - value
        cds_last.end = cds_last.end + value


def process_chromosome(gff_filename: str=None, fasta_chr: fasta_parser.Fasta=None, chr_id: str="", features: list=[], basename_out: str="", out_formats: list=[".pfasta", ".nfasta"], elongation: int=None) -> dict:
    """
    Return a dictionary with chromosome name as key, and gff_parser.Chromosome instance as value
    if asked GFF elements that are CDS only
    and / or 
    write on the fly fasta / gff file for asked GFF elements that are not CDS.

    Args:
        gff_filename (str, optional): GFF filename. Defaults to None.
        fasta_chr (str, optional): fasta_parser.Fasta instance. Defaults to None.
        asked_chrs (list, optional): Chromosome names to be processed. By default, all chromosomes are processed (default: []).
        features (list, optional): List of GFF features to process. Empty list means all features will be treated. Defaults to [].
        basename_out (str, optional): Output file basename. Defaults to "".
        out_formats (list, optional): List of output formats to be generated. Outputs can be .fna, .pfasta or both.
        asked_chrs (list, optional): Chromosome names not to be processed. By default, all chromosomes are processed (default: []).

    Returns:
        dict: Dictionary with chromosome name as key, and gff_parser.Chromosome instance as value
              with CDS GFF elements only.  
    """
    start_time = time.time()

    process_name = multiprocessing.current_process().name
    print(f"{process_name} -  Starting to process chromosome {chr_id}")

    features_to_keep = "|".join(features)

    # queue used to process CDS particular (CDS can't be written 'on the fly')
    queue = CDSQueue()

    if elongation:
        basename_out += "_elongated"
        out_formats = [".nfasta", ".gff"]
    
    try:
        out_files = { _ext:open(basename_out+_ext, "w") for _ext in out_formats }

        with open(gff_filename, 'r') as gff_file:
            # retrieve en-of-file value
            eof = gff_file.seek(0, 2)

            # place cursor to the position of the 1st given chromosome line 
            gff_file.seek(GFF_DESCR[chr_id], 0)

            # read line as long as it is not starting with '#' and chr name = chr_id
            while True:
                line = gff_file.readline()
                # skip lines starting with a comment mark
                if line.startswith("#"):
                    continue

                # quit the loop in chromosome name is not the chromosome of interest
                chr_name = line.split('\t')[0]
                if chr_name != chr_id:
                    break

                element_type = line.split('\t')[2]
                
                # skip all gff elements not matching features_to_keep:
                if not re.findall(features_to_keep, element_type):
                    continue

                try:
                    gff_element = gff_parser.GffElement(gff_line=line, fasta_chr=fasta_chr)
                except KeyError as e:
                    print(f"Warning: {chr_name} seems to be absent from the fasta file. Continue...")
                

                # process non CDS elements
                if gff_element.type != "CDS":
                    if elongation:
                        gff_element.start = gff_element.start - elongation
                        gff_element.end = gff_element.end + elongation
                    
                    # writing non CDS elements
                    if ".pfasta" in out_files:
                        out_files[".pfasta"].write(gff_element.get_fastaline())
                    if ".nfasta" in out_files:
                        out_files[".nfasta"].write(gff_element.get_fastanuc_line())
                    if ".gff" in out_files:
                        out_files[".gff"].write(gff_element.get_gffline())

                # CDS elements must be treated differently (CDS must be merged for multi-exonic proteins)
                else:
                    queue.add(cds=gff_element)
                    if queue.cds_completed:
                        # print(", ".join([x.name for x in queue.cds_completed]))
                        if elongation:
                            queue.elongate(value=elongation)

                        # writing non CDS elements
                        if ".pfasta" in out_files:
                            out_files[".pfasta"].write(queue.get_fasta(_type="proteic"))
                        if ".nfasta" in out_files:
                            out_files[".nfasta"].write(queue.get_fasta(_type="nucleic"))
                        if ".gff" in out_files:
                            out_files[".gff"].write(gff_element.get_gffline())                        

                if gff_file.tell() == eof:
                    break

        if queue.cds_list:
            queue.cds_completed = queue.cds_list
            if elongation:
                queue.elongate(value=elongation)

            # writing non CDS elements
            if ".pfasta" in out_files:
                out_files[".pfasta"].write(queue.get_fasta(_type="proteic"))
            if ".nfasta" in out_files:
                out_files[".nfasta"].write(queue.get_fasta(_type="nucleic"))
            if ".gff" in out_files:
                out_files[".gff"].write(gff_element.get_gffline())
  
    finally:
        for _, _file in out_files.items():
            _file.close()

    print(f"{process_name} -  Chromosome {chr_id} successfully processed in {round(time.time()-start_time, 2)} seconds")
                                

def main():
    args = get_args()
    genomic_fna = args.fna
    genomic_gff = args.gff
    chromosomes_to_process = args.chrs
    chromosomes_to_exclude = args.chr_exclude
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

    print("Parsing gff file...")
    if not GFF_DESCR:
        set_gff_descr(genomic_gff)

    processes = []
    for i, chr_id in enumerate(fasta_hash):
        if chromosomes_to_process and chr_id not in chromosomes_to_process:
            continue
        elif chromosomes_to_exclude and chr_id in chromosomes_to_exclude:
            continue

        basename = f"{basename_out}_{chr_id}_{i}"

        p = multiprocessing.Process(target=process_chromosome, args=(genomic_gff, fasta_hash[chr_id], chr_id, features, basename, out_formats, elongation,))
        p.chr_id = chr_id
        processes.append(p)

    for p in processes:
        p.start()        


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

    Genecode = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
            }

    std_table = table.standard_dna_table

    for k, v in std_table.forward_table.items():
        if Genecode[k] == v:
            print(k, v, Genecode[k])