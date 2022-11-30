import argparse
import logging
import os
from pathlib import Path
import re
import sys
from threading import Thread
from typing import List
import time

import psutil

from packages.orftrack.lib import gff_parser, fasta_parser, inspect

import multiprocessing


logger = logging.getLogger()

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
        default=[],
        help="Chromosomes to be excluded (By definition is None)"
    )

    parser.add_argument(
        "-cpus",
        type=int,
        required=False, 
        nargs="?",
        default=multiprocessing.cpu_count(),
        help="Maximum number of CPUs to use. By default, all available CPUs are used."
    )

    args = parser.parse_args()

    return args


class IndexFasta(Thread):
    def __init__(self, filename):
        super().__init__()
        self.filename = filename
        self.result = None

    def run(self):
        print(f"{self.name} - Indexing fasta file...")
        start_time = time.time()
        self.result = fasta_parser.parse(fasta_filename=self.filename)
        print(f"{self.name} - Fasta file parsed in {round(time.time()-start_time, 2)} seconds")


class IndexGFF(Thread):
    def __init__(self, filename):
        super().__init__()
        self.filename = filename
        self.result = None

    def run(self):
        print(f"{self.name} - Indexing GFF file...")
        start_time = time.time()
        self.result = set_gff_descr(gff_fname=self.filename)
        print(f"{self.name} - GFF file parsed in {round(time.time()-start_time, 2)} seconds")



class CDSQueue:
    """
    Object used to process CDS, i.e. make them ready to be written, nearly "on-the-fly".


    Workflow attributes
    -------------------
    cds_list: List[gff_parser.GffElement] = []  # stores parent-related CDS unready
    stored_parent: str = ""  # stored last read parent
    cds_completed: List[gff_parser.GffElement] = []  # stores parent-related CDS ready to be merged
    

    Workflow explanation
    --------------------

    Every time a CDS is read, it is added to the CDSQueue instance.
    When a CDS is added, we first check if the precedently stored parent 
    is the same as the parent of the currently read/added CDS:
        - if both parents are the same, we add the current CDS to cds_list and clear cds_completed content
        - if both parents are different:
            1. we copy cds_list elements into the cds_completed list
            2. we update the stored parent value
            3. we clear the cds_list content
            4. we add the current CDS to cds_list

    To know when CDS is ready to be written (i.e fully merged in case of multi-exon proteins),
    we just have to check if queue.cds_completed is empty or not. If not empty, CDS can we written.

    
    Workflow illustration
    --------------------

    actual CDS      stored.parent   update()?   stored.parent   cds_completed status            cds_list status
    ----------      -------------   ---------   -------------   --------------------            ---------------
    cds1-parentA    None            YES         parentA         []                              [cds1-parentA]
    cds1-parentB    parentA         YES         parentB         [cds1-parentA]                  [cds1-parentB]
    cds2-parentB    parentB         NO          parentB         []                              [cds1-parentB, cds2-parentB]
    cds1-parentC    parentB         YES         parentC         [cds1-parentB, cds2-parentB]    [cds1-parentC]
    ...


    Workflow usage pseudo-code
    --------------------------
    queue = CDSQueue()

    for line in gff_file:
        element_type = line.split('\t')[2]
        gff_element = gff_parser.GffElement(gff_line=line, fasta_chr=fasta_chr)

        # we add to the queue every CDS element
        if gff_element.type == "CDS":
            queue.add(cds=gff_element)

        # we check if CDS are ready to be written (i.e. merged in case of multi-exonic proteins)
        if queue.cds_completed:
            # print(", ".join([x.name for x in queue.cds_completed]))
            write(queue.get_fasta())

    """
    def __init__(self) -> None:
        self.cds_list: List[gff_parser.GffElement] = []  # stores parent-related CDS unready
        self.stored_parent: str = ""  # stored last read parent
        self.cds_completed: List[gff_parser.GffElement] = []  # stores parent-related CDS ready to be merged

    def add(self, cds: gff_parser.GffElement=None):
        # when CDS parent has not been met yet, we update the state of CDSQueue instance
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


def memory_usage():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)

    return mem


def set_gff_descr(gff_fname):
    gff_indexes = {}
        
    with open(gff_fname, 'rb') as gff_file:
        line = gff_file.readline().decode(encoding='utf-8')
        while line:
            if not line.startswith('#'):
                name = line.split('\t')[0]
                pos_chr = gff_file.tell()-len(line)
                if name not in gff_indexes:
                    gff_indexes[name] = pos_chr
                
            line = gff_file.readline().decode(encoding='utf-8')

    return gff_indexes


def process_chromosome(gff_filename: str=None, fasta_chr: fasta_parser.Fasta=None, chr_id: str="", features: list=[], basename_out: str="", out_formats: list=[".pfasta", ".nfasta"], elongation: int=None, gff_indexes: dict={}) -> dict:
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
            gff_file.seek(gff_indexes[chr_id], 0)

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
        
        # The last CDS is still in cds_list
        if queue.cds_list:
            
            # must copy cds_list in cds_completed since queue.get_fasta() applies only on cds_completed
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
                                


def get_chunks(l: list, chunks_size: int=10):     
    for i in range(0, len(l), chunks_size):
        yield l[i:i+chunks_size]


def validate_chromosomes(to_include: list=[], to_exclude: list=[], in_gff: list=[], in_fasta: list=[]) -> list:
    # get chromosomes present both in the GFF file and the fasta file
    chrs_common = inspect.check_chrids(chrs_gff=sorted(in_gff), chrs_fasta=sorted(in_fasta))

    # get chromosomes to be treated
    if not to_include:
        chr_ids = sorted([x for x in chrs_common if x not in to_exclude])
    else:
        chr_ids = sorted([x for x in to_include if x not in to_exclude and x in chrs_common])

    for chromosome in to_include:
        if chromosome not in chr_ids:
            logger.warning(f"Warning: chromosome {chromosome} not found in either the fasta or GFF file. It will be skipped...")

    return chr_ids


def index_genomes(fasta_file: str="", gff_file: str=""):
    indexing_threads = [IndexFasta(filename=fasta_file), IndexGFF(filename=gff_file)]
    for t in indexing_threads:
        t.start()

    for t in indexing_threads:
        t.join()

    return [x.result for x in indexing_threads]


def main():
    args = get_args()
    genomic_fna = args.fna
    genomic_gff = args.gff
    chromosomes_to_process = args.chrs
    chromosomes_to_exclude = args.chr_exclude
    features = args.features_include
    elongation = args.elongate
    cpus = args.cpus

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

    # get useful indexes of fasta & gff files 
    fasta_hash, gff_indexes = index_genomes(fasta_file=genomic_fna, gff_file=genomic_gff)


    # chromosomes must be consistent between the fasta file and gff file
    # if a same chromosome is given in both chr_includes and chr_excludes, this chromosome will be included.
    valid_chromosomes = validate_chromosomes(to_include=chromosomes_to_process, to_exclude=chromosomes_to_exclude, in_gff=gff_indexes.keys(), in_fasta=fasta_hash.keys())
    print("Chromosomes that will be processed: {}".format(", ".join(valid_chromosomes)))

    chunks = get_chunks(l=valid_chromosomes, chunks_size=cpus)

    i = 0
    for chunk in chunks:

        # create and start processed
        processes = []

        for chr_id in chunk:
            basename = f"{basename_out}_{chr_id}_{i}"

            p = multiprocessing.Process(target=process_chromosome, args=(genomic_gff, fasta_hash[chr_id], chr_id, features, basename, out_formats, elongation, gff_indexes,))
            p.chr_id = chr_id
            processes.append(p)
            i += 1

        for p in processes:
            p.start()

        for p in processes:
            p.join()





if __name__ == '__main__':
    sys.exit(main())

    genomic_fna = "/home/nchenche/projects/orfmine_workdir/examples/database/Scer.fna"
    genomic_gff = "/home/nchenche/projects/orfmine_workdir/examples/database/Scer.gff"
