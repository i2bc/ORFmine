import logging
import os
from pathlib import Path
from shutil import copyfileobj as shutil_copyfileobj
from threading import Thread
import time
from typing import List, Union

import psutil

from orfmine.orftrack.lib import gff_parser, fasta_parser, inspect


class IndexFasta(Thread):
    def __init__(self, filename, codon_table_id: str="1"):
        super().__init__()
        self.filename = filename
        self.codon_table_id = codon_table_id
        self.result = None

    def run(self):
        print(f"{self.name} - Indexing fasta file...")
        start_time = time.time()
        self.result = fasta_parser.parse(fasta_filename=self.filename, codon_table_id=self.codon_table_id)
        print(f"{self.name} - Fasta file parsed in {round(time.time()-start_time, 2)} seconds")


class IndexGFF(Thread):
    def __init__(self, filename):
        super().__init__()
        self.filename = filename
        self.result = None

    def run(self):
        print(f"{self.name} - Indexing GFF file...")
        start_time = time.time()
        self.result = get_indexes_gff(gff_fname=self.filename)
        print(f"{self.name} - GFF file parsed in {round(time.time()-start_time, 2)} seconds")


class Concatenate(Thread):
    def __init__(self, inputs, output):
        super().__init__()
        self.inputs = inputs
        self.output = output

    def run(self):
        print(f"{self.name} - Merging intermediate output files...")
        start_time = time.time()
        concatfiles(inputs=self.inputs, output=self.output)
        print(f"{self.name} - Merging done in {round(time.time()-start_time, 2)} seconds")


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
    def __init__(self, check_stop_end: bool=False, check_stop_intra: bool=False) -> None:
        self.cds_list: List[gff_parser.GffElement] = []  # stores parent-related CDS unready
        self.stored_parent: str = ""  # stored last read parent
        self.cds_completed: List[gff_parser.GffElement] = []  # stores parent-related CDS ready to be merged; CDS in the list are added as completed only when the current parent is different from the stored/previous one
        self.check_stop_end = check_stop_end
        self.check_stop_intra = check_stop_intra
        self.check = True if check_stop_end or check_stop_intra else False

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

        if self.check:
            is_valid = self.has_protein_valid_stop()
            if not is_valid:
                self.cds_completed.clear()

        self.cds_list.clear()


    def has_protein_valid_stop(self):
        """Check the validity of a protein (merged CDS, aka completed) according to the presence of stop codons.
        
        If a protein has a stop codon (*):
            - at the end of its sequence (if self.check_stop_end) -> VALID if stop found at the end, INVALID otherwise
            - anywhere else in the sequence (if self.check_stop_intra) -> VALID if no stop found in the sequence, INVALID otherwise

        Return:
            bool: True if valid, False otherwise
        """
        if not self.cds_completed:
            return False
        
        aa_sequence = ''.join([cds.translate() for cds in self.cds_completed])
        if not aa_sequence:
            print(f"Error: amino acid sequence of protein {self.cds_completed[0].id_} is empty")
            return False
        
        if self.check_stop_end and aa_sequence[-1] != "*":
            return False
        
        if self.check_stop_intra and "*" in aa_sequence[:-1]:
            return False
        
        return True
    
    def get_fasta(self, _type="nucleic"):
        try:
            header = f">{self.cds_completed[0].id_}_mRNA\n"
        except:
            print(f"Exception, list error: {self.cds_completed}")

        if _type == "nucleic":
            return header + ''.join([cds.sequence() for cds in self.cds_completed]) + '\n'
        else:
            return header + ''.join([cds.translate() for cds in self.cds_completed]) + '\n'

    def get_gfflines(self):
        return "".join([cds.get_gffline() for cds in sorted(self.cds_completed, key=lambda x: x.start)])

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

    def build_fake_gff(self, elongate=0):
        gene_start = 1
        gene_end = sum([len(cds.sequence()) for cds in self.cds_completed]) + 2 * elongate
        cds_start = elongate + 1
        cds_end   = gene_end - elongate
        cds_template = self.cds_completed[0]

        gene_line = "{}_mRNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cds_template.id_, "elongated", "gene", gene_start, gene_end, ".", "+", ".", "ID=" + cds_template.id_ + "_gene")
        mrna_line = "{}_mRNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cds_template.id_, "elongated", "mRNA", gene_start, gene_end, ".", "+", ".","ID=" + cds_template.id_ + "_mRNA; Parent=" + cds_template.id_ + "_gene")
        utr5_line = "{}_mRNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cds_template.id_, "elongated", "five_prime_UTR", gene_start, elongate, ".", "+", ".","ID=" + cds_template.id_ + "_5UTR; Parent=" + cds_template.id_ + "_mRNA")
        cds_line = "{}_mRNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cds_template.id_, "elongated", "CDS", cds_start, cds_end, ".", "+", ".","ID=" + cds_template.id_ + "_CDS; Parent=" + cds_template.id_ + "_mRNA")
        utr3_line = "{}_mRNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cds_template.id_, "elongated", "three_prime_UTR", cds_end + 1, gene_end, ".", "+", ".","ID=" + cds_template.id_ + "_3UTR; Parent=" + cds_template.id_ + "_mRNA")
                
        return gene_line + mrna_line + utr5_line + cds_line + utr3_line


def memory_usage():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)

    return mem


def get_indexes_gff(gff_fname):
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


def get_chunks(l: list, chunks_size: int=10):     
    for i in range(0, len(l), chunks_size):
        yield l[i:i+chunks_size]


def validate_chromosomes(to_include: list=[], to_exclude: list=[], in_gff: list=[], in_fasta: list=[]) -> list:
    logger = logging.getLogger(__name__)

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


def index_genomes(fasta_file: str, gff_file: str, codon_table_id: str="1"):
    indexing_threads: List[Union[IndexGFF, IndexFasta]]

    indexing_threads = [IndexFasta(filename=fasta_file, codon_table_id=codon_table_id), IndexGFF(filename=gff_file)]
    for t in indexing_threads:
        t.start()

    for t in indexing_threads:
        t.join()

    return [x.result for x in indexing_threads]


def concatfiles(inputs: List[Union[str, Path]], output: Union[str, Path]):
    with open(output, 'w') as wfd:
        for f in inputs:
            with open(f,'r') as fd:
                shutil_copyfileobj(fd, wfd)
            Path(f).unlink()
            
        
def merge_outfiles(inpath: Union[str, Path], infiles: List[Union[str, Path]], extensions: List[str], outbasename: str):
    threads: List[Concatenate] = []
    for ext in extensions:
        output = inpath / f"{outbasename}{ext}"
        files_ext = [ f"{base_file}{ext}" for base_file in infiles if Path(f"{base_file}{ext}").exists() ]

        threads.append(Concatenate(inputs=files_ext, output=output))

    for t in threads:
        t.start()

    for t in threads:
        t.join()
