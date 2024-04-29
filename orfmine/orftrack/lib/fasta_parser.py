# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:59:14 2020

@author: nicolas
"""
# import random
from json import load as json_load
from pathlib import Path
from pkg_resources import resource_filename as pkg_resource_filename

# from orfmine.orftrack.lib import logHandler
from orfmine.utilities.lib.logging import get_logger


logger = get_logger(name=__name__)


PATH_TO_CODON_TABLES = pkg_resource_filename('orfmine.utilities', 'data')


STANDARD_CODON_TABLE = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
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


class Fasta:
    """
    Class used to read/get informations of a chromosome from a genomic fasta file.
    """

    base_complement = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
        'R': 'Y',
        'Y': 'R',
        'W': 'N',
        'S': 'N',
        'M': 'N',
        'K': 'N',
        'B': 'N',
        'D': 'N',
        'H': 'N',
        'V': 'N',
        'X': 'N',
    }

    def __init__(self, codon_table_id: str="1"):
        """

        All key index file positions necessary to read a chromosome from a genomic fasta file
        without the need to store it in memory.

        fasta_fname (str): fasta filename
        chr (str): chromosome id
        curpos_start (int): index file position of the first nucleotide
        curpos_end (int): index file position of the last nucleotide
        seq_len (int): length of a nucleotide line in the fasta file
        line_len (int): length of a nucleotide line with extra-characters
        off_char (int): difference between seq_len and line_len
        codon_table_id (str): codon table id consistent with the NCBI codon tables

        """
        self.logger = get_logger(name=f"{__name__}.{self.__class__.__name__}")

        # self.logger = logHandler.Logger(name=f"{__name__}.Fasta")

        self.fasta_fname = None
        self.chr = None
        self.curpos_start = None
        self.curpos_end = None
        self.seq_len = None
        self.line_len = None
        self.off_char = None
        self.nucid_max = None
        self.codon_table = self.get_codon_table(codon_table_id)

    
    def get_codon_table(self, codon_table_id: str):
        codon_table_name = Path(PATH_TO_CODON_TABLES) / f"codon_table_{codon_table_id}.json"
        try:
            with open(Path(PATH_TO_CODON_TABLES / codon_table_name)) as table_file:
                codon_table = json_load(table_file)
        except:
            self.logger.error(f"Unable to initialize the codon table id from {codon_table_name}")
            self.logger.info("The standard codon table will be used")
            codon_table = STANDARD_CODON_TABLE

        return codon_table

    def init_nucid_max(self):
        """

        Returns the last nucleotide ID position of the chromosome.

        Returns:
            object: int

        """
        if not self.nucid_max:
            len_pos = self.curpos_end - self.curpos_start + 1
            n_line = int(len_pos / self.line_len)

            self.nucid_max = n_line * self.seq_len + len_pos % self.line_len

    def sequence(self, start=1, end=5, phase=0, strand='+'):
        """

        Returns the sequence corresponding to the given coordinates.
        The start position must be at least 1. The end position must not exceed nucid_max.

        Args:
            start (int): start position of the desired sequence
            end (int): end position of the desired sequence (int)
            phase (int): number of nucleotide to "remove" in the first codon (0, 1 or 2)
            (http://gmod.org/wiki/GFF#Nesting_Features)
            strand (str): strand of the nucleotide sequence ('+' or '-')

        Returns:
            object: str

        """
        start = start if start >= 1 else 1
        end = end if end <= self.nucid_max else self.nucid_max

        line_start = self.get_line_nucindex(index=start)
        line_end = self.get_line_nucindex(index=end)
        lines = self.get_lines(_from=line_start, to=line_end)

        len_sequence = end - start + 1
        start_pos = start % self.seq_len if start % self.seq_len else self.seq_len
        end_pos = start_pos - 1 + len_sequence

        if strand == '+':
            return ''.join(lines)[start_pos - 1 + phase:end_pos].upper()
        else:
            return self.reverse_complement(''.join(lines)[start_pos - 1:end_pos - phase]).upper()

    def get_line_nucindex(self, index=1):
        """

        Defines the line position where to find a given nucleotide index

        Args:
            index: nucleotide number (int)

        Returns:
            line_idx: line position where to find the nucleotide (int)

        """
        line_idx_raw = index / self.seq_len
        line_idx = int(line_idx_raw) if line_idx_raw % 1 == 0 else int(line_idx_raw) + 1

        return line_idx

    def get_lines(self, _from=1, to=1):
        """

        Returns all lines included between the given indexes

        Args:
            _from: index of the desired starting line (int)
            to: index of the desired last line (int)

        Returns:
            lines: a list of strings

        """
        n_lines = to - _from + 1
        lines = []
        with open(self.fasta_fname, 'rb') as fasta_file:
            fasta_file.seek(self.curpos_start + self.line_len * (_from-1))
            for n in range(n_lines):
                lines.append(fasta_file.readline().decode().strip())

        return lines

    def reverse_complement(self, sequence: str) -> str:
        """

        Returns the reverse complementary sequence of a given nucleotide sequence

        Arguments:
            - sequence: nucleotide sequence (str)

        Returns:
            object: str

        """
        sequence = list(sequence)
        sequence.reverse()

        return ''.join([self.base_complement[x.upper()] for x in sequence])

    def translate(self, start=1, end=10, strand='+', phase=0):
        """

        Translates a nucleotide sequence from its coordinates

        Arguments:
            start (int): start position of the desired sequence
            end (int): end position of the desired sequence (int)
            phase (int): number of nucleotide to "remove" in the first codon (0, 1 or 2)
            (http://gmod.org/wiki/GFF#Nesting_Features)
            strand (str): strand of the nucleotide sequence ('+' or '-')

        Returns:
            - object: str

        """

        if isinstance(phase, int):
            offset = (3 - ((end - start + 1 - phase) % 3)) % 3

            if strand == '+':
                end = end + offset if end + offset <= self.nucid_max else end
            elif strand == '-':
                start = start - offset if start - offset > 0 else start

        else:
            phase = 0

        sequence = self.sequence(start=start, end=end, phase=phase, strand=strand)

        protein = ''
        codons = (sequence[i:i + 3] for i in range(0, len(sequence), 3))

        for codon in codons:
            if codon.upper() not in self.codon_table.keys():
                protein += 'X'
            else:
                protein += self.codon_table[codon.upper()]

        return protein

    def index_resume(self):
        """

        A formatted string describing all key index positions stored.

        Returns:
            object: str

        """
        tab = '{:12}' * 7
        return tab.format(self.chr, self.curpos_start, self.curpos_end,
                          self.seq_len, self.line_len, self.off_char, self.nucid_max)


def parse(fasta_filename, codon_table_id: str="1"):
    """
    Reads a fasta file and stores key index positions to parse it without the need to store the whole file in memory.

    Args:
        fasta_filename: fasta filename (str)

    Returns:
        A list of FastaIndex instances.

    """
    chr_indexes = []
    with open(fasta_filename, 'rb') as fasta_file:
        for line in fasta_file:
            if line.startswith(b'>'):
                chr_index = Fasta(codon_table_id)
                chr_index.fasta_fname = fasta_filename
                chr_index.chr = line.decode().strip().split('>')[-1].split()[0]
                logger.debug('Reading chromosome: ' + chr_index.chr)
                chr_index.curpos_start = fasta_file.tell()

                seqline = fasta_file.readline()
                chr_index.line_len = len(seqline)
                chr_index.seq_len = len(seqline.decode().strip())
                chr_index.off_char = chr_index.line_len - chr_index.seq_len
                fasta_file.seek(-len(seqline), 1)

                if chr_indexes:
                    chr_indexes[-1].curpos_end = chr_index.curpos_start - len(line) - chr_index.off_char - 1
                    chr_indexes[-1].init_nucid_max()

                chr_indexes.append(chr_index)

        chr_indexes[-1].curpos_end = fasta_file.tell() - chr_indexes[-1].off_char - 1
        chr_indexes[-1].init_nucid_max()

    # for chr_index in chr_indexes:
    #     print(chr_index.index_resume())

    return {x.chr: x for x in chr_indexes}
