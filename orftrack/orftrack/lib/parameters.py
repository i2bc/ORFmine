# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:37:10 2020

@author: nicolas
"""

import os,re
import argparse
from orftrack.lib import logHandler

logger = logHandler.Logger(name=__name__)

class Param:
    """
    Wrapper of all input files and arguments
    """
    def __init__(self, args):
        self.fasta_fname = args.fna
        self.gff_fname = args.gff

        self.chr = []
        self.chr_exclude = []
        if args.chr:
            self.chr = list(set(args.chr))
        if args.chr_exclude:
            self.chr_exclude = list(set(args.chr_exclude))

        self.types_only = []
        self.types_except = []
        if args.types_only:
            self.types_only = list(set(['CDS'] + args.types_only))
            self.types_except = []
        elif args.types_except:
            self.types_except = list(set(args.types_except))
            self.types_only = []

        self.o_include = args.o_include
        self.o_exclude = args.o_exclude
        self.orf_len = args.orf_len
        self.co_ovp = args.co_ovp

        self.bool_types = args.bool_types
        self.bool_chrs = args.bool_chrs
        self.is_frag = args.bool_isfrag
        self.o_fasta = args.bool_ofasta

        self.outpath = args.out if args.out.endswith('/') else args.out + '/'
        os.makedirs(self.outpath, exist_ok=True)
        self.default_basename = 'mapping_orf_'
        self.default_mainname = '.'.join(os.path.basename(self.fasta_fname).split('.')[:-1])
        self.outfile = self.outpath + self.default_basename + self.default_mainname

    def description(self):
        """

        A formatted string describing all key index positions stored.

        Returns:
            object: str

        """
        chr_asked = 'None' if not self.chr else ', '.join(self.chr)
        chr_exclude = 'None' if not self.chr_exclude else ', '.join(self.chr_exclude)
        logger.info('')
        logger.info('Parameters description:')
        logger.info('- fasta filename: ' + self.fasta_fname)
        logger.info('- gff filename: ' + self.gff_fname)
        logger.info('- chr: ' + chr_asked)
        logger.info('- chr_exclude: ' + chr_exclude)
        logger.info('- types_only: ' + ', '.join(self.types_only))
        logger.info('- types_except: ' + ', '.join(self.types_except))
        logger.info('- o_include: ' + ', '.join(self.o_include))
        logger.info('- o_exclude: ' + ', '.join(self.o_exclude))
        logger.info('- orf_len: ' + str(self.orf_len))
        logger.info('- co_ovp : ' + str(self.co_ovp))
        logger.info('- outfile: ' + self.outfile)
        logger.info('- bool_types: ' + str(self.bool_types))
        logger.info('- bool_chrs: ' + str(self.bool_chrs))
        logger.info('- isfrag: ' + str(self.is_frag))
        logger.info('- ofasta: ' + str(self.o_fasta))
        logger.info('')


def get_args():
    """
    Returns:
        Parameters
    """

    parser = argparse.ArgumentParser(description='Genomic mapping of pseudo-ORF')
    parser._action_groups.pop()  # Remove original optional argument group

    mandatory_arguments = parser.add_argument_group('Mandatory arguments')
    optional_arguments = parser.add_argument_group('Optional arguments')

    mandatory_arguments.add_argument("-fna", required=True, nargs="?",
                        help="Genomic fasta file (.fna) ")
    mandatory_arguments.add_argument("-gff", required=True, nargs="?",
                        help="GFF annotation file (.gff)")

    optional_arguments.add_argument("-chr", required=False, nargs="+", type=str, default=[],
                        help="List of seqID(s) - 1st column in the GFF file, generally chromosome or contig ID - \
                          to be treated by ORFtrack (all seqIDs are treated by default). \
                          The seqIDs must be separated by a space: -chr NC_001148.4 NC_001139.3")

    optional_arguments.add_argument("-chr_exclude", required=False, nargs="+", type=str, default=[],
                        help="List of seqID(s) you want to exclude (None by default.)")
    optional_arguments.add_argument("-types_only", required=False, nargs="+", default=[],
                        help="List of feature type(s) to use as reference(s) ('CDS' is included by default).")
    optional_arguments.add_argument("-types_except", required=False, nargs="+", default=['gene', 'exon'],
                        help="List of feature type(s) to not consider as reference(s) ('gene' and 'exon' by default).")

    # optional_arguments.add_argument("-o_include", required=False, nargs="+", default=['all'],
    #                     help="Type feature(s) and/or Status attribute(s) desired to be written in the output (all by default).")
    optional_arguments.add_argument("-o_include", required=False, nargs="+", default=['all'],
                        help=argparse.SUPPRESS)
    # optional_arguments.add_argument("-o_exclude", required=False, nargs="+", default=[],
    #                     help="Type feature(s) and/or Status attribute(s) desired to be excluded (None by default).")
    optional_arguments.add_argument("-o_exclude", required=False, nargs="+", default=[],
                        help=argparse.SUPPRESS)

    optional_arguments.add_argument("-orf_len", required=False, nargs="?", default=60, type=int,
                        help="Minimum number of coding nucleotides required to define a sequence between two consecutive stop codons\
                         as an ORF sequence (60 coding nucleotides by default).")
    optional_arguments.add_argument("-co_ovp", required=False, nargs="?", default=0.7, type=float,
                        help="Overlapping cutoff definition. An ORF sequence is considered overlapping with a genomic feature if \
                         its sequence overlaps with a fraction equals or above co_ovp (By default, co_ovp=0.7).")

    'Cutoff defining the minimum ORF sequence fraction'
    optional_arguments.add_argument("-out", required=False, nargs="?", default='/database/', type=str,
                        help="Output directory")

    optional_arguments.add_argument('--show-types', action='store_true', default=False,
                        dest='bool_types',
                        help='Print all feature types by seqIDs')
    optional_arguments.add_argument('--show-chrs', action='store_true', default=False,
                        dest='bool_chrs',
                        help='Print all seqIDs')

    # optional_arguments.add_argument('--frag-cds', action='store_true', default=False,
    #                     dest='bool_isfrag',
    #                     help='Generates fragments for CDS extremities that respect orf_len parameter.')
    optional_arguments.add_argument('--frag-cds', action='store_true', default=False,
                        dest='bool_isfrag',
                        help=argparse.SUPPRESS)

    optional_arguments.add_argument('--ofasta', action='store_true', default=False,
                        dest='bool_ofasta',
                        help='Writes amino acid and nucleic fasta files for ORFs')

    args = parser.parse_args()
    prefix = "/database/"
    if not re.match(prefix, args.fna):
        args.fna = prefix + args.fna
    if not re.match(prefix, args.gff):
        args.gff = prefix + args.gff
    if 'args.out' in locals():
        if not re.match(prefix, args.out):
            args.out = prefix + args.out

    return args
