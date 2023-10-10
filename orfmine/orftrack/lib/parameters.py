# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:37:10 2020

@author: nicolas
"""


import argparse
import configparser
from pathlib import Path
from typing import List, Union
import sys

from orfmine.utilities.container import add_container_args
from orfmine.utilities.lib.logging import get_logger


class Parameter:
    """
    A class to encapsulate all input files and arguments for a genomic mapping with orftrack.

    Attributes
    ----------
    fna : Union[str, Path], optional
        Path to a genomic fasta file (.fna), by default None
    gff : Union[str, Path], optional
        Path to a GFF annotation file (.gff), by default None
    chr : List, optional
        List of seqID(s) to be treated by ORFtrack, by default empty list
    chr_exclude : List, optional
        List of seqID(s) to exclude, by default empty list
    types_only : List, optional
        List of feature type(s) to use as reference(s), by default empty list
    types_except : List, optional
        List of feature type(s) to not consider as reference(s), by default ['gene', 'exon']
    o_include : List, optional
        List of feature type(s) to include, by default ['all']
    o_exclude : List, optional
        List of feature type(s) to exclude, by default empty list
    orf_len : int, optional
        Minimum number of coding nucleotides required to define a sequence between two consecutive stop codons, by default 60
    co_ovp : float, optional
        Overlapping cutoff definition, by default 0.7
    bool_types : bool, optional
        Whether to print all feature types by seqIDs, by default False
    bool_chrs : bool, optional
        Whether to print all seqIDs, by default False
    bool_isfrag : bool, optional
        (Undocumented), by default False
    bool_ofasta : bool, optional
        Whether to write amino acid and nucleic fasta files for ORFs, by default False
    codon_table : str, optional
        Codon table ID to be used for considered chromosomes, by default '1'
    out : Union[str, Path], optional
        Output directory, by default './'

    """

    def __init__(
        self,
        fna: Union[str, Path]=None,
        gff: Union[str, Path]=None,
        chr: List=[],
        chr_exclude: List=[],
        types_only: List=[],
        types_except: List=['gene', 'exon'],
        o_include: List=['all'],
        o_exclude: List=[],
        orf_len: int=60,
        co_ovp: float=0.7,
        bool_types: bool=False,
        bool_chrs: bool=False,
        bool_isfrag: bool=False,
        bool_ofasta: bool=False,
        codon_table: str='1',
        out: Union[str, Path]='.',
    ):
        self.logger = get_logger(name=__name__)

        self.fasta_fname = fna
        self.gff_fname = gff

        self.chr = [] if chr is None else list(set(chr))
        self.chr_exclude = [] if chr_exclude is None else list(set(chr_exclude))

        self.types_only = []
        self.types_except = []
        if types_only:
            self.types_only = list(set(['CDS'] + types_only))
            self.types_except = []
        elif types_except:
            self.types_except = list(set(types_except))
            self.types_only = []

        self.o_include = o_include
        self.o_exclude = [] if o_exclude is None else o_exclude
        self.orf_len = orf_len
        self.co_ovp = co_ovp

        self.bool_types = bool_types
        self.bool_chrs = bool_chrs
        self.is_frag = bool_isfrag
        self.o_fasta = bool_ofasta

        self.codon_table_id = codon_table

        self.outpath = Path(out)
        self.outpath.mkdir(exist_ok=True, parents=True)
        
        self.default_basename = 'mapping_orf'
        self.default_mainname = Path(self.fasta_fname).stem
        self.outfile = str(self.outpath / f"{self.default_basename}_{self.default_mainname}")


    def description(self):
        """
        A formatted string describing all key index positions stored.

        # for arg, value in vars(args).items():
        #     print(f"{arg}: {value}")


        Returns:
            object: str
        """
        chr_asked = 'None' if not self.chr else ', '.join(self.chr)
        chr_exclude = 'None' if not self.chr_exclude else ', '.join(self.chr_exclude)

        self.logger.info('')
        self.logger.info('Parameters description:')
        self.logger.info('- fasta filename: ' + self.fasta_fname)
        self.logger.info('- gff filename: ' + self.gff_fname)
        self.logger.info('- codon table id: ' + self.codon_table_id)
        self.logger.info('- chr: ' + chr_asked)
        self.logger.info('- chr_exclude: ' + chr_exclude)
        self.logger.info('- types_only: ' + ', '.join(self.types_only))
        self.logger.info('- types_except: ' + ', '.join(self.types_except))
        self.logger.info('- o_include: ' + ', '.join(self.o_include))
        self.logger.info('- o_exclude: ' + ', '.join(self.o_exclude))
        self.logger.info('- orf_len: ' + str(self.orf_len))
        self.logger.info('- co_ovp : ' + str(self.co_ovp))
        self.logger.info('- outfile: ' + str(self.outfile))
        self.logger.info('- bool_types: ' + str(self.bool_types))
        self.logger.info('- bool_chrs: ' + str(self.bool_chrs))
        self.logger.info('- isfrag: ' + str(self.is_frag))
        self.logger.info('- ofasta: ' + str(self.o_fasta))
        self.logger.info('')


def set_parameters(args):
    if args.config:
        config_values = parse_config(args.config)
    else:
        config_values = {}

    provided_args = get_provided_args(parser=get_parser(), args=args, ignore_args=["--config"])
    for key, value in provided_args.items():
        config_values[key] = value

    params = Parameter(**config_values)

    return params


def get_provided_args(parser: argparse.ArgumentParser, args: argparse.Namespace, ignore_args: List=[]):
    provided_args = {}
    for action in parser._actions:
        if any(arg in sys.argv for arg in action.option_strings if arg not in ignore_args):
            provided_args[action.dest] = getattr(args, action.dest)
            
    return provided_args



def get_parser():
    parser = argparse.ArgumentParser(description='Genomic mapping of pseudo-ORF')
    parser._action_groups.pop()  # Remove original optional argument group

    mandatory_arguments = parser.add_argument_group('Mandatory arguments')
    optional_arguments = parser.add_argument_group('Optional arguments')

    mandatory_arguments.add_argument(
        "--fna", "-F",
        required=True,
        nargs="?",
        help="Genomic fasta file (.fna)"
    )

    mandatory_arguments.add_argument(
        "--gff", "-G",
        required=True,
        nargs="?",
        help="GFF annotation file (.gff)"
    )

    optional_arguments.add_argument(
        "--codon-table", "-T",
        required=False,
        type=str,
        default="1",
        help="Codon table ID to be used for considered chromosomes. \
            The ID must be consistent with the NCBI codon table IDs \
            (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). \
            By default, the standard codon table (id = 1) is used. Defaults to '1'"
    )

    optional_arguments.add_argument(
        "--chr", "-C",
        required=False,
        nargs="+",
        type=str,
        default=[],
        help="List of seqID(s) - 1st column in the GFF file, generally chromosome or contig ID - \
            to be treated by ORFtrack (all seqIDs are treated by default). \
            The seqIDs must be separated by a space: -chr NC_001148.4 NC_001139.3"
    )

    optional_arguments.add_argument(
        "--chr-exclude",
        required=False,
        nargs="+",
        type=str,
        default=[],
        help="List of seqID(s) you want to exclude (None by default.)"
    )

    optional_arguments.add_argument(
        "--types-only",
        required=False,
        nargs="+",
        default=[],
        help="List of feature type(s) to use as reference(s) ('CDS' is included by default)."
        )
    
    optional_arguments.add_argument(
        "--types-except",
        required=False,
        nargs="+",
        default=['gene', 'exon'],
        help="List of feature type(s) to not consider as reference(s) ('gene' and 'exon' by default)."
    )

    optional_arguments.add_argument(
        "--o-include",
        required=False,
        nargs="+",
        default=['all'],
        help=argparse.SUPPRESS
    )

    optional_arguments.add_argument(
        "--o-exclude",
        required=False,
        nargs="+",
        default=[],
        help=argparse.SUPPRESS
    )

    optional_arguments.add_argument(
        "--orf-len",
        required=False,
        default=60,
        type=int,
        help="Minimum number of coding nucleotides required to define a sequence between two consecutive stop codons\
            as an ORF sequence (60 coding nucleotides by default)."
    )

    optional_arguments.add_argument(
        "--co-ovp",
        required=False,
        nargs="?",
        default=0.7,
        type=float,
        help="Overlapping cutoff definition. An ORF sequence is considered overlapping with a genomic feature if \
            its sequence overlaps with a fraction equals or above co_ovp (By default, co_ovp=0.7)." 
    )

    optional_arguments.add_argument(
        "--out", "-O",
        required=False,
        nargs="?",
        default='.',
        type=str,
        help="Output directory ('.' by default)."
    )

    optional_arguments.add_argument(
        '--show-types',
        action='store_true',
        default=False,
        dest='bool_types',
        help='Print all feature types by seqIDs'
    )

    optional_arguments.add_argument(
        '--show-chrs',
        action='store_true',
        default=False,
        dest='bool_chrs',
        help='Print all seqIDs'
    )

    optional_arguments.add_argument(
        '--frag-cds',
        action='store_true',
        default=False,
        dest='bool_isfrag',
        help=argparse.SUPPRESS
    )

    optional_arguments.add_argument(
        '--ofasta',
        action='store_true',
        default=False,
        dest='bool_ofasta',
        help='Writes amino acid and nucleic fasta files for ORFs'
    )
    optional_arguments.add_argument(
        '--config',
        default=False,
        type=str,
        help='Path to a INI configuration file for orftrack'
    )

    parser = add_container_args(parser=parser)
    
    return parser


def get_args():
    """
    Returns:
        Parameters
    """

    parser = get_parser()
    args = parser.parse_args()
    
    return args


def parse_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    # Helper function to parse values and return None if empty
    def get_option(section, option, is_split=False):
        value = config.get(section, option)
        if value.strip() in ['', '""', "''"]:
            value = None

        if not value:
            return None
                
        return value.strip() if not is_split else value.split()

    # Parse mandatory arguments
    fna = get_option('mandatory', 'fna')
    gff = get_option('mandatory', 'gff')

    # Parse optional arguments
    codon_table = get_option('optional', 'codon_table')
    chr = get_option('optional', 'chr', is_split=True)
    chr_exclude = get_option('optional', 'chr_exclude', is_split=True)
    types_only = get_option('optional', 'types_only', is_split=True)
    types_except = get_option('optional', 'types_except', is_split=True)
    co_ovp = config.getfloat('optional', 'co_ovp')
    orf_len = config.getint('optional', 'orf_len')
    out = get_option('optional', 'out')
    show_types = config.getboolean('optional', 'show_types')
    show_chrs = config.getboolean('optional', 'show_chrs')
    ofasta = config.getboolean('optional', 'ofasta')
    o_include = get_option('optional', 'o_include', is_split=True)
    o_exclude = get_option('optional', 'o_exclude', is_split=True)
    frag_cds = config.getboolean('optional', 'frag_cds')

    # Store configuration values in a dictionary
    config_values = {}
    config_values['fna'] = fna
    config_values['gff'] = gff
    config_values['codon_table'] = codon_table
    config_values['chr'] = chr
    config_values['chr_exclude'] = chr_exclude
    config_values['types_only'] = types_only
    config_values['types_except'] = types_except
    config_values['co_ovp'] = co_ovp
    config_values['orf_len'] = orf_len
    config_values['out'] = out
    config_values['bool_types'] = show_types
    config_values['bool_chrs'] = show_chrs,
    config_values['bool_ofasta'] = ofasta
    config_values['o_include'] = o_include
    config_values['o_exclude'] = o_exclude
    config_values['bool_isfrag'] = frag_cds

    return config_values

