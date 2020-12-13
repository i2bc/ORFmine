# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import sys
from orfmap.lib import orfmap
from orfmap.lib import logHandler
from orfmap.lib import fasta_parser
from orfmap.lib import gff_parser
from orfmap.lib import parameters
from orfmap.lib import tools
import time


def main():
    # gets arguments
    # start_time = time.time()

    param = parameters.Param(args=parameters.get_args())

    if param.bool_chrs:
        tools.get_infos(_input=param.gff_fname, option='chrs')
        sys.exit(0)
    elif param.bool_types:
        tools.get_infos(_input=param.gff_fname, option='types')
        sys.exit(0)

    logger = logHandler.Logger(name=__name__, outpath=param.outpath)
    logo(logger)
    param.description()

    # sys.exit(0)

    # parses fasta & gff by chromosomes
    logger.title('# Parsing GFF and fasta input files #')
    fasta_hash = fasta_parser.parse(fasta_filename=param.fasta_fname)
    gff_data = gff_parser.parse(param=param, fasta_hash=fasta_hash, chr_asked=param.chr, chr_exclude=param.chr_exclude)

    # sys.exit(0)

    # ORFs mapping (scans genome for stop-to-stop sequences and assigns them a status)
    logger.title('# Mapping ORFs (stop-to-stop codons) #')
    orfmap.mapping(gff_data=gff_data, param=param)

    # Print a brief summary of ORFs mapping
    orfmap.summary(gff_outfile=param.outfile+'.gff')

    # logger.title("-- Execution time: {} seconds --".format(round((time.time() - start_time), 2)))


def logo(logger):
    logger.info('')
    logger.info('   ___    ___   ___   __  __               ')
    logger.info('  / _ \  | _ \ | __| |  \/  |  __ _   _ __ ')
    logger.info(' | (_) | |   / | _|  | |\/| | / _` | | ')
    logger.info('  \___/  |_|_\ |_|   |_|  |_| \__,_| | .__/')
    logger.info('                                     |')


if __name__ == '__main__':

    sys.exit(main())
