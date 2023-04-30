# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import sys
import time

from orfmine.orftrack.lib import mapping
from orfmine.orftrack.lib import logHandler
from orfmine.orftrack.lib import fasta_parser
from orfmine.orftrack.lib import gff_parser
from orfmine.orftrack.lib import parameters
from orfmine.orftrack.lib import tools


def main():
    # gets arguments
    start_time = time.time()

    params = parameters.set_parameters(args=parameters.get_args())
    # params = parameters.Param(args=parameters.get_args())

    if params.bool_chrs:
        tools.get_infos(_input=params.gff_fname, option='chrs')
        sys.exit(0)
    elif params.bool_types:
        tools.get_infos(_input=params.gff_fname, option='types')
        sys.exit(0)

    logger = logHandler.Logger(name=__name__, outpath=params.outpath)
    logo(logger)
    params.description()


    # parses fasta & gff by chromosomes
    logger.title('# Parsing GFF and fasta input files #')
    fasta_hash = fasta_parser.parse(fasta_filename=params.fasta_fname, codon_table_id=params.codon_table_id)
    gff_data = gff_parser.parse(param=params, fasta_hash=fasta_hash)

    # ORFs mapping (scans genome for stop-to-stop sequences and assigns them a status)
    logger.title('# Mapping ORFs (stop-to-stop codons) #')
    mapping.mapping(gff_data=gff_data, param=params)

    # Print a brief summary of ORFs mapping
    mapping.summary(gff_outfile=params.outfile+'.gff')

    logger.title("-- Execution time: {} seconds --".format(round((time.time() - start_time), 2)))


def logo(logger):
    logger.info('')
    logger.info('   ___    ____    _____   _                           _     ')
    logger.info('  / _ \  |  _ \  |  ___| | |_   _ __    __ _    ___  | | __ ')
    logger.info(' | | | | | |_) | | |_    | __| | \'__|  / _` |  / __| | |/ / ')
    logger.info(' | |_| | |  _ <  |  _|   | |_  | |    | (_| | | (__  |   <  ')
    logger.info('  \___/  |_| \_\ |_|      \__| |_|     \__,_|  \___| |_|\_\ ')
    logger.info('')
                                                            

if __name__ == '__main__':

    sys.exit(main())
