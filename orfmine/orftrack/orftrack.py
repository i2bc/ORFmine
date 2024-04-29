# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""
from argparse import Namespace
import logging
import sys
import time
from orfmine import DOCKER_IMAGE

from orfmine.orftrack.lib.fasta_parser import parse as fasta_parser
from orfmine.orftrack.lib.gff_parser import parse as gff_parser
from orfmine.orftrack.lib import mapping
from orfmine.orftrack.lib.parameters import get_args, set_parameters
from orfmine.orftrack.lib.tools import get_infos
from orfmine.utilities.container import ContainerCLI
from orfmine.utilities.lib.logging import get_logger, set_root_logger


logger = get_logger(name=__name__)


def run_orftrack(args):

    params = set_parameters(args=get_args())

    if params.bool_chrs:
        get_infos(_input=params.gff_fname, option='chrs')
        exit(1)
    elif params.bool_types:
        get_infos(_input=params.gff_fname, option='types')
        exit(1)

    # show set of used paramters 
    params.description()

    # parses fasta & gff by chromosomes
    logger.info('Parsing fasta file...')
    fasta_hash = fasta_parser(fasta_filename=params.fasta_fname, codon_table_id=params.codon_table_id)
    logger.info('Parsing gff file...')
    gff_data = gff_parser(param=params, fasta_hash=fasta_hash)


    # ORFs mapping (scans genome for stop-to-stop sequences and assigns them a status)
    logger.info('Mapping ORFs (stop-to-stop codons)...')
    mapping.mapping(gff_data=gff_data, param=params)

    # Print a brief summary of ORFs mapping
    mapping.summary(gff_outfile=params.outfile+'.gff')


def run_orftrack_containerized(args: Namespace):

    # list of flags related to input files
    input_args = ["--fna", "--gff"]
    # flag related to output path/file

    # update default config with from config file, if given; add input file if present
    if args.config:
        input_args += ["--config"]

    # instantiate containerCLI handler
    cli = ContainerCLI(
            input_args=input_args,
            output_arg="--out",
            args=args,
            image_base=DOCKER_IMAGE,
            prog="orftrack",
            container_type="docker" if args.docker else "singularity",
            dev_mode=args.dev,
            package_binding={"orfmine": "/home/orfuser/orfmine/orfmine"}
        )
    
    cli.show()
    if not args.dry_run:
        cli.run()


def main():

    # gets arguments
    args = get_args()

    if args.docker or args.singularity:
        run_orftrack_containerized(args=args)
    else:
        start_time = time.time()

        # set root logger
        set_root_logger(outpath=args.out, root="orfmine.orftrack", filename="orftrack.log", level=logging.INFO)
        run_orftrack(args=args)

        logger.info("-- Execution time: {} seconds --".format(round((time.time() - start_time), 2)))


if __name__ == '__main__':
    main()
