import sys
from orfmine.orftrack.lib import logHandler
from orfmine.utilities.lib.logging import get_logger


logger = get_logger(name=__name__)


def check_types(gff_data=None, types=None):

    unconsistent_types = []
    all_types = get_types(gff_data)

    for _type in types:
        if _type not in all_types:
            unconsistent_types.append(_type)

    if unconsistent_types:
        logger.error('Wrong argument type(s) has(have) been given:')
        for wrong_type in unconsistent_types:
            logger.error(' - ' + wrong_type)
        logger.error('')

        logger.error('You can choose amongst the types listed below:')
        for valid_type in all_types:
            logger.error(' - ' + valid_type)
        logger.error('')
        sys.exit(1)


def get_types(gff_data):
    all_types_2d = (gff_data[x].get_types() for x in (sorted(gff_data)))
    all_types_flatten = (sorted(set([val for sublist in all_types_2d for val in sublist])))

    return all_types_flatten


def check_chrids(chrs_gff: list = None, chrs_fasta: list = None):
    chr_common = set(chrs_gff).intersection(chrs_fasta)

    if chr_common:
        if len(chr_common) != len(chrs_fasta):
            table_chrs(chrs_gff, chrs_fasta)
            warning_log = 'All chromosomes are not shared between GFF and fasta files.' + \
                          'The process will continue with shared chromosome IDs only.'
            logger.warning(warning_log)
        else:
            table_chrs(chrs_gff, chrs_fasta)
            logger.info('All chromosomes are shared between GFF and fasta files.')
        return sorted(chr_common)
    else:
        table_chrs(chrs_gff, chrs_fasta)
        logger.error('Chromosomes are not consistent between GFF and fasta files.')
        sys.exit(1)


def table_chrs(chrs_gff, chrs_fasta):

    table_header = ["Chromosome ids", "in GFF", "in fasta"]
    spacer_len = 20
    table_border = spacer_len * len(table_header) * '-'
    row_format = ('{:>' + str(spacer_len) + '}') * (len(table_header))
    logger.info(row_format.format(*table_header))
    logger.info(table_border)

    all_chrs = sorted(set(chrs_gff + chrs_fasta))
    for chromosome in all_chrs:
        if chromosome in chrs_gff and chromosome in chrs_fasta:
            logger.info(row_format.format(chromosome, 'X', 'X'))
        elif chromosome in chrs_gff and chromosome not in chrs_fasta:
            logger.info(row_format.format(chromosome, 'X', '-'))
        elif chromosome in chrs_fasta and chromosome not in chrs_gff:
            logger.info(row_format.format(chromosome, '-', 'X'))
    logger.info('')
