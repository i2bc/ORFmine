import os

def write_orf(orf, param=None):
    """
    Writes fasta and gff files for orfs.

    Args:
        all_orfs: list of lib.gff_parser.Gff_element instances
        param: instance of lib.parameters.Param

    Returns:
        None

    """

    if os.path.exists(param.outfile + '.gff'):
        os.remove(param.outfile + '.gff')
    if os.path.exists(param.outfile + '.fa'):
        os.remove(param.outfile + '.fa')

    with open(param.outfile + '.gff', "a+") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(param.fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(param.gff_fname))
        out_gff.write(header)
        with open(param.outfile + '.fa', "a+") as out_fasta:
            for orf in all_orfs:
                out_gff.write(orf.get_gffline())
                out_fasta.write(orf.get_fastaline())



def write_orfs(all_orfs: list, param=None):
    """
    Writes fasta and gff files for orfs.

    Args:
        all_orfs: list of lib.gff_parser.Gff_element instances
        param: instance of lib.parameters.Param

    Returns:
        None

    """

    with open(param.outfile + '.gff', "w") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(param.fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(param.gff_fname))
        out_gff.write(header)
        with open(param.outfile + '.fa', "w") as out_fasta:
            for orf in all_orfs:
                out_gff.write(orf.get_gffline())
                out_fasta.write(orf.get_fastaline())


def write_orf(orf, param=None):
    """
    Writes fasta and gff files for orfs.

    Args:
        all_orfs: list of lib.gff_parser.Gff_element instances
        param: instance of lib.parameters.Param

    Returns:
        None

    """

    with open(param.outfile + '.gff', "w") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(param.fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(param.gff_fname))
        out_gff.write(header)
        with open(param.outfile + '.fa', "w") as out_fasta:
            for orf in all_orfs:
                out_gff.write(orf.get_gffline())
                out_fasta.write(orf.get_fastaline())
