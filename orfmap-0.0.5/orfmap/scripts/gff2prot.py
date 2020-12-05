import os
import sys
import argparse
from orfmap.lib import gff_parser
from orfmap.lib import fasta_parser


def main():
    args = get_args()
    genomic_fna = args.fna
    genomic_gff = args.gff
    outpath = args.out if args.out.endswith('/') else args.out + '/'
    os.makedirs(outpath, exist_ok=True)
    if args.name:
        basename_out = outpath + args.name
    else:
        basename_out = outpath + os.path.basename('.'.join(genomic_fna.split('.')[:-1])) + '_proteins'

    fasta_hash = fasta_parser.parse(fasta_filename=genomic_fna)
    gff_data = parse_gff(gff_filename=genomic_gff, fasta_hash=fasta_hash)

    with open(basename_out + '.faa', 'w') as faa_file:
        with open(basename_out + '.fna', 'w') as fna_file:
            for chr_name in sorted(gff_data):
                proteins = gff_data[chr_name].group_cds()
                for protein in sorted(proteins):
                    header = '>' + protein + ':' + gff_data[chr_name].id_ + '\n'
                    fasta_aa = header + ''.join([cds.translate() for cds in proteins[protein]]) + '\n'
                    fasta_nuc = header + ''.join([cds.sequence() for cds in proteins[protein]]) + '\n'
                    faa_file.write(fasta_aa)
                    fna_file.write(fasta_nuc)


def parse_gff(gff_filename=None, fasta_hash=None):
    gff_data = {}
    with open(gff_filename, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('#'):
                chr_name = line.split('\t')[0]
                if chr_name not in gff_data:
                    gff_data[chr_name] = gff_parser.Chromosome(id_=chr_name, fasta_chr=fasta_hash[chr_name])
                    chromosome = gff_data[chr_name]
                    chromosome.source = line.split('\t')[1]

                element_type = line.split('\t')[2]
                if element_type == 'CDS':
                    chromosome.add(gff_element=gff_parser.GffElement(gff_line=line, fasta_chr=fasta_hash[chr_name]))

    return gff_data


def get_args():
    """
    Returns:
        Parameters
    """

    parser = argparse.ArgumentParser(description='Returns amino acid and nucleic fasta files from genomic data')
    parser.add_argument("-fna", required=True, nargs="?",
                        help="Genomic fasta file (.fna) ")
    parser.add_argument("-gff", required=True, nargs="?",
                        help="GFF annotation file (.gff)")
    parser.add_argument("-out", required=False, nargs="?", default='./', type=str,
                        help="Output directory")
    parser.add_argument("-name", required=False, nargs="?", default=None, type=str,
                        help="Basename for output files")

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    sys.exit(main())
