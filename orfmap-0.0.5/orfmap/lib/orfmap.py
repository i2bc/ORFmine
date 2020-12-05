import os
from orfmap.lib import logHandler
from orfmap.lib import gff_parser

logger = logHandler.Logger(name=__name__)


def mapping(gff_data, param):

    if os.path.exists(param.outfile + '.gff'):
        os.remove(param.outfile + '.gff')
    if os.path.exists(param.outfile + '.faa'):
        os.remove(param.outfile + '.faa')
    if os.path.exists(param.outfile + '.fna'):
        os.remove(param.outfile + '.fna')

    with open(param.outfile + '.gff', "a+") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(param.fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(param.gff_fname))
        out_gff.write(header)

        if param.o_fasta:
            out_fasta = open(param.outfile + '.faa', "a+")
            out_nucleic = open(param.outfile + '.fna', "a+")
        else:
            out_fasta = None
            out_nucleic = None

        for chr_id in sorted(gff_data):
            logger.info('Reading chromosome {} ...'.format(chr_id))
            gff_chr = gff_data[chr_id]

            logger.info(' - ORF mapping and assignment')
            get_orfs(gff_chr=gff_chr, param=param, outfiles=[out_gff, out_fasta, out_nucleic])
            logger.info('')

        if param.o_fasta:
            out_fasta.close()
            out_nucleic.close()


def get_orfs(gff_chr, param, outfiles: list):
    max_subsequence_length = 1999998
    out_gff = outfiles[0]
    out_fasta = outfiles[1]
    out_nucleic = outfiles[2]
    orf_len = param.orf_len + 6
    pos = 0

    # loops on each possible frame (the negative frame is defined in "frame_rev")
    for frame in range(3):
        subsequences = (gff_chr.sequence(start=i, end=i+max_subsequence_length) for i in
                        range(1+frame, gff_chr.end, max_subsequence_length))

        frame_rev = (gff_chr.end % 3 - frame) % 3
        start_pos = None
        start_pos_rev = None
        end_pos_rev = None

        for n, sequence in enumerate(subsequences):

            logger.debug('    - reading frame {} of subsequence {}'.format(str(frame), str(n)))
            codons = (sequence[i:i + 3].upper() for i in range(0, len(sequence), 3) if len(sequence[i:i + 3]) == 3)

            start_pos = start_pos if start_pos else frame + 1

            for pos, codon in enumerate(codons, start=int(n*max_subsequence_length/3)):
                if codon in ['TAG', 'TGA', 'TAA']:
                    end_pos = pos * 3 + 1 + 2 + frame
                    if end_pos - start_pos + 1 >= orf_len:
                        orf = build_orf(gff_chr=gff_chr, strand='+', frame=frame, coors=(start_pos, end_pos),
                                        param=param)
                        write_outputs(out_fasta=out_fasta, out_gff=out_gff, out_nucleic=out_nucleic, orf=orf, param=param)

                    start_pos = end_pos - 2

                elif codon in ['CTA', 'TCA', 'TTA']:
                    if start_pos_rev is None:
                        start_pos_rev = pos * 3 + 1 + frame
                    else:
                        end_pos_rev = pos * 3 + 1 + 2 + frame
                        if end_pos_rev - start_pos_rev + 1 >= orf_len:
                            orf = build_orf(gff_chr=gff_chr, strand='-', frame=frame_rev, coors=(start_pos_rev, end_pos_rev),
                                            param=param)
                            write_outputs(out_fasta=out_fasta, out_gff=out_gff, out_nucleic=out_nucleic, orf=orf, param=param)

                        start_pos_rev = end_pos_rev - 2

        # adds coordinates of ORF in negative strand extremity
        if end_pos_rev:
            start_pos_rev = end_pos_rev - 2
            end_pos_rev = pos * 3 + 1 + 2 + frame
            if end_pos_rev - start_pos_rev + 1 >= orf_len:
                orf = build_orf(gff_chr=gff_chr, strand='-', frame=frame_rev, coors=(start_pos_rev, end_pos_rev),
                                param=param, extremity=True)
                write_outputs(out_fasta=out_fasta, out_gff=out_gff, out_nucleic=out_nucleic, orf=orf, param=param)


def build_orf(gff_chr, strand, frame, coors, param, extremity=False):
    start_pos = coors[0]
    end_pos = coors[1]

    orf = gff_parser.GffElement(fasta_chr=gff_chr.fasta_chr)
    orf.seqid = gff_chr.id_
    orf.source = gff_chr.source
    orf.strand = strand
    orf.frame = frame

    if strand == '+':
        if start_pos == frame + 1:
            orf.start = start_pos
        else:
            orf.start = start_pos + 3
        orf.end = end_pos

    else:
        orf.start = start_pos
        orf.end = end_pos - 3 if not extremity else end_pos

    orf.run_assignment(elements=gff_chr.get_elements(coors=orf.get_coors()),
                       param=param, is_fragment=False)

    if param.is_frag:
        if orf.type == "c_CDS" and orf.ovp_phased:
            orf.fragment_phased_cds(orf_len=param.orf_len)
            if orf.suborfs:
                for suborf in orf.suborfs:
                    suborf.run_assignment(elements=gff_chr.get_elements(coors=suborf.get_coors()),
                                          param=param, is_fragment=True)

    return orf


def write_outputs(out_fasta, out_gff, out_nucleic, orf, param):
    if is_orf_asked(orf=orf, param=param):
        out_gff.write(orf.get_gffline())
        if param.o_fasta:
            out_fasta.write(orf.get_fastaline())
            out_nucleic.write(orf.get_fastanuc_line())
    if orf.suborfs:
        for suborf in orf.suborfs:
            if is_orf_asked(orf=suborf, param=param):
                out_gff.write(suborf.get_gffline())
                if param.o_fasta:
                    out_fasta.write(orf.get_fastaline())
                    out_nucleic.write(orf.get_fastanuc_line())


def is_orf_asked(orf=None, param=None):
    if 'all' in param.o_include:
        if not param.o_exclude:
            return True
        else:
            if not is_orf_exclude(orf=orf, exclude=param.o_exclude):
                return True
            else:
                return False
    else:
        if is_orf_include(orf=orf, include=param.o_include):
            if not is_orf_exclude(orf=orf, exclude=param.o_exclude):
                return True
            else:
                return False
        else:
            return False


def is_orf_include(orf=None, include=None):
    return True in [x in [orf.type, orf.status] for x in include]


def is_orf_exclude(orf=None, exclude=None):
    return True in [x in [orf.type, orf.status] for x in exclude]
