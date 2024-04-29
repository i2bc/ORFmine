#!/usr/bin/env python3

""" Some util functions related to biological DNA/RNA sequences """


# Constants:
DNA = b"ATGC"
RNA = b"AUGC"
DNA_C = b"TACG"
RNA_C = b"UACG"


def reverse_complement(sequence):
    """ Return the reverse complement of 'sequence'.
    Sequence Can be DNA (T) or RNA (U) string.
    """
    to_rev_comp = sequence.upper()
    if "U" in to_rev_comp:
        return to_rev_comp.translate(
            bytes.maketrans(RNA, RNA_C))[::-1]
    else:
        return to_rev_comp.translate(
            bytes.maketrans(DNA, DNA_C))[::-1]

def translate(sequence):
    """ Translate the nucleotide sequence to protein sequence with standard
    genetic code.
    """
    if len(sequence) % 3 != 0:
        print("Not a multiple of 3, sequence can't be translated: %s" % (
            sequence))
        return ""

    checked_seq = sequence.upper()
    gen_code = _genetic_code_dictionary()
    if "T" in checked_seq:
        # Use RNA version of the sequence to translate
        checked_seq = checked_seq.replace("T", "U")

    codons = [checked_seq[i:i+3] for i in range(0, len(checked_seq), 3)]
    protein_aminoacids = [gen_code[codon] for codon in codons]
    return "".join(protein_aminoacids)

def _genetic_code_dictionary():
    """ Return a dictionary with codons as key and the associated 1 letter code
    of amino acid as value.
    """
    stop_symbol = "*"
    codon_map = {
        "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
        "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
        "UAU":"Y", "UAC":"Y", "UAA":stop_symbol, "UAG":stop_symbol,
        "UGU":"C", "UGC":"C", "UGA":stop_symbol, "UGG":"W",
        "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
        "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
        "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
        "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
        "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    }
    return codon_map
