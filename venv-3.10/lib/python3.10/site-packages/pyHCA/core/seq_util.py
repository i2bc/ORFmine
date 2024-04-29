#!/usr/bin/env python
""" sequences handling functions
"""

import os, sys, argparse, string
from Bio import Seq
from Bio.SubsMat import MatrixInfo

__author__ = "Tristan Bitard-Feildel, Guillem Faure"
__licence__= "MIT"
__version__ = 0.1
__email__ = "t.bitard-feildel [you know what] upmc.fr.de"
__institute__ = "UPMC"

# problem of the SeqIO module: not memory efficient

__all__ = ["itercodon", "six_frames", "transform_seq"]

def transform_seq(seq):
    return seq.replace("*", "").replace("-", "").replace("?", "").replace("!","").replace(".", "")

def has_illegal_char(seq):
    illegal_char = set(["*", "-", "?", "!", "."])
    for c in seq:
        if c in illegal_char:
            return True
    return False

def has_gap_char(seq):
    illegal_char = set(["-", "."])
    for c in seq:
        if c in illegal_char:
            return True
    return False

def check_if_msa(sequences):
    sizes = set()
    illegal_char = False
    if isinstance(sequences, dict):
        for rec in sequences:
            seq = sequences[rec]
            sizes.add(len(seq))
            if has_gap_char(seq):
                illegal_char = True
    elif isinstance(sequences, list):
        for seq in sequences:
            sizes.add(len(seq))
            if has_gap_char(seq):
                illegal_char = True
    elif isinstance(sequences, str):
        # only one sequence
        illegal_char = False
    else:
        raise ValueError("Unknown argument type passed to check_if_msa(), {}").format(type(sequences))
    is_msa = False
    if sizes != set():
        ofsame_size = len(sizes) == 1 
        if illegal_char and not ofsame_size:
            print("Warning, multiple sequences of different lengths found with MSA symbols. Sequences of MSA must be of the same length", file=sys.stderr)
            sys.exit(1)
        elif illegal_char and ofsame_size:
            is_msa = True
    return is_msa

def itercodon(seq, frame, offset, table, reverse=False):
    stop = 0
    if not reverse:
        for i in xrange(frame, len(seq)-offset, 3):
            subseq = str(seq.seq)[i:i+3]
            assert(len(subseq)%3==0),(str(seq))
            aa = Seq.translate(subseq, table)
            yield i, aa
        if i+3 != len(seq):
            subseq = seq[i+3:] + "N"*(3-offset)
            assert(len(subseq)%3==0)
            aa = Seq.translate(subseq, table)
            yield i, aa
    else:
        for i in xrange(len(seq), offset, -3):
            # the reverse complement
            subseq = Seq.reverse_complement(str(seq.seq)[i-3:i])
            assert(len(subseq)%3==0)
            aa = Seq.translate(subseq, table)
            yield i, aa
        if offset:
            subseq = Seq.reverse_complement("N"*(3-offset) + str(seq.seq)[:offset])
            assert(len(subseq)%3==0)
            aa = Seq.translate(subseq, table)
            yield i, aa
    
#modulo = 0
#frame 0 0
#frame 1 -2
#frame 2 -1

#modulo = 1
#frame 0 -1
#frame 1  0
#frame 2 -2 

#modulo = 2
#frame 0 -2
#frame 1 -1
#frame 2  0

def six_frames(seq, table=1):
    """ Return the six frame translation of Seq protein object
    """
    offset = len(seq) % 3
    if len(seq) % 3 == 0:
        offset = [0, 2, 1]
    elif len(seq) % 3 == 1:
        offset = [1, 0, 2]
    else:
        offset = [2, 1, 0]
    start = 0
    reverse, prev = False, False
    for strand in [1, -1]:
        if strand < 0:
            reverse = True
        for frame in range(3):
            subprot = ""
            for nuc_idx, aa in itercodon(seq, frame, offset[frame], table, reverse):
                #if frame == 1 and strand == 1:
                    #print (start, nuc_idx, aa)
                if aa == "*" and subprot:
                    yield strand, frame, start, subprot
                    subprot = ""
                    prev = True
                elif prev:
                    start = nuc_idx
                    prev = False
                else:
                    subprot += aa
            if subprot:
                yield strand, frame, start, subprot
    
    
def compute_offset_pos(msa_seq, seq):
    """ from a sequence without gap position, computes the corresponding position in a MSA
    
    Parameters
    ==========
    msa_seq : string
        the sequence from the MSA 
    seq : string
        the protein sequence
        
    Return
    ======
    offset_msa2seq : int
        the position from the MSA to the sequence
    offset_seq2msa : int
        the position from the sequence to the MSA
    """
    msa_characters = set(["-", "?", "!", "*", "."])
    offset_msa2seq = dict()
    offset_seq2msa = dict()
    k = 0
    pos = 0
    for k in range(len(msa_seq)):
        isgap = True
        if msa_seq[k] not in msa_characters:
            isgap = False
            offset_seq2msa[pos] = k
            pos += 1    
        offset_msa2seq[k] = (pos - 1, isgap)
    return offset_msa2seq, offset_seq2msa


def compute_conserved_positions(dfasta, dmsa, score_type=0, matrix_name="blosum62"):
    """ compute conservation of a column relative to a sequence position of a msa
    score_type, 0: identity score, 1: binarized similarity (1 if sim > 0 else 0)
    """ 
    dconserv_per_prot = dict()
    dconserv_per_col = dict()
    positions_msa2prot = dict()
    
    matrix = getattr(MatrixInfo, matrix_name)
    items = list(matrix.items())
    matrix.update(((b,a),val) for (a,b),val in items)
    
    records = list()
    seq_idx = dict()
    for rec in dfasta:
        records.append(rec)
        msa2seq, seq2msa = compute_offset_pos(dmsa.get(rec, dfasta[rec]), dfasta[rec])
        seq_idx[rec] = msa2seq
        dconserv_per_prot[rec] = [0] * len(dfasta[rec])
        positions_msa2prot[rec] = dict()
        
    nb_seq =  len(records)
    if nb_seq > 0:
        nb_cols = len(dmsa[records[0]])
        tot = (nb_seq * (nb_seq - 1)) / 2
        for c in range(nb_cols):
            dconserv_per_col[c] = 0
            for k in range(len(records)-1):
                record_k = records[k]
                i, isgapi = seq_idx[record_k][c]
                seqk = dmsa[record_k]
                if seqk[c] != "-":
                    # no gap in sequence k
                    for l in range(k+1, len(records)):
                        record_l = records[l]
                        j, isgapj = seq_idx[record_l][c]
                        seql = dmsa[record_l]
                        if seql[c] != "-":
                            if score_type == 0:
                                score = 1 if (seqk[c] == seql[c]) else 0
                            else: # score_type == 1:
                                score = 1 if matrix[(seqk[c], seql[c])] > 0 else 0
                            dconserv_per_prot[record_k][i] += score
                            dconserv_per_prot[record_l][j] += score
                            dconserv_per_col[c] += score
            dconserv_per_col[c] /= tot
        # normalize score by number of sequence
        for rec in records:
            for i in range(len(dconserv_per_prot[rec])):
                dconserv_per_prot[rec][i] /= (nb_seq-1)
    return dconserv_per_prot, dconserv_per_col, seq_idx

def compute_hca_conserved_positions(dfasta, dmsa):
    """ compute conservation of a column relative to a hydrophobic cluster position
    """ 
    from .HCA import HCA
    dconserv_per_prot = dict()
    dconserv_per_col = dict()
    positions_msa2prot = dict()
    
    # hca transformation
    seq_bin = dict()
    for prot in dfasta:
        hca = HCA(seq=dfasta[prot])
        clusters = hca.get_clusters()
        binseq = [0] * len(dfasta[prot])
        for clust in clusters:
            for i in range(clust.start, clust.stop):
                binseq[i] = 1
        seq_bin[prot] = binseq
    
    records = list()
    seq_idx = dict()
    for rec in dfasta:
        records.append(rec)
        msa2seq, seq2msa = compute_offset_pos(dmsa.get(rec, dfasta[rec]), dfasta[rec])
        seq_idx[rec] = msa2seq
        dconserv_per_prot[rec] = [0] * len(dfasta[rec])
        positions_msa2prot[rec] = dict()
        
    nb_seq =  len(records)
    if nb_seq > 0:
        nb_cols = len(dmsa[records[0]])
        tot = (nb_seq * (nb_seq - 1)) / 2
        for c in range(nb_cols):
            dconserv_per_col[c] = 0
            for k in range(len(records)-1):
                record_k = records[k]
                i, isgapi = seq_idx[record_k][c]
                seqk = dmsa[record_k]
                if seqk[c] != "-":
                    # no gap in sequence k
                    for l in range(k+1, len(records)):
                        record_l = records[l]
                        j, isgapj = seq_idx[record_l][c]
                        seql = dmsa[record_l]
                        if seql[c] != "-":
                            score = 1 if (seq_bin[record_k][i] == seq_bin[record_l][j] == 1) else 0
                            dconserv_per_prot[record_k][i] += score
                            dconserv_per_prot[record_l][j] += score
                            dconserv_per_col[c] += score
            dconserv_per_col[c] /= tot
        # normalize score by number of sequence
        for rec in records:
            for i in range(len(dconserv_per_prot[rec])):
                dconserv_per_prot[rec][i] /= (nb_seq-1)
    return dconserv_per_prot, dconserv_per_col, seq_idx
