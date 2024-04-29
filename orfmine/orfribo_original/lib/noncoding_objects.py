#!/usr/bin/env python3

""" A. Lopes project: The non coding regions associated objects.

This module contains all objects relative to intergenic regions / ORFs, called
IGR and IGORF respectively.
"""

# Local
from .constants import *
from .biobjects import Feature
from .sequtils import translate, reverse_complement
from .utils import other_frames, find_stretches

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"

# Constants
STOPS = ["TAG", "TAA", "TGA"]
STOPS_REV = ["TTA", "CTA", "TCA"]

IG_COLS = {
    "IGR":"#2eb82e",
    "+":{
        0:"#ff4d4d",
        1:"#ff4d4d",
        2:"#ff4d4d"
    },
    "-":{
        0:"#3366ff",
        1:"#3366ff",
        2:"#3366ff"
    }
}

class Igr():
    """ IGR strands for InterGenic Region. """
    # Constructor
    def __init__(self, seqid, left_genomic, right_genomic, left_strand="",
                 right_strand="", left_frame=None, right_frame=None):
        """ constructor v1: strands are -1 ('-') and +1 ('+') """
        self.seqid = seqid
        self.left_genomic = left_genomic
        self.right_genomic = right_genomic
        self.left_strand = left_strand
        self.right_strand = right_strand
        self.left_frame = left_frame
        self.right_frame = right_frame
        self.length = (self.right_genomic - self.left_genomic) + 1

        self.ID = "{seqid}_{left}-{right}".format(
            seqid=self.seqid,
            left=self.left_genomic,
            right=self.right_genomic
        )
        # Sequence can be init latter calling self.get_sequence method
        self.sequence = None
        self.childs = []        # Will be a list of Igorf objects

    # Private methods
    def __lt__(self, other):
        """ sorted(igr_list) will sort igrs on their left_genomic. """
        return self.left_genomic < other.left_genomic

    # Public methods
    def tuple_view(self):
        """ Return a tupleview of the IGR. """
        return (self.left_genomic, self.right_genomic)

    def gffline(self):
        """ Return the gffline corresponding to the IGR """
        attributes = "ID={oid};color={col};length={ilen}".format(
            oid=self.ID,
            col=IG_COLS["IGR"],
            ilen=self.length
        )

        fields = [
            self.seqid,  # seqid
            "SOURCE",  # source
            "igr",  # type
            str(self.left_genomic),  # start (min genomic)
            str(self.right_genomic),  # stop (max genomic)
            ".",  # score
            "+",  # strand
            ".",  # phase
            attributes  # attributes
        ]
        return GFF_SEP.join(fields)

    def get_sequence(self, fapick):
        """ Get the IGR sequence, strand +, from the fasta pickle. """
        seq = ""
        try:
            seqid_sequence = fapick["sequences"][self.seqid]
        except KeyError as ke:
            print("%s not found in fapick file" % (ke))
            return seq
        seq = seqid_sequence[self.left_genomic-1 : self.right_genomic].upper()
        return seq


class Igorf(Feature):
    """ Initialize IGORFs from Igorf detection results. """
    def __init__(self, gffline):
        """ gffline required. """
        # Init with Feature constructor to access genome arithmetics
        super().__init__(gffline=gffline)

    # Public methods
    def frames_coverage(self, bamfile, shift=SHIFT, side=SIDE, kmers=KMERS,
                        rm_zero=True):
        """ Return a dictionary with the coverage of the 3 frames.
        rm_zero is used to remove 0 in other frames for a given position.
        """
        covs = {0:[], 1:[], 2:[]}
        pos_cov = self.pos_cov(bamfile, shift=shift, side=side, kmers=kmers)
        for (idx, cov) in enumerate(pos_cov["coverages"]):
            frame = idx % 3
            covs[frame].append(cov)
            if not rm_zero:
                for o_frame in other_frames(frame):
                    covs[o_frame].append(0)
        return covs


    def stretches(self, bamfile, shift=SHIFT, side=SIDE,
                  kmers=KMERS, frames=False, threshold=1,
                  min_len=10):
        """ Return the stretch of position with coverages < threshold.
        If frames is True, the stretches will be seeked in each frame coverage.
        """
        if not frames:
            pos_cov = self.pos_cov(bamfile, shift=shift, side=side, kmers=kmers)
            all_stretches = find_stretches(pos_cov["coverages"],
                                           threshold=threshold,
                                           min_len=min_len)
            return all_stretches
        else:
            frames_stretches = {0:[], 1:[], 2:[]}
            frames_covs = self.frames_coverage(bamfile,
                                               shift=shift,
                                               side=side,
                                               kmers=kmers)
            for frame in frames_covs:
                f_cov = frames_covs[frame]
                f_stretches = find_stretches(f_cov,
                                             threshold=threshold,
                                             min_len=min_len)
                frames_stretches[frame] = f_stretches
            return frames_stretches


    def score_coverage(self, bamfile, shift=SHIFT, side=SIDE,
                       kmers=KMERS, frames=False, threshold=1):
        """ The score_coverage is defined by the % of positions with a coverage
        >= threshold. If frames is True, this score will be computed for each
        frame separately.
        """
        if not frames:
            pos_cov = self.pos_cov(bamfile, shift=shift, side=side, kmers=kmers)
            count_gt_threshold = sum(1 for val in pos_cov["coverages"]
                                     if val >= threshold)
            score = (count_gt_threshold * 100) / len(pos_cov["coverages"])
            return round(score, 2)
        else:
            scores = {0:0, 1:0, 2:0}
            # rm_zero option is used to be consistent with 1/3
            # of positions for each frames
            frames_covs = self.frames_coverage(bamfile,
                                               shift=shift,
                                               side=side,
                                               kmers=kmers,
                                               rm_zero=True)
            for frame in frames_covs:
                f_cov = frames_covs[frame]
                count_gt_threshold = sum(1 for val in f_cov if val >= threshold)
                score = (count_gt_threshold * 100) / len(f_cov)
                scores[frame] = score
            return scores



class IgorfDetector():
    """ Store IGORFs attributes for the detection. """
    def __init__(self, igr_parent, start, stop, strand,
                 increment, sequence, absolute_frame):
        """ Init Igorf as 'strict' and 'IG' """
        # Init attributes
        self.igr_parent = igr_parent
        self.start = start
        self.stop = stop
        self.strand = strand
        self.increment = increment
        self.sequence = sequence
        self.absolute_frame = absolute_frame

        self.length = None
        self._compute_length()

        self.ID = "{seqid}_{st}_{sta}-{sto}_{fra}".format(
            seqid=self.igr_parent.seqid,
            st=self.strand,
            sta=self.start,
            sto=self.stop,
            fra=self.absolute_frame
        )
        # Attribute to set while detecting igorfs
        self.overlaps_coding_blocks = False
        self.location = "IG" # default = intergenic, no overlap
        self.other_igorfs_overlapped = []
        self.status = "strict" # IG or overlap_size < minlen of ORFs.
        self.keep_me = True    # if length > minlen of ORFs
        self.left_ov_size = 0
        self.right_ov_size = 0

    # Private methods
    def _compute_length(self):
        self.length = abs(self.stop - self.start) + 1

    def _adjust_length(function):
        """ Decorator to update length when update_right/left are called. """
        def update(self, *args, **kwargs):
            function(self, *args, **kwargs)
            self._compute_length()
        return update

    # Public methods
    def gffline(self):
        """ Return the GFF line correponding to the igorf. """

        if self.right_ov_size > 0 or self.left_ov_size > 0:
            overlap = "True"
        else:
            overlap = "False"

        attributes = "ID={oid};Parent={igr};color={col};left_ov_len={rlen};right_ov_len={llen};Location={loc};Status={status};Overlap={ov}".format(
            oid=self.ID,
            igr=self.igr_parent.ID,
            col=IG_COLS[self.strand][self.absolute_frame],
            rlen=self.right_ov_size,
            llen=self.left_ov_size,
            loc=self.location,
            status=self.status,
            ov=overlap
        )

        fields = [
            self.igr_parent.seqid, # seqid
            "SOURCE",              # source
            "igorf",               # type
            str(min(self.start, self.stop)), # start (min genomic)
            str(max(self.start, self.stop)), # stop (min genomic)
            ".",                        # score
            self.strand,                # strand
            str(self.absolute_frame),         # phase
            attributes                  # attributes
        ]

        return GFF_SEP.join(fields)

    @_adjust_length
    def update_left(self, fapick, location, minlen):
        """ Update the left of the sequence if overlapping
        with coding block.
        """
        seqid_seq = fapick["sequences"][self.igr_parent.seqid]
        self.location = location
        self.overlaps_coding_blocks = True

        if self.strand == "+":
            stops = STOPS
            idx = self.start - 1
        else:
            stops = STOPS_REV
            idx = self.stop - 1

        seq = ""
        while True:
            next_codon = seqid_seq[idx-3:idx]
            if next_codon in stops:
                seq = next_codon + seq
                break
            else:
                seq = next_codon + seq
                idx -= 3
        self.left_ov_size = len(seq)

        ## DEFINE STATUS HERE ##
        if self.length < minlen:
            self.status = "soft"

        if self.strand == "-":
            seq = reverse_complement(seq) # 3', contains stop
            self.sequence = self.sequence + seq
            self.stop -= len(seq)
        else:
            seq = seq[3:] # no stop as 5' of ORF
            self.sequence = seq + self.sequence
            self.start -= len(seq)

    @_adjust_length
    def update_right(self, fapick, location, minlen):
        """ Update the right of the sequence if overlapping
        with coding block.
        """
        seqid_seq = fapick["sequences"][self.igr_parent.seqid]
        self.location = location
        self.overlaps_coding_blocks = True

        # +1 already taken in account as annotation 1based and python idx 0based
        if self.strand == "+":
            stops = STOPS
            idx = self.stop
        else:
            stops = STOPS_REV
            idx = self.start

        seq = ""
        while True:
            next_codon = seqid_seq[idx:idx+3]
            #print(next_codon)
            if next_codon in stops:
                seq += next_codon
                break
            else:
                seq += next_codon
                idx += 3

        self.right_ov_size = len(seq)

        ## DEFINE STATUS HERE ##
        if self.length < minlen:
            self.status = "soft"

        if self.strand == "-":
            seq = reverse_complement(seq)
            seq = seq[3:] # Because strand - == 5', stop codon not included
            self.sequence = seq + self.sequence
            self.start += len(seq)
            # sequence is 5' => 3' oriented, so only need after the stop codon.
        else:
            self.sequence = self.sequence + seq
            self.stop += len(seq)
