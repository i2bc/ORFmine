#!/usr/bin/env python3

""" Biological entities handling.

biobjects are classes related to biological entities
such as Features, Transcripts and Genes.

"""

# Stdlib

# Local
from constants import *
from gffutils import get_tags_and_values, attributes_string_from_tags, check_rank_in_id, read_all_fields
from sambam import get_stranded_region_coverage_reduced, get_bank_size, get_stranded_region_counts
from sequtils import reverse_complement, translate

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"


class StraGeRa():
    """ StraGeRa: STRAnded GEnome RAnge.

    This class aims at handling a continuous nucleotide sequence
    of a genome. It may be defined by 4 attributes:
    > seqid (more used as chromosome)
    > start (5' end)
    > end (3' end)
    > strand ('+' or '-')

    Since a StraGeRa object will always be 5' => 3' oriented, the
    consistency of start/end and strand is automatically tested.
    Example:

        strag = StraGeRa("chrI", 10, 50, "-")
        print(strag.start)
        # output: 50, not 10.

    In contrast with 'start'/'end' (5' => 3' orientation), the 'left'/'right'
    attributes respectively refer to the minimum and maximum genomic positions
    covered by the StraGeRa. Those attributes are used in the GFF3/BED formats
    to define limits of the features.
    Example:
        strag = StraGeRa("chrI", 10, 50, "-")
        print(strag.left)
        # output: 10.
    """

    def __init__(self, seqid, start, end, strand):
        """ Constructor """
        self.seqid = seqid
        self.strand = strand
        self.left = min(start, end)
        self.right = max(start, end)

        try: # Check start and end types
            self.start = int(start)
            self.end = int(end)
        except ValueError as ve:
            raise ValueError(
                """StraGeRa 'start' and 'end' must be integers:
                {r}{ve}{n} can't be converted into <int>.""".format(
                    r=REIT,
                    n=NA,
                    ve=ve,
                ))
        try: # Strand must be '+' or '-'
            self.increment = STRANDS_MAP[strand]["one"]
        except KeyError as ke:
            raise KeyError(
                """{ke} is not a legal strand. '+' or '-' expected.""".format(
                    **locals()))
        self._check_start_end_order() # Check start/end and strand consistency

    # Private methods
    def __lt__(self, other):
        """ StraGeRa are sorted by their start position (5'). """
        return self.start < other.start

    def __repr__(self):
        """ Convenient visualisation on printing """
        return "<%s.%s seqid=%s, start=%s, end=%s, strand=%s>" % (
            __name__,
            self.__class__.__name__,
            self.seqid,
            self.start,
            self.end,
            self.strand
            )

    def _check_start_end_order(self):
        """ Prevent range() to fail keeping start, end
        and increment consistent.
        """
        if (self.increment == -1 and self.start < self.end) \
           or (self.increment == 1 and self.start > self.end):
            self.start, self.end = self.end, self.start

    # Public methods
    def length(self):
        """ Return the length of the StraGeRa. """
        length = abs(self.end - self.start) + 1
        return length

    def tupleview(self, strand="int"):
        """ Return a tuple view of the StraGeRa with:
        (start, end, increment/strand)
        MIN_COORDS is a namedtuple defined in constants.py.

        I think I need to remove it, and use StraGeRa objects instead?
        """
        if strand == "int":
            tupview = MIN_COORDS(self.start, self.end, self.increment)
        elif strand == "symbol":
            tupview = MIN_COORDS(self.start, self.end, self.strand)
        else: # If not valid strand, the symbol is returned
            tupview = MIN_COORDS(self.start, self.end, self.strand)
            print("""{r}{strand}{na} is not a valid strand option.
            'int' or 'symbol' excepted, 'symbol' was used instead.
            """.format(
                r=REIT,
                n=NA,
                strand=strand
            ))
        return tupview

    def adjust_limits(self, n):
        """ Return the adjused limits start/end depending on the number
        of positions required by 'n'. The increment is usefull
        to conserve the orientation (can be 1 or -1).
        """
        length = self.length()
        if abs(n) > length or not n:
            return (self.start, self.end)

        shift = length - abs(n)
        if n > 0:
            adj_start = self.start
            adj_end = self.end - (shift * self.increment)
        elif n < 0:
            adj_start = self.start + (shift * self.increment)
            adj_end = self.end

        return (adj_start, adj_end)

    def covered_positions(self, n=N_0):
        """ Return a generator of all positions covered by
        the StraGeRa (start and end are included).
        'n' can be used to select a sub-range. If n > 0, the
        positions are selected from the 5'. If n < 0, from the 3'.
        If n > self.length() or n = 0, n is ignored and all positions
        are returned.
        """
        adj_start, adj_end = self.adjust_limits(n)
        for pos in range(adj_start, adj_end + self.increment, self.increment):
            yield pos

    def pos_cov(self, bamfile, shift=SHIFT, side=SIDE, kmers=KMERS, n=N_0):
        """ Return a dictionary with the coverage informations with the keys:
        - 'positions' = all the genomic positions covered by StraGeRa (5'=>3')
        - 'coverages' = # of reads on each positions, same order as 'positions'
        - 'details' = list of dictionaries with the Kmer count details
        - 'bank_size' = total number of read mapped in the bamfile
        """
        adj_start, adj_end = self.adjust_limits(n)
        pos_cov = get_stranded_region_coverage_reduced(
            self.seqid,
            self.strand,
            adj_start,
            adj_end,
            bamfile,
            shift=shift,
            side=side,
            kmers=kmers
        )
        ordered_coverages = {   # 5' => 3' ordered
            "positions":[],
            "coverages":[],
            "details":[],
            "bank_size":get_bank_size(bamfile)
        }
        for pos in range(adj_start, adj_end + self.increment, self.increment):
            p_cov = pos_cov[pos] # p_cov contains details of Kmers
            ordered_coverages["positions"].append(pos)
            ordered_coverages["coverages"].append(p_cov[TOTAL_KEYNAME])
            ordered_coverages["details"].append(p_cov) # Keep kmers distrib
        return ordered_coverages

    def counts(self, bamfile, shift=SHIFT, side=SIDE, kmers=KMERS):
        """ Return a dictionary with the number of reads (each kmer detailed)
        on the StraGeRa positions.
        """
        counts = get_stranded_region_counts(
            self.seqid,
            self.strand,
            self.start,
            self.end,
            bamfile,
            shift=shift,
            side=side,
            kmers=kmers
        )
        return counts

    def sequence(self, fapick, n=N_0, translate_to_protein=False):
        """ Return the 5' -> 3' oriented sequence of the StraGeRa.
        `fapick` = <dict> with seqid as key and sequence as value.
        Use 'fasta2pickle.py' program to generate a fapick from Fasta.
        """
        adj_start, adj_end = self.adjust_limits(n)
        try:
            reference_sequence = fapick["sequences"][self.seqid]
        except KeyError as ke:
            print("%s seqid not found in fapick dictionary" % (ke))
            return ""

        if self.strand == "-":
            sequence = reference_sequence[adj_end-1 : adj_start].upper()
            stragera_seq = reverse_complement(sequence)
        else:
            stragera_seq = reference_sequence[adj_start-1 : adj_end].upper()

        if translate_to_protein:
            return translate(stragera_seq)
        else:
            return stragera_seq

    def mask(self, fapick):
        """ Mask the StraGeRa sequence in the fapick (fasta in pickle).
        /!\ Warning, fapick is modified inplace!
        """
        # -1 to adjust python 0 index to annotation 1 index
        genomic_left = min(self.start, self.end) - 1
        genomic_right = max(self.start, self.end) - 1

        seqid_seq = fapick["sequences"][self.seqid]

        stretch = MASK_BASE * self.length()
        seqid_seq = seqid_seq[:genomic_left] \
                    + stretch \
                    + seqid_seq[genomic_right+1 : ]
        fapick["sequences"][self.seqid] = seqid_seq

    # Need plot function here
    # Make the DATAFRAME => then use kmers or total vbar_figure
    def bokeh_figure(self, shift=SHIFT, side=SIDE, fapick=None,
                     only_total=False):
        """ Coverage representation of the StraGeRa. """
        # To be implemented: Need to be able to print sequence
        # + annot + kmers details or not.


class Feature(StraGeRa):
    """ Feature class can be seen as a line of GFF file.
    A Feature inherits from StraGeRa to access all genome
    arithmetic functions.
    """
    def __init__(self, feature_id="NO_ID", gffline=None, **kwargs):
        """ A Feature can be initialized in two ways:
        1) A line of GFF3 formatted file
        2) Manually from kwargs
        Feature_id is ONLY used for manual initialization. With gffline,
        the feature_id is initialised from attributes field ('ID' tag).
        """
        # Initialisation by gffline is prioritary
        if gffline:
            self._init_from_gffline(gffline)
        else:
            self._init_manually(feature_id, kwargs)

        super().__init__(self.seqid, self.start, self.end, self.strand)
        self.tags = get_tags_and_values(self.attributes)

        try: # Check ID of the feature, which is supposed to be MANDATORY
            feat_id = self.tags[ID_TAG]
        except KeyError: # If no ID, generate an unique ID
            feat_id = self._generate_id()
        self.set_tag(ID_TAG, feat_id)

        try:
            self.rank = int(self.get_tag(RANK_TAG))
        except TypeError: # get_tag returned None
            self.rank = check_rank_in_id(self.get_tag(ID_TAG))

    # Private methods
    def __repr__(self):
        """ Convenient visualisation on printing """
        return "<%s.%s ID=%s seqid=%s, start=%s, end=%s, strand=%s>" % (
            __name__, self.__class__.__name__,
            self.ID, self.seqid,
            self.start, self.end, self.strand)

    def _generate_id(self):
        """ Generate an unique ID for the feature. """
        feat_id = "%s_%s_%s_%s_%s" % (self.ftype, self.start, self.end,
                                      self.strand, self.seqid)
        return feat_id

    def _init_from_gffline(self, gffline):
        """ Init the feature from a GFF line. """
        self.seqid, self.source, self.ftype, self.start, self.end, self.score,\
            self.strand, self.phase, self.attributes = read_all_fields(gffline)

    def _init_manually(self, feature_id, kwargs):
        """ Init manually can be useful, in exon creation for instance. """
        self.seqid = kwargs.pop("seqid", "noseqid")
        self.source = kwargs.pop("source", "nosource")
        self.ftype = kwargs.pop("ftype", "noftype")
        self.start = kwargs.pop("start", 0)
        self.end = kwargs.pop("end", 0)
        self.score = kwargs.pop("score", ".")
        self.strand = kwargs.pop("strand", "+")
        self.phase = kwargs.pop("phase", ".")
        self.attributes = kwargs.pop("attributes", "ID=%s" % (feature_id))

    # Properties
    @property
    def ID(self):
        """ Convenient way to access ID tag. """
        return self.tags[ID_TAG]

    @ID.setter
    def ID(self, value):
        """ set_tag function always used to check the entered tag. """
        self.set_tag(ID_TAG, value)

    # Public Methods
    def get_tag(self, tag_name):
        """ Return the value of the tag_name. If not found, return None. """
        try:
            return self.tags[tag_name]
        except KeyError as ke:
            return None

    def set_tag(self, tag_name, tag_value):
        """ Convenient way to set a tag value and check if the tag is legal. """
        if tag_name in self.tags:
            self.tags[tag_name] = tag_value
            self.update_attributes()
        else:
            print("""%s is not a legal tag (avail tags are %s%s%s).
            Use %sadd_tag()%s to add it.""" % (
                tag_name,
                ORIT,
                ",".join(self.tags.keys()),
                NA,
                BLIT,
                NA))

    def add_tag(self, tag_name, tag_value):
        if tag_name not in self.tags:
            self.tags[tag_name] = tag_value
            self.update_attributes()
        else:
            print("%s already in the tags. Use %sset_tag()%s to modify it." % (
                tag_name,
                BLIT,
                NA))

    def rm_tag(self, tag_name):
        try:
            del(self.tags[tag_name])
        except KeyError as ke:
            print("""%s is not a legal tag (avail tags are %s%s%s).
            nothing done.""" % (
                ke,
                ORIT,
                ",".join(self.tags.keys()),
                NA))

    def update_attributes(self):
        """ Update self.attributes taking tag values in self.tags. """
        self.attributes = attributes_string_from_tags(self.tags)

    def gffline(self):
        """ Return the gff line (without \n) corresponding to the feature. """
        # To be up-to-date if self.tag was changed without set/add methods
        self.update_attributes()

        gff_line = GFF_SEP.join([
            self.seqid,
            self.source,
            self.ftype,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.phase,
            self.attributes
        ])
        return gff_line


# class Region():
#     """ Base class for regions.
#
#     The Regions are the different transcripts part like:
#     CDS, UTRs and exons. These regions can contain 1 or more
#     features. If the GFF file provide the features order,
#     the ordering is enforced by these declaration. IF it is not,
#     features will be ranked by theu genomic position, and using
#     the transcript strand, their order will be automatically defined.
#     """
#
#     def __init__(self, features, transcript_strand):
#         """ A region is defined from a list of features objects and
#         the strand of the transcript it belongs to.
#         """
#         self.features = features
#         self.transcript_strand = transcript_strand
#         self._sort_by_rank() # Ensure well ordered features
#
#         # Define start and end (5' => 3')
#         self.start = self.features[0].start
#         self.end = self.features[-1].end
#         self.left = min(self.features).left
#         self.right = max(self.features).right
#
#     # Private methods
#     def __len__(self):
#         """ Return the number of features in the region. """
#         return len(self.features)
#
#     def __iter__(self):
#         """ Make this class iterable on the list of features. """
#         return (f for f in self.features)
#
#     def __repr__(self):
#         """ Convenient visualisation on printing. """
#         return "<%s.%s start=%s, end=%s, tr_strand=%s>" % (
#             __name__,
#             self.__class__.__name__,
#             self.start, self.end, self.transcript_strand)
#
#     def _sort_by_rank(self):
#         """ Sort the features by rank (updates self.features order). """
#         if len(self.features) == 1: # 1 feature, no sorting
#             return
#         elif self.features[0].rank > 0: # ranks were found in GFF
#             self.features = sorted(self.features, key=attrgetter("rank"))
#         elif self.features[0].rank == 0: # ranks were not found, sort by start
#             reverse = False
#             if self.transcript_strand == "-":
#                 reverse = True
#             self.features = sorted(self.features, reverse=reverse)
#
#     # Public methods
#     def length(self):
#         """ Return the length in nucleotide of the Region. """
#         if self.features:
#             recs_lens = [feature.length() for feature in self.features]
#             region_length = sum(recs_lens)
#             return region_length
#         else:
#             return 0
#
#     def split_sequence(self, fapick):
#         """ Return a list with the sequence of each feature of the region. """
#         seqs = [f.sequence(fapick) for f in self.features]
#         return seqs
#
#     def mask(self, fapick):
#         """ Mask the whole region in the fasta pickle. """
#         for feature in self.features:
#             feature.mask(fapick)
#
#     def sequence(self, fapick, n=N_0):
#         """ Return the sequence of the Region from fapick sequences.
#         Sequence is 5' => 3' oriented.
#         """
#         sequence = ""
#         split_seq = self.split_sequence(fapick)
#
#         if not n or abs(n) >= self.length():
#             sequence = "".join(split_seq)
#         elif n > 0:
#             missing = n
#             for feat_seq in split_seq:
#                 if missing >= len(feat_seq) and missing > 0:
#                     sequence += feat_seq
#                     missing -= len(feat_seq)
#                 else:
#                     sequence += feat_seq[:missing]
#                     break
#         else: # n < 0
#             missing = abs(n)
#             for feat_seq in split_seq[::-1]:
#                 if missing >= len(feat_seq) and missing > 0:
#                     sequence = feat_seq + sequence
#                     missing -= len(feat_seq)
#                 else:
#                     sequence = feat_seq[missing*-1 : ] + sequence
#                     break
#         return sequence
#
#     def counts(self, bamfile, shift=12, side=5, kmers=[TOTAL_KEYNAME]):
#         """ Return the total count of reads on the Region. """
#         if self.features:
#             features_counts = [feature.counts(bamfile, kmers=kmers,
#                                              shift=shift, side=side)
#                                for feature in self.features]
#             summed_counts = sum_keys(features_counts)
#             return summed_counts
#         else:
#             return {}
#
#     def pos_cov(self, bamfile, shift=SHIFT, side=SIDE, n=N_0, kmers=KMERS):
#         """ Return the pos_cov dictionary for the region. """
#         # Manually init as StraGeRa pos_cov, need automatic way
#         region_pos_cov = {
#             "positions":[],
#             "coverages":[],
#             "details":[],
#             "bank_size":0
#         }
#         # Keys to be extended while parsing all the generated dicts
#         EXT_KEYS = [
#             "coverages",
#             "positions",
#             "details",
#         ]
#
#         if not n or abs(n) > self.length():
#             for feature in self.features:
#                 pos_cov = feature.pos_cov(bamfile, kmers=kmers,
#                                           shift=shift, side=side)
#                 extend_list_keys(region_pos_cov, pos_cov, EXT_KEYS)
#         elif n > 0:
#             missing = n
#             for feature in self.features:
#                 if missing >= feature.length():
#                     pos_cov = feature.pos_cov(bamfile, kmers=kmers,
#                                               shift=shift, side=side)
#                     missing -= feature.length()
#                     extend_list_keys(region_pos_cov, pos_cov, EXT_KEYS)
#                 else:
#                     pos_cov = feature.pos_cov(bamfile, kmers=kmers, n=missing,
#                                               shift=shift, side=side)
#                     extend_list_keys(region_pos_cov, pos_cov, EXT_KEYS)
#                     break
#         else: # n < 0, reverse list order and insert new lists
#             missing = abs(n)
#             for feature in self.features[::-1]:
#                 if missing >= feature.length():
#                     pos_cov = feature.pos_cov(bamfile, kmers=kmers,
#                                               shift=shift, side=side)
#                     missing -= feature.length()
#                     insert_list_keys(region_pos_cov, pos_cov, EXT_KEYS)
#                 else:
#                     pos_cov = feature.pos_cov(bamfile, kmers=kmers,
#                                               n=missing * -1, # 3' side
#                                               shift=shift, side=side)
#                     insert_list_keys(region_pos_cov, pos_cov, EXT_KEYS)
#                     break
#
#         region_pos_cov["bank_size"] = pos_cov["bank_size"]
#         # region_pos_cov now correctly populated with 5' => 3' order kept
#         return region_pos_cov
#
#     def tupleview(self):
#         """ Return the tupleview. """
#         tupleviews = [feature.tupleview() for feature in self.features]
#         return tupleviews
#
#     def detect_introns(self):
#         """ Return limits (5' => 3') of the introns in the region. """
#         strand_one = STRANDS_MAP[self.transcript_strand]["one"]
#
#         intron = namedtuple("intron", ["start", "end"])
#         introns_limits = []
#         for idx, feature in enumerate(self.features):
#             try:
#                 next_feat = self.features[idx + 1]
#             except IndexError:
#                 # last element of self.features is the current feature
#                 return introns_limits
#             intron_start = feature.end + strand_one
#             intron_end = next_feat.start - strand_one
#             introns_limits.append(intron(intron_start, intron_end))
#         return intron_limits
#
#     # Need plot function here => Only 'mature', it means all the features
#     # will be concatenated
#     # Make the DATAFRAME => then use kmers or total vbar_figure
#     def bokeh_figure(self, shift=SHIFT, side=SIDE, only_total=False):
#         """ """



# class Cds(Region):
#     """ Cds class has several attributes corresponding to biological CDSes.
#     """
#     def __init__(self, features, transcript_strand):
#         """ Same init as the Region class. """
#         super().__init__(features, transcript_strand)
#
#     def __repr__(self):
#         """ Convenient visualisation. """
#         return "<%s.%s start=%s, end=%s, tr_strand=%s>" % (
#             __name__,
#             self.__class__.__name__,
#             self.start, self.end, self.transcript_strand)
#
#     # Public methods
#     def length(self, unit="nt"):
#         """ length can be returned in nucleotide (nt) or codon (codon). """
#         nt_length = super().length()
#         if unit == "nt":
#             return nt_length
#         elif unit == "codon":
#             codon_len = nt_length / 3
#             if not codon_len.is_integer():
#                 print("%sWarning%s: CDS length must be a multiple of 3." % (
#                     ORIT, NA))
#             return int(codon_len)
#         else:
#             print("%s%s%s is not a length unit ('codon' or 'nt')." % (
#                 REIT, unit, NA))
#             return -1 # Not legal length for codons
#
#     def pos_cov(self, bamfile, shift=SHIFT, side=SIDE, n=N_0,
#                 kmers=KMERS, tr=None):
#         """ Overide the Region pos_cov because CDS can be plotted alone. When
#         it happens, the gen_rel_map of the transcript should be updated.
#         """
#         cds_pos_cov = super().pos_cov(bamfile, kmers=kmers,
#                                       n=n, shift=shift, side=side)
#         if tr:
#             tr.gen_rel_map = get_genomic_relative_map(cds_pos_cov["positions"])
#         return cds_pos_cov
#
#     def start_codon(self, fapick):
#         """ Return the sequence of the start codon. """
#         return super().sequence(fapick, n=3)
#
#     def stop_codon(self, fapick):
#         """ Return the sequence of the stop codon. """
#         return super().sequence(fapick, n=-3)
#
#     def allframes_reads(self, bamfile, shift=SHIFT, side=SIDE):
#         """ Return a <dict> with the number and percentage
#         of reads in each phase. No n selection implemented. Can it be useful ?
#         """
#         frames = {0:{}, 1:{}, 2:{}}
#         pos_cov = self.pos_cov(bamfile, shift=shift, side=side)
#         for idx, covs in enumerate(pos_cov["details"]):
#             pos_frame = idx % 3
#             for kmer in covs: # Total is counted as a kmer
#                 try:
#                     frames[pos_frame][kmer] += covs[kmer]
#                 except KeyError:
#                     frames[pos_frame][kmer] = covs[kmer]
#         return frames
#
#     def frames_coverages(self, bamfile, shift=SHIFT, side=SIDE, kmers=KMERS):
#         """ Return a <dict> with a coverage list for each frame. """
#         frames_covs = {0:[], 1:[], 2:[]}
#         pos_cov = self.pos_cov(bamfile, shift=shift, side=side, kmers=kmers)
#         for (idx, pos), covs in zip(enumerate(pos_cov["coverages"]),
#                                     pos_cov["details"]):
#             kmers_cov = [pos_cov["details"][idx][k] for k in kmers]
#             total_cov = sum(kmers_cov)
#             pos_frame = idx % 3
#             frames_covs[pos_frame].append(total_cov)
#         return frames_covs
#
#     def start_around(self, bamfile, meta_start,
#                      shift=SHIFT, side=SIDE, kmers=KMERS,
#                      inner=INNER, outer=OUTER):
#         """ Fill the meta_start with the coverage of the positions
#         in the OUTER-INNER interval. Used for meta genes.
#         """
#         strand_one = STRANDS_MAP[self.transcript_strand]["one"]
#         five_prime = self.start - (outer * strand_one)
#         three_prime = self.start + (inner * strand_one)
#         # Use the first feature seqid.
#         seqid = self.features[0].seqid
#         get_env_reduced_pos(
#             seqid, self.transcript_strand, five_prime, three_prime,
#             bamfile, meta_start, self.start,
#             shift=shift, side=side, kmers=kmers)
#
#     def stop_around(self, bamfile, meta_stop,
#                     shift=SHIFT, side=SIDE,
#                     inner=INNER, outer=OUTER, kmers=KMERS):
#         """ Fill the meta_start with the coverage of the positions
#         in the OUTER-INNER interval. Used for meta genes.
#         """
#         strand_one = STRANDS_MAP[self.transcript_strand]["one"]
#         five_prime = self.end - (inner * strand_one)
#         three_prime = self.end + (outer * strand_one)
#         seqid = self.features[0].seqid
#         get_env_reduced_pos(
#             seqid, self.transcript_strand, five_prime, three_prime,
#             bamfile, meta_stop, self.end,
#             shift=shift, side=side, kmers=kmers)


# class Transcript(Feature):
#     """ Base class for Transcripts.
#
#     Transcripts are defined from their parent (gene) id and name.
#     A transcript is made of regions.
#     """
#
#     def __init__(self, gene_id, gene_name, record=None,
#                  transcript_id=None, gffline=None):
#         """ Constructor from GffRecord (loaders.py) or gffline. """
#         # Transcript attributes
#         self._regions = {}
#         self.gene_name = gene_name
#         self.gen_rel_map = None # Will be init from pos_cov functions
#
#         self.start = None
#         self.end = None
#         if record: # only record initialisation builds the complete hierarchy
#             self._init_from_record(record, gene_id, gene_name,
#                                    transcript_id=transcript_id)
#             # From dictionary of features, init the regions as attributes
#             self._init_regions() # also update self.start and self.end
#
#         elif gffline:
#             super().__init__(gffline=gffline)
#         else:
#             print("At least one of record or gffline initialisation"
#                   " methods has to be selected to initialise a transcript."
#                   " Abort.")
#             raise SystemExit
#
#     def _init_from_record(self, record, gene_id, gene_name, transcript_id=None):
#         """ init from record, used in Gene class.  """
#         # Record can be a 'gene' or a 'transcript'
#         super().__init__(gffline=record.s)
#
#         if transcript_id: # Record is a gene without transcript defined
#             self.set_tag(ID_TAG, transcript_id)
#             # a gene doesn't have parent, declared explicitely here
#             self.add_tag(PARENT_TAG, gene_id)
#
#         # Init the feature objects for childs (CDS, UTRs, ...)
#         for child_id in record.c:
#             child_rec = record.c[child_id]
#             feature = Feature(gffline=child_rec.s)
#             check_utr_label(feature)
#             try: # populate each region with it's features
#                 self._regions[feature.ftype].append(feature)
#             except KeyError:
#                 self._regions[feature.ftype] = [feature]
#
#     # Private methods
#     def __repr__(self):
#         return "<%s.%s tr_id=%s, gene_id=%s, strand=%s>" % (
#             __name__,
#             self.__class__.__name__,
#             self.ID, self.get_tag(PARENT_TAG), self.strand)
#
#     def _init_regions(self):
#         """ Init the regions as attributes. """
#         for region_name in self._regions:
#             features = self._regions[region_name]
#             if region_name == CDS:
#                 region = Cds(features, self.strand)
#             else:
#                 region = Region(features, self.strand)
#             setattr(self, region_name, region)
#         if not hasattr(self, EXON):
#             self._init_exons()
#
#     def _init_exons(self):
#         """ If no exons are annotated in the GFF, they are automatically
#         generated if CDS is defined as attribute. Otherwise, no EXON region
#         available, which can generate troubles for Transcript's methods.
#         """
#         if not hasattr(self, CDS):
#             print("No CDS detected for the transcript: %s,"
#                   " exons can't be initialized. Skip." % (
#                       self.ID))
#             return
#
#         exons_tuples = []
#         for reg in ORDERED_REGIONS:
#             try:
#                 region = getattr(self, reg)
#             except AttributeError:
#                 continue
#             for feature in region: # Features are already sorted 5'->3'
#                 exons_tuples.append(feature.tupleview(strand="symbol"))
#
#         if exons_tuples:
#             exons_features = []
#             merge_mixt_exons(exons_tuples) # Update exons_tuples
#             for idx, exon in enumerate(exons_tuples):
#                 exon_number = idx + 1 # 1 based index for exon number
#                 exon_feat = self._create_exon(exon, exon_number)
#                 exons_features.append(exon_feat)
#         else:
#             print("Empty list of exons, should not append, at least CDS")
#             return
#
#         region = Region(exons_features, self.strand) # EXON region init
#         setattr(self, EXON, region)
#         # Adjut the transcript start and end to the exons
#         self.start = region.start
#         self.end = region.end
#
#     def _create_exon(self, exon_tuple, exon_number):
#         """ From an exon_tuple, create an EXON Feature. """
#         exon_id = "%s:%s:%s" % (EXON, self.ID, exon_number)
#         exon_attributes = "{id_tag}={id_val};{parent_tag}={parent_val}"\
#                           ";{rank_tag}={rank_val}".format(
#                               id_tag=ID_TAG,
#                               parent_tag=PARENT_TAG,
#                               rank_tag=RANK_TAG,
#                               id_val=exon_id,
#                               parent_val=self.ID,
#                               rank_val=exon_number
#                           )
#         exon_feat = Feature(
#             exon_id,
#             start=exon_tuple.start,
#             end=exon_tuple.end,
#             strand=exon_tuple.strand,
#             ftype=EXON,
#             seqid=self.seqid,
#             attributes=exon_attributes
#         )
#         return exon_feat
#
#     # Public methods
#     def introns_limits(self):
#         """ Return a list of tuples with the introns limits (5' => 3'). """
#         introns_limits = self.exon.detect_introns()
#         return introns_limits
#
#     def introns_sequences(self, fapick):
#         """ Return the sequences of each introns.
#         The order returned is 5' => 3' in the transcript.
#         """
#         i_seqs = []
#         i_lims = self.introns_limits()
#         for intron in i_lims:
#             i_stragera = StraGeRa(intron.start, intron.end,
#                                   self.strand, self.seqid)
#             i_seqs.append(i_stragera.sequence(fapick))
#         return i_seqs
#
#     def counts(self, bamfile, shift=SHIFT, side=SIDE, kmers=KMERS):
#         """ Return a dictionary with the count in each features for each
#         kmer in 'kmers'.
#         """
#         counts = {}
#         for reg_name in ORDERED_REGIONS:
#             if hasattr(self, reg_name):
#                 region = getattr(self, reg_name)
#                 counts[reg_name] = region.counts(bamfile, shift=shift,
#                                                  side=side, kmers=kmers)
#             else:
#                 counts[reg_name] = {}
#         return counts
#
#     def mature_sequence(self, fapick, n=N_0):
#         """ Return exon region sequence. """
#         mature_seq = self.exon.sequence(fapick, n=n)
#         # Need to add exception if some transcripts havn't got Exons...
#         return mature_seq
#
#     def pre_sequence(self, fapick, n=N_0):
#         """ Return the pre-mature sequence of the transcript.
#         Results will be consistents only for linear genes.
#         """
#         tr_seq = self.sequence(fapick, n=n)
#         return tr_seq
#
#     def mature_pos_cov(self, bamfile, shift=SHIFT, side=SIDE,
#                        kmers=KMERS, n=N_0):
#         """ Return the mature mRNA positions and coverage from EXON region. """
#         tr_pos_cov = self.exon.pos_cov(bamfile, shift=shift, side=side,
#                                        kmers=kmers, n=n)
#         # Init correspondances dictionary for relative and genomic positions
#         self.gen_rel_map = get_genomic_relative_map(tr_pos_cov["positions"])
#         return tr_pos_cov
#
#     def pre_pos_cov(self, bamfile, shift=SHIFT, side=SIDE,
#                     kmers=KMERS, n=N_0):
#         """ Return the pre-mRNA (with intron) positions and coverages. """
#         tr_pos_cov = self.pos_cov(bamfile, shift=shift, side=side,
#                                   n=n, kmers=kmers)
#         self.gen_rel_map = get_genomic_relative_map(tr_pos_cov["positions"])
#         return tr_pos_cov
#
#     # Need a method to have the mature length...!
#     def mature_length(self):
#         """ """




# class Gene():
#     """ Base class for genes
#     Intentionnaly, Gene is not considered as a Feature. This for
#     several reasons. The strongest is that a Gene limits can be defined
#     by one or more transcripts (isoforms). Moreover, EVERY gene has at
#     least one transcript. So Gene class is based on a list of transcripts
#     to compute genome arithmetics. The Gene start and end are defined by
#     the larger limits, even if two transcripts have to be combined.
#     """
#     def __init__(self, name, record, transcript_keyword=TRANSCRIPT):
#         """ initialised from a GffRecord object (loader.py).
#         Transcript_keyword is correspond to the transcript identifier in
#         the third column of GFF file.
#         """
#         self.ID = record.i
#         self.seqid = record.k
#         self.strand = record.sd
#         # add the source ? check time + memory
#         self.Name = name
#         self.tags = get_tags_and_values(record.s.split(GFF_SEP)[ATTRIBUTES])
#         self.transcripts = []
#
#         self._init_hierarchy(record, transcript_keyword)
#
#     def _init_hierarchy(self, record, transcript_keyword):
#         """ Init hierarchy with childs of the record
#         (transcripts, regions, features).
#         """
#         childs_types = [record.c[recid].t for recid in record.c]
#         if transcript_keyword in childs_types:
#             # init each transcript
#             for trid in record.c:
#                 tr_rec = record.c[trid]
#                 transcript = Transcript(record.i, self.Name, record=tr_rec)
#                 self.transcripts.append(transcript)
#         else: # No annotated transcript. Create one with a new trid
#             trid = "%s_%s" % (transcript_keyword, self.ID)
#             transcript = Transcript(self.ID,
#                                     self.Name,
#                                     record=record,
#                                     transcript_id=trid)
#             self.transcripts.append(transcript)
#
#     def __repr__(self):
#         return "<%s.%s ID=%s, Name=%s, seqid=%s>" % (
#             __name__,
#             self.__class__.__name__,
#             self.ID, self.Name, self.seqid)
#
#     # Properties
#     @property
#     def start(self):
#         """ The start is the first start of all transcripts. """
#         larger_5prime = 0
#         for transcript in self.transcripts:
#             if not larger_5prime:
#                 larger_5prime = transcript.start
#             else:
#                 if self.strand == "+":
#                     if transcript.start < larger_5prime:
#                         larger_5prime = transcript.start
#                 else: # strand == "-"
#                     if transcript.start > larger_5prime:
#                         larger_5prime = transcript.start
#         return larger_5prime
#
#     @property
#     def end(self):
#         """ The end is the last end of all transcripts. """
#         larger_3prime = 0
#         for transcript in self.transcripts:
#             if not larger_3prime:
#                 larger_3prime = transcript.end
#             else:
#                 if self.strand == "+":
#                     if transcript.end > larger_3prime:
#                         larger_3prime = transcript.end
#                 else:
#                     if transcript.end < larger_3prime:
#                         larger_3prime = transcript.end
#         return larger_3prime
#
#     # Public method
#     def get_tag(self, tag_name):
#         """ Return the value associated with the tag_name. """
#         try:
#             return self.tags[tag_name]
#         except KeyError as ke:
#             print("%s not in the tags of %s gene" % (
#                 ke, self.ID))
#             return None
#
#     def get_transcript(self, transcript_id):
#         """ Return the Transcript object with transcript_id.
#         It not found, return None.
#         """
#         for tr in self.transcripts:
#             if tr.ID == transcript_id:
#                 return tr
#         print("%s not found in %s transcripts list." % (
#             transcript_id, self.ID))
#         return None
#
#     def select_transcripts(self, transcripts_to_select):
#         """ Return a generator with the selected transcripts (Transcript). """
#         # Need a new version, where a dictionary of geneid:[tr_to_select]
#         # is provided, and for all geneid, extract the wanted transcripts.
#         found_idxs = []
#         for (idx, tr) in enumerate(self.transcripts):
#             if tr.ID in transcripts_to_select:
#                 found_idxs.append(idx)
#                 yield tr
#             else:
#                 continue
#         if len(found_idxs) < len(transcripts_to_select):
#             not_found = [transcripts_to_select[idx]
#                          for (idx, _) in enumerate(transcripts_to_select)
#                          if idx not in found_idxs]
#             print("Transcripts not found: %s" % (", ".join(not_found)))
#
#     def appris_transcripts(self):
#         """ Return a generator to iterate only on transcripts (Transcript)
#         with 'appris_principal' tag. Usually it allows a significant
#         filter of spurious transcripts.
#         """
#         for tr in self.transcripts:
#             appris_tag = tr.get_tag("tag")
#             if appris_tag and APPRIS_P_LEV in appris_tag:
#                 yield tr
#             else:
#                 pass
#
#     def gffline(self):
#         """ Return the gffline (without ending "\n") of the gene. """
#         attributes = attributes_string_from_tags(self.tags)
#         source = "BYHAND" # Only source have to be checked in order to keep
#         # the orignial source...! (from record perhaps ?)
#         score = "."
#         phase = "."
#         start = min(self.start, self.end)
#         end = max(self.start, self.end)
#         gff_line = GFF_SEP.join([
#             self.seqid,
#             source,
#             GENE,
#             str(start),
#             str(end),
#             score,
#             self.strand,
#             phase,
#             attributes
#         ])
#         return gff_line
#
#     def get_transcript_figure(
#             self, bamfile, transcript_id=None, shift=SHIFT, side=SIDE,
#             kmers=KMERS, state="mature", only_total=False, n=N_0,
#             colors=None, glyphs=["annot"], fapick=None,
#             scaling=False, rpm=False, split_signal=False, total_color="black"):
#         """ Plot the selected transcript of the gene. If no transcript
#         provided, the first transcript of the list is plotted.
#         """
#         if not self.transcripts:
#             print("No transcript defined for %s" % (self.ID))
#             return
#         if not transcript_id:
#             tr_to_plot = self.transcripts[0]
#         else:
#             tr_to_plot = self.get_transcript(transcript_id)
#             if not tr_to_plot:
#                 print("%s not in %s transcripts list." % (
#                     transcript_id, self.ID))
#                 raise SystemExit
#
#         if ("nt" in glyphs or "codon" in glyphs) and not fapick:
#             print("Fasta pickle (fapick=) is required to display"
#                   " codon or nucleotide glyphs")
#             return
#
#         # pos_cov and sequence are generated accordingly to the transcript state
#         sequence = ""
#         if state == "mature":
#             pos_cov = tr_to_plot.mature_pos_cov(bamfile, shift=shift, side=side,
#                                                 kmers=kmers, n=n)
#             if fapick:
#                 sequence = tr_to_plot.mature_sequence(fapick, n=n)
#         elif state == "pre":
#             pos_cov = tr_to_plot.pre_pos_cov(bamfile, shift=shift, side=side,
#                                              kmers=kmers, n=n)
#             if fapick:
#                 sequence = tr_to_plot.pre_sequence(fapick, n=n)
#         elif state == "cds":
#             pos_cov = tr_to_plot.CDS.pos_cov(bamfile, shift=shift, side=side,
#                                              kmers=kmers, n=n, tr=tr_to_plot)
#             if fapick:
#                 sequence = tr_to_plot.CDS.sequence(fapick, n=n)
#         else:
#             print("%s state is not legal. 'mature', 'pre' or 'cds' "
#                   "excepted." % (state))
#             return
#         # Scalings
#         y_label = "Ribosome occupancy"
#         if rpm:
#             pos_cov["coverages"] = rpm_normalisation(pos_cov["coverages"],
#                                                      pos_cov["bank_size"])
#             y_label += " (RPM)"
#         if scaling:
#             pos_cov["coverages"] = scale_huge_peaks(pos_cov["coverages"])
#             y_label += " (scaled)"
#
#         title = "%s - %s" % (self.ID, tr_to_plot.ID)
#         if not split_signal:
#             tr_figure = figure_from_poscov(pos_cov,
#                                            tr_to_plot,
#                                            glyphs=glyphs,
#                                            sequence=sequence,
#                                            title=title,
#                                            only_total=only_total,
#                                            color=total_color)
#             tr_figure.yaxis.axis_label = y_label
#             tr_figure.xaxis.axis_label = "Genomic positions (%s)" % (self.seqid)
#         else:
#             tr_figure = figure_splitted_signal(pos_cov,
#                                                tr_to_plot,
#                                                glyphs=glyphs,
#                                                sequence=sequence,
#                                                title=title,
#                                                only_total=only_total,
#                                                color=["blue", "green", "red"])
#             # Link plots + ranges ? or ranges already done?
#         return tr_figure
