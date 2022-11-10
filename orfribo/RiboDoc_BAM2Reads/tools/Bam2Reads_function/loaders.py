#!/usr/bin/env python3

""" Module containing several loaders for ones of the most common file formats
in bioinformatics: GFF, FASTA, BED...
"""

# Local
from constants import *
from biobjects import Feature
from noncoding_objects import Igorf
from gffutils import gff_tag
from gffutils import REGEX_ID, REGEX_PARENT

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"


class Gff():
    """ Handle GFF3 loading and generate Genes objects. """
    def __init__(self, gff_file, name_tag=None, transcript_keyword=TRANSCRIPT,
                 all_as_high=False):
        """ Constructor:
        If 'name_tag' not provided, the name is picked from 'Name' tag.
        Otherwise, the name_tag is used to get gene name.
        'transcript_keyword' is also important to fit the GFF3 hierarchy.
        """
        self.gff_file = gff_file
        self.name_tag = name_tag
        self.transcript_keyword = transcript_keyword

        # Privates, used from generators
        self._highests = {}   # highest features: no parent found
        self._records = {}
        self._idnames = {}
        self._namesid = {}

        # Load info
        self.load()
        # if not all_as_high: #could exclude
        #     self.build_hierarchy()

        self.build_idname_map()

    def load(self):
        """ Load the gfflines in GffRecord. """
        no_id = 0
        print("Loading %s" % (self.gff_file))
        with open(self.gff_file, 'r') as f:
            for line in f:
                if not line.startswith(GFF_COM) and not line.isspace():
                    record = GffRecord(line.strip())
                    try:
                        self._records[record.i] = record
                    except AttributeError: # ID is missing, create one
                        no_id += 1
                        record.i = "NO_ID_%s" % (no_id)
                    finally:
                        self._records[record.i] = record
        if no_id:
            print("%s lines without 'ID=' tag (random IDs generated)" % (no_id))
        print("%s features loaded" % (len(self._records)))

    # def build_hierarchy(self): #could exclude
    #     """ Connect features with Parent tag. """
    #     print("Building hierarchy")
    #     for rec_id in self._records:
    #         rec = self._records[rec_id]
    #         try:
    #             parent = self._records[rec.p] # Parent exists
    #             parent.c[rec_id] = rec        # Add rec to its parent's childs
    #         except KeyError:
    #             self._highests[rec_id] = rec
    #         finally:
    #             parent = None   # Reset needed in case of 'except'

    def build_idname_map(self):
        """ Build 2 dictionaries: id:names and names:id. """
        for feat_id in self._highests:
            missing = None
            feat_rec = self._highests[feat_id]
            if self.name_tag:
                feat_name = gff_tag(self.name_tag, feat_rec.s)
                if not feat_name:
                    missing = self.name_tag
            else:               # Default, NAME_TAG = 'Name'
                feat_name = gff_tag(NAME_TAG, feat_rec.s)
                if not feat_name:
                    missing = NAME_TAG
            if missing:
                feat_name = "%s_tag_missing" % (missing)
            # finally
            self._idnames[feat_id] = feat_name
            self._namesid[feat_name] = feat_id

    def all_features(self, cast_into="Feature"): #could exclude CDS/feature
        """ Return a generator to iterate over all highest features.
        Highest features can be casted in several biobjects, depending
        on the needs. Igorf added for A.Lopes project.
        """
        for feat_id in self._records:
            feat_rec = self._records[feat_id]
            feature = Feature(gffline=feat_rec.s)
            if cast_into == "Igorf":
                yield Igorf(feat_rec.s)
            # elif cast_into == "Feature":
            #     yield feature
            # elif cast_into == "Cds":
            #     cds = Cds([feature], feat_rec.sd)
            #     yield cds
            else:
                print("highest features can be casted to Feature,"
                      " Cds and Igorf only.")
                # Transcrit and gene can be implemented, but is that usefull?
                #as the feat_rec is only one genome portion.
                return

    # def genes(self): #could exclude
    #     """ Return a generator with all genes (Gene) in the GFF. """
    #     for feat_id in self._highests:
    #         feat_rec = self._highests[feat_id]
    #         if feat_rec.t == GENE:
    #             gene_name = self._idnames[feat_id]
    #             gene = Gene(gene_name,
    #                         feat_rec,
    #                         transcript_keyword=self.transcript_keyword)
    #             yield gene


    # def select_gene_names_id(self, selected_genes): #could exclude
    #     """ Return the Gene object for the selected IDs/Names. """
    #     for g in selected_genes:
    #         if g in self._namesid:
    #             gene_id = self._namesid[g]
    #             gene_name = g
    #         elif g in self._idnames:
    #             gene_id = g
    #             gene_name = self._idnames[g]
    #         else:
    #             print("%s%s%s gene id/name not found" % (REIT, g, NA))
    #             continue
    #         record = self._highests[gene_id]
    #         gene = Gene(gene_name,
    #                     record,
    #                     transcript_keyword=self.transcript_keyword)
    #         yield gene
    #
    # def select_gene_tags(self, **kwargs): #could exclude
    #     """ Select genes depending on the tags and their values. """
    #     for feat_id in self._highests:
    #         keep_me = True
    #         feat_rec = self._highests[feat_id]
    #         if feat_rec.t == GENE:
    #             gene = Gene(self._idnames[feat_id],
    #                         feat_rec,
    #                         transcript_keyword=self.transcript_keyword)
    #             for (tag_name, tag_value) in kwargs.items():
    #                 try:
    #                     # Need a version where, if AT LEAST 1 tag is ok,
    #                     # the gene is kept ?
    #                     if gene.tags[tag_name] != tag_value:
    #                         keep_me = False
    #                         break
    #                 except KeyError:
    #                     continue # Inexistant tag
    #             if not keep_me:
    #                 continue
    #             else:
    #                 yield gene
    #
    #
    # def select_feature(self, feature_identifiers): #could exclude
    #     """ Select the features of self._records only if they match
    #     the feature_identifiers list.
    #     """
    #     selected_features = []
    #     # Return Features ?...? StraGeRa ?
    #     # Or implement something to load all feature as high,
    #     # it means all feature as StraGeRa ? => Like this, all
    #     # genomic arithmetic is available.
    #
    # def features_summary(self): #could exclude
    #     """ Return a summary of the features present in the GFF. """
    #     # Need to be implement, pretty print of a summary
    #     # of each features type + feature level


class GffRecord():
    """ A minimal object to save memory and time while loading. """
    __slots__ = ["i", "p", "c", "s", "t", "sd", "k"]

    def __init__(self, s):
        """ minimal view of gff line. """
        self.s = s              # gff line.strip()
        self.c = {}             # childs
        self.p = self.get_parent() # parent
        self.t = self.get_type()   # type
        self.sd = self.get_strand() # strand
        self.k = self.get_seqid()   # seqid
        try:
            self.i = self.get_id() # ID
        except AttributeError:
            pass # No self.i declared => attributeError raised in Gff loader

    def get_parent(self):
        """ Return the Parent of the record. """
        attributes = self.s.split(GFF_SEP)[-1]
        res = REGEX_PARENT.search(attributes)
        try:
            return res.group(1)
        except AttributeError:
            return None

    def get_id(self):
        attributes = self.s.split(GFF_SEP)[-1]
        res = REGEX_ID.search(attributes)
        return res.group(1)

    # Short methods to access GFF fields to avoid testing using gff_field method
    def get_type(self):
        return self.s.split(GFF_SEP)[TYPE]

    def get_seqid(self):
        return self.s.split(GFF_SEP)[SEQID]

    def get_strand(self):
        return self.s.split(GFF_SEP)[STRAND]
