#!/usr/bin/env python3

""" Contains functions to parse GFF3 formatted file.

Fields are tab separated, and attributes ';' separated.
These constants are imported from constants module.
"""

# Stdlib
import re

# Local
from .constants import \
    GFF_SEP, GFF_COM, TAG_SEP, VAL_SEP, \
    PARENT_TAG, ID_TAG, UTR, \
    SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTES

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"

# Constants
REGEX_PARENT = re.compile("%s=([^;]+);?" % (PARENT_TAG))
REGEX_ID = re.compile("%s=([^;]+);?" % (ID_TAG))


def read_all_fields(gffline):
    """ Return the value of each field from gffline """
    split_line = gffline.strip().split(GFF_SEP)
    fields = (field for field in split_line)
    return tuple(fields)

def get_tags_and_values(attributes):
    """ Return all the tag and their value in attributes. """
    tags = {}
    tags_values = (
        (tagval.split(VAL_SEP)[0], tagval.split(VAL_SEP)[1])
        for tagval in attributes.split(TAG_SEP))
    for tagval in tags_values:
        tag_name = tagval[0]
        tag_value = tagval[1]
        tags[tag_name] = tag_value
    return tags

def gff_tag(tag_name, gffline):
    """ Return the value of the given tag from a gffline.
    None is returned if tag_name is not found.
    """
    regex = re.compile("%s=([^;]+);?" % (tag_name))
    res = regex.search(gffline)
    try:
        return res.group(1)
    except AttributeError:
        return None

def attributes_string_from_tags(tags):
    """ From a tags dictionary of tag_name: tag_values,
    return a string following GFF3 attributes field specifications.
    """
    attributes = []
    for (tag_name, tag_value) in tags.items():
        tag_str = "%s=%s" % (tag_name, str(tag_value))
        if tag_name == ID_TAG:
            attributes.insert(0, tag_str)
        elif tag_name == PARENT_TAG:
            attributes.insert(1, tag_str)
        else:
            attributes.append(tag_str)
    return TAG_SEP.join(attributes)

# def gff_field(field, gffline): #could exclude
#     """ Return the field value from gffline already strip(). """
#     split_line = gffline.split(GFF_SEP)
#     if field == "seqid":
#         return split_line[SEQID]
#     elif field == "strand":
#         return split_line[STRAND]
#     elif field == "type":
#         return split_line[TYPE]
#     elif field == "start":
#         return split_line[START]
#     elif field == "end":
#         return split_line[END]
#     elif field == "attributes":
#         return split_line[ATTRIBUTES]
#     elif field == "source":
#         return split_line[SOURCE]
#     elif field == "phase":
#         return split_line[PHASE]
#     elif field == "score":
#         return split_line[SCORE]
#     else:
#         return ""


def check_rank_in_id(feature_id):
    """ Check if the rank is in ID tag (feature_type:transcript_id:rank).
    If yes, return this number. 0 otherwise.
    """
    splitted_id = feature_id.split(":")
    try:
        last_element = int(splitted_id[-1])
    except ValueError:
        last_element = None

    if len(splitted_id) > 1 and last_element:
        return last_element
    else:
        return 0

# def check_utr_label(feature):#could exclude
#     """ If feature.ftype == UTR, UTR5 or UTR3 have to be defined.
#     Attempt to get it from feature.ID. If not possible, UTR is kept
#     but plots will not contains UTRs.
#     Need to implement a function to detect UTR5 or UTR3.
#     """
#     if feature.ftype != UTR:
#         return
#     else:
#         utr_label = feature.ID.split(":")[0]
#         if UTR in utr_label:
#             feature.ftype = utr_label
#             return
#         else:
#             return
