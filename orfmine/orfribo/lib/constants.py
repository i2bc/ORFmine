""" Constants of the library. """

# Stdlib
from collections import namedtuple

## GFF ##
# Separator
GFF_SEP = "\t"
GFF_COM = "#"
TAG_SEP = ";"
VAL_SEP = "="
# Tags
PARENT_TAG = "Parent"
ID_TAG = "ID"
NAME_TAG = "Name"
RANK_TAG = "exon_number"
# Fields
SEQID = 0
SOURCE = 1
TYPE = 2
START = 3
END = 4
SCORE = 5
STRAND = 6
PHASE = 7
ATTRIBUTES = 8
# Types
GENE = "gene"
TRANSCRIPT = "transcript"
CDS = "CDS"
EXON = "exon"
UTR5 = "UTR5"
UTR3 = "UTR3"
UTR = "UTR"
# APPRIS
APPRIS_P_ALL = "appris_principal" # level of appris_p doesnt matter
APPRIS_LEVEL = 1
APPRIS_P_LEV = "%s_%s" % (APPRIS_P_ALL, str(APPRIS_LEVEL))
#Region order
ORDERED_REGIONS = [UTR5, CDS, UTR3]

## Alignments ##
# strands map
STRANDS_MAP = {
    "+":{
        "one":1,
        "ali":0,
        "ali_rev":16
    },
    "-":{
        "one":-1,
        "ali":16,
        "ali_rev":0
    },
    # "." handle the features like chromosome.
    # It can be usefull to work with them
    ".":{
        "one":1
    }
}
# other
TOTAL_KEYNAME = "Total"
MASK_BASE = "N"
SHIFT = 12
SIDE = 5
N_0 = 0
KMERS = [TOTAL_KEYNAME]
INNER = 100
OUTER = 20
MIN_COORDS = namedtuple("tupleview", ["start", "end", "strand"])
## Colors ##
NA="\033[0m"
ORIT="\033[3;33m"
REIT="\033[3;31m"
GRIT="\033[3;32m"
BLIT="\033[3;34m"
BLITUN="\033[4;34m"
## Frames ##
# First keys are the shift of phases
FRAMES_CORREL = {
    0:{
        0:0,
        1:1,
        2:2
    },
    1:{
        0:2,
        1:0,
        2:1
    },
    2:{
        0:1,
        1:2,
        2:0
    }
}
