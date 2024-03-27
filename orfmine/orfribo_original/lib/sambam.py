#!/usr/bin/env python3

""" Module to handle SAM/BAM files.

The notion of 'alignment strand' is introduced here.
For RIBO-SEQ data, the flags corresponding to the
alignment strands are 0 ('+') and 16 ('-').
"""

# Stdlib
import os

# Third
import pysam

# Local
from .constants import *

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"

def get_bank_size(bamfile):
    """ Return the total number of reads mapped in the bamfile. """
    bam = init_pysam_bam(bamfile)
    return bam.mapped


def get_stranded_region_coverage_reduced(
        chromosome, strand, start_1_based, end_1_based, bamfile,
        shift=SHIFT, side=SIDE, kmers=KMERS):
    """ Return a dictionary with the coverage of each position in the given
    region. The reads are reduced to only one position, defined by the 'shift'
    and the 'side'. 'kmers' is used to selected specific sizes of reads.
    """
    not_found_refs = {}
    chromosome = chromosome.strip()
    # start MUST be LESS that end.
    start = min(start_1_based, end_1_based) - 1
    end = max(start_1_based, end_1_based) - 1
    try:
        ali_strand = STRANDS_MAP[strand]["ali"]
    except KeyError as ke:
        print("%s is not a legal strand ('+' or '-' expected)" % (
            ke))
        return {}
    reduced_index = get_shift_idx(shift, side, ali_strand)
    # Only stranded reads will be kept
    pos_cov = {
         ali_strand:{}
    }
    # +1 to be 1based index, and end + 1 once again to inlude the end
    for pos in range(start+1, end + 2):
        pos_cov[ali_strand][pos] = {TOTAL_KEYNAME:0}
        # pos is 1 based index

    bam = init_pysam_bam(bamfile)
    try:
        for read in bam.fetch(chromosome, start, end + 1):
            reduced_pos = read.positions[reduced_index]
            try: # Each concerned pos is already init
                pos_cov[read.flag][reduced_pos + 1][TOTAL_KEYNAME] += 1
                try:
                    pos_cov[read.flag][reduced_pos + 1][read.rlen] += 1
                except KeyError as ke: #New key for the current kmer size
                    pos_cov[read.flag][reduced_pos + 1][read.rlen] = 1
            except KeyError as ke:
                # Bad strand or reduced_pos not in [start, end] interval
                continue
    except ValueError: #Bad chromosome
        try:
            not_found_refs[chromosome] += 1
        except KeyError:
            not_found_refs[chromosome] = 1

    if not_found_refs:
        print("References not found in %s:" % (bamfile))
        for k in not_found_refs:
            print("%s: %s reads" % (k, not_found_refs[k]))

    if kmers != [TOTAL_KEYNAME]: #Kmer selection
        return select_kmers(kmers, pos_cov[ali_strand])
    else:
        return pos_cov[ali_strand]


def select_kmers(kmers, pos_cov):
    """ Select the kmers in `kmers` for each positions
    in pos_cov. Return a new dictionary with TOTAL adjusted.
    """
    selected_kmers = {}
    for pos in pos_cov:
        selected_kmers[pos] = {TOTAL_KEYNAME:0}
        current_pos_cov = pos_cov[pos]
        for kmer in kmers:
            if kmer in pos_cov[pos]:
                kmer_cov = current_pos_cov[kmer]
                selected_kmers[pos][TOTAL_KEYNAME] += kmer_cov
                if kmer != TOTAL_KEYNAME:
                    try:
                        selected_kmers[pos][kmer] += kmer_cov
                    except KeyError as ke:
                        selected_kmers[pos][kmer] = kmer_cov

    return selected_kmers


def get_stranded_region_counts(chromosome, strand,
                     start_1_based, end_1_based, bamfile,
                     shift=12, side=5, kmers=[TOTAL_KEYNAME]):
    """ Return a dictionary with the count of each kmers
    between start and end. Only stranded reads are kept.
    """
    start = min(start_1_based, end_1_based) - 1
    end = max(start_1_based, end_1_based) - 1
    ali_strand = STRANDS_MAP[strand]["ali"]
    reduced_index = get_shift_idx(shift, side, ali_strand)

    if kmers == [TOTAL_KEYNAME]:
        counts = count_all_kmers(chromosome,
                                 ali_strand,
                                 start,
                                 end,
                                 bamfile,
                                 reduced_index)
    else:
        counts = count_selected_kmers(chromosome,
                                      ali_strand,
                                      start,
                                      end,
                                      bamfile,
                                      reduced_index,
                                      kmers)
    return counts[ali_strand]

def count_selected_kmers(chromosome, ali_strand,
                         start, end, bamfile, reduced_index, kmers):
    """ Populate the kmers_counts with the selected kmers and the total. """
    kmers_counts = {
        ali_strand:{
            TOTAL_KEYNAME:0
        }
    }
    not_found_refs = []
    for kmer in kmers:
        kmers_counts[ali_strand][kmer] = 0

    bam = init_pysam_bam(bamfile)
    try:
        for read in bam.fetch(chromosome, start, end + 1): # +1 as end excluded
            reduced_pos = read.positions[reduced_index]
            # start and end already have the -1 shift applied
            if reduced_pos >= start and reduced_pos <= end:
                try:
                    kmers_counts[read.flag][read.rlen] += 1
                    kmers_counts[read.flag][TOTAL_KEYNAME] += 1
                except KeyError as ke:  # Bad Kmer or bad strand, skip read
                    continue
            else: # read was fetched as it overlaps start-end region, not needed
                continue
    except ValueError as ve:
        # Invalid reference (chromosome)
        if ve not in not_found_refs:
            not_found_refs.append(ve)
    return kmers_counts

def count_all_kmers(chromosome, ali_strand, start, end, bamfile, reduced_index):
    """ Populate all_kmers_counts with all the Kmers found in the bamfile. """
    all_kmers_counts = {
        ali_strand:{
            TOTAL_KEYNAME:0
        }
    }
    not_found_refs = []
    bam = init_pysam_bam(bamfile)

    try:
        for read in bam.fetch(chromosome, start, end + 1):
            reduced_pos = read.positions[reduced_index]
            if reduced_pos >= start and reduced_pos <= end:
                try:
                    all_kmers_counts[read.flag][TOTAL_KEYNAME] += 1
                except KeyError as ke: # Bad strand
                    continue
                try:
                    all_kmers_counts[read.flag][read.rlen] += 1
                except KeyError as ke:  # Kmer was not initialised
                    all_kmers_counts[read.flag][read.rlen] = 1
            else: # Idem as above, read fetched but not selected here
                continue
    except ValueError as ve:
        # Invalid reference (chromosome)
        if ve not in not_found_refs:
            not_found_refs.append(ve)
    return all_kmers_counts


def get_shift_idx(shift, side, ali_strand):
    """ Return the reduced position index with shift applied.
    shift = 0 means 5' or 3' most positions. Side can be 3 or 5,
    and ali_strand 0 or 16.
    """
    if shift == 0:
        if side == 5:
            if ali_strand == 0:
                return 0
            else:
                return -1
        else:
            if ali_strand == 0:
                return -1
            else:
                return 0
    elif shift > 0:
        if side == 5:
            if ali_strand == 0:
                return shift
            else:
                return (shift * -1) - 1
        else:
            if ali_strand == 0:
                return (shift * -1) - 1
            else:
                return shift
    else:
        print("Shift must be >= 0")
        return 0


def init_pysam_bam(bamfile_path, auto_sort=False):
    """ Check if the bamfile has an index. If not, the file is indexed.
    If it can be indexed, it will be sorted.
    Returns a pysam.AlignmentFile.
    """
    ali_file = pysam.AlignmentFile(bamfile_path, "rb")
    output_bam = bamfile_path.split(".")[0] + "_sorted_bystrali.bam"
    try:
        ali_file.check_index()
        return ali_file
    except ValueError: # File not indexed, try to do it
        try:
            print("Indexing %s" % (bamfile_path))
            pysam.index(bamfile_path)
            return ali_file
        except pysam.SamtoolsError as se: # File is not sorted
            if not auto_sort and not os.path.isfile(output_bam):
                print("%s is not sorted, it can't be indexed." % (bamfile_path))
                return ali_file
            elif auto_sort and os.path.isfile(output_bam):
                print("%s is already sorted: %s" % (bamfile_path, output_bam))
                return pysam.AlignmentFile(output_bam, "rb")
            elif auto_sort:
                try:
                    print("Sorting %s" % (bamfile_path))
                    pysam.sort("-o", output_bam, bamfile_path)
                    print("Indexing %s" % (output_bam))
                    pysam.index(output_bam)
                    return pysam.AlignmentFile(output_bam, "rb")
                except pysam.SamtoolsError as se:
                    print("Error: %s, file can't be sorted." % (se))
                    return ali_file


def get_env_reduced_pos(reference, strand, five_prime_1b, three_prime_1b,
                        bamfile, meta_codon, ref_pos_1b,
                        shift=SHIFT, side=SIDE, kmers=KMERS):
    """ Automatically fill the meta_codon with the relative position between
    five_prime1b and three_prime_1b.
    """
    not_found_refs = []

    fetch_start = min(five_prime_1b, three_prime_1b) - 1
    fetch_stop = max(five_prime_1b, three_prime_1b) - 1
    ref_pos = ref_pos_1b - 1
    ali_strand = STRANDS_MAP[strand]["ali"]
    one_strand = STRANDS_MAP[strand]["one"]
    reduced_index = get_shift_idx(shift, side, ali_strand)
    bam = init_pysam_bam(bamfile)
    try_dic = {ali_strand:0}
    #print(fetch_start, fetch_stop)
    try:
        for read in bam.fetch(reference, fetch_start, fetch_stop + 1):
            #print(read.flag, read.positions)
            try:
                try_dic[read.flag] += 1
            except KeyError:
                continue
            #if read.flag != ali_strand:
            #    continue
            reduced_pos = read.positions[reduced_index]
            relative_pos = (reduced_pos - ref_pos) * one_strand
            try:
                meta_codon[relative_pos][read.rlen] += 1
            except KeyError: # Not in meta range or bad kmer
                continue
            meta_codon[relative_pos][TOTAL_KEYNAME] += 1
    except ValueError as ve:
        if ve not in not_found_refs:
            not_found_refs.append(ve)
