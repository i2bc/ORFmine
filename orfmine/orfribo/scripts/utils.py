#!/usr/bin/env python3

""" Some util functions """

# Third
import pandas as pd
from numpy import mean, median

# Local
from constants import MIN_COORDS, FRAMES_CORREL, TOTAL_KEYNAME




def find_stretches(coverage, threshold=1, min_len=10):
    """ Detect stretches of values < threshold in coverage (list).
    Stretches are discarded if they are < min_len.
    """
    stretches = []

    tmp_stretch = []
    stretch_start = 0
    for (idx, val) in enumerate(coverage):
        if val < threshold:
            if not tmp_stretch:
                stretch_start = idx
            tmp_stretch.append(val)
        else:
            if len(tmp_stretch) >= min_len:
                stretches.append((stretch_start, len(tmp_stretch)))
            tmp_stretch = []
            stretch_start = 0
    # If the end of the list if a stretch:
    if len(tmp_stretch) >= min_len:
        stretches.append((stretch_start, len(tmp_stretch)))
    return stretches

def other_frames(frame):
    """ Return the 2 other frames. """
    frame_map = {
        0:[1, 2],
        1:[0, 2],
        2:[0, 1]
    }
    return frame_map[frame]

def sum_keys(dics):
    """ From all dict in 'dics', all the unique keys are
    seeked. Then, the count of all these keys is aggregated
    in one dictionary by a sum.
    ALL KEYS MUST HAVE VALUE OF TYPE <INT>.
    """
    uniq_keys = set([k for d in dics for k in d])
    agg_dict = dict.fromkeys(uniq_keys, 0)
    for dic in dics:
        for key in dic:
            agg_dict[key] += dic[key]
    return agg_dict

def extend_list_keys(to_dic, from_dic, keys):
    """ Extend the 'keys' of 'to_dic' with the
    value found in from_dic.
    """
    for key in keys:
        try:
            to_dic[key].extend(from_dic[key])
        except KeyError as ke:
            print("%s not found when trying to extend." % (ke))

def insert_list_keys(to_dic, from_dic, keys):
    """ Insert at the beginning of the 'keys' lists of 'to_dic' the lists in
    'from_dic' keys.
    """
    for key in keys:
        try:
            to_dic[key][0:0] = from_dic[key]
        except KeyError as ke:
            print("%s not found when trying to insert list." % (ke))

def get_genomic_relative_map(positions):
    """ From a list of genomic positions, generate a dictionary
    with two entries: genomic keys and relative keys to have
    their correspondance.
    """
    gen_rel_map = {"genomic_keys":{}, "relative_keys":{}}

    for rel, gen in enumerate(positions):
        gen_rel_map["genomic_keys"][gen] = rel
        gen_rel_map["relative_keys"][rel] = gen

    return gen_rel_map


def merge_mixt_exons(tuples):
    """ Merge element in the list tuples if they are separated by only 1nt.
    [(10, 20, +), (21, 30, +)] => [(10, 30, +)]
    # Use namedtuples for tupview to make this code better
    """
    for i, first in enumerate(tuples):
        for j, second in enumerate(tuples[i+1:], i+1):
            if abs(first.end - second.start) == 1:
                tuples[i] = MIN_COORDS(first.start,
                                       tuples.pop(j).end,
                                       first.strand)
            else:
                continue
            return merge_mixt_exons(tuples)


def interval_around(from_position, left, right):
    """ Return list with relative and absolute position of the interval
    around from position. """
    positions = {
        "relative":[],
        "absolute":[]
    }
    for i,j in zip(range(left * -1, right + 1),
                   range(from_position-left, from_position+right+1)):
        positions["relative"].append(i)
        positions["absolute"].append(j)
    return positions

def get_value(details_dict, kmer):
    """ Extract the value associated with the kmer in a dictionary. """
    try:
        return details_dict[kmer]
    except KeyError:
        return 0

def get_kmers(details_dict, found_kmers):
    """ Populate the found_kmers list if a new kmer is found as key
    of details_dict.
    """
    for kmer in details_dict:
        if kmer not in found_kmers:
            found_kmers.append(kmer)
        else:
            continue

def add_kmers_columns(pos_cov_df, colname="details"):
    """ Define a column for all kmers detected in the 'colname' column,
    containing dictionary with kmers as key and their counts as value.
    """
    found_kmers = []
    pos_cov_df[colname].apply(lambda x: get_kmers(x, found_kmers))
    for kmer in found_kmers:
        pos_cov_df[str(kmer)] = pos_cov_df[colname].apply(
            lambda x: get_value(x, kmer))
    return found_kmers


def init_poscov():
    """ Return a poscov dictionary. """
    return {"positions":[], "coverages":[], "details":[]}


def poscov_to_dataframe(pos_cov):
    """ Return a dataframe version of pos_cov dictionary.
    For each kmer, a new column is created.
    """
    df = pd.DataFrame(pos_cov)
    found_kmers = add_kmers_columns(df)
    return (df, found_kmers)


def scale_huge_peaks(coverages, scaling_factor="10*mean"):
    """From a list with coverages, the huge pics (> scaling_factor)
    are scaled.
    """
    # Counter of scaled values
    number_of_scaled = 0

    # Remove na to compute median
    all_no_0 = [val for val in coverages if val > 0]

    # max_value = 10 * mean(all_no_0)
    max_value = eval(scaling_factor + "(all_no_0)")

    scaled_values = []
    for original_val in coverages:
        if original_val > max_value:
            scaled_values.append(max_value)
            number_of_scaled += 1
        else:
            scaled_values.append(original_val)

    print(" -> Scaling \033[3;32mENABLED\033[0m \
    (\033[3;33m%s\033[0m pics scaled to \033[3;33m%s\033[0m)" % (
        number_of_scaled, max_value))

    return scaled_values


def rpm_normalisation(coverages, bank_size):
    """ Compute the RPM normalisation"""
    print(" -> RPM \033[3;32mENABLED\033[0m")

    rpm_values = []
    for original_val in coverages:
        rpm_val = (original_val / bank_size) * 1000000
        rpm_values.append(rpm_val)
    return rpm_values


def poscov_frames_signal(frames_signals, ref_frame, pos_cov):
    """ Populate the frames_signals list of coverages for each frame. """
    max_cov = 0
    for ((idx, pos), details, cov) in zip(enumerate(pos_cov["positions"]),
                                          pos_cov["details"],
                                          pos_cov["coverages"]):
        pos_frame = idx % 3
        # Adjust the frame to have frame 0 = ATG frame
        adj_frame = FRAMES_CORREL[ref_frame][pos_frame]

        for frame in frames_signals:
            frames_signals[frame]["positions"].append(pos)
            if frame == adj_frame:
                frames_signals[frame]["coverages"].append(cov)
                frames_signals[frame]["details"].append(details)
                if cov >= max_cov: #Can use Total...so coverages=total..rm cov?
                    max_cov = cov
            else:
                frames_signals[frame]["coverages"].append(0)
                frames_signals[frame]["details"].append({TOTAL_KEYNAME:0})
    return max_cov
