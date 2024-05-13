#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:04:13 2019

@author: c.papadopoulos
"""
# todo add headers
import argparse
from datetime import datetime
import re

from orfmine.orfribo.lib.loaders import Gff
#import orfmine.orfribo.lib.biobjects 

#from loaders import Gff
#import biobjects

def get_args():
    """
    :return: parameters
    """
    parser = argparse.ArgumentParser(description='Bam2Read')
    parser.add_argument("-gff",
                        type=str,
                        required=True,
                        nargs="?",
                        help="GFF annotation file of your regions of interest")
    parser.add_argument("-bam",
                        type=str,
                        required=True,
                        nargs="?",
                        help="FASTA file containing the genome sequence")
    parser.add_argument("-features_include",
                        type=str,
                        required=False,
                        nargs="*",
                        default=['all'],
                        help="Annotation features to be considered (By definition is all)")
    parser.add_argument("-features_exclude",
                        type=str,
                        required=False,
                        nargs="?",
                        default=["None"],
                        help="Annotation features not to be considered (By definition is None)")
    parser.add_argument("-outpath",
                        type=str,
                        # action='store',
                        required=True,
                        nargs="?",
                        help="The path to save the dictionary output")
    parser.add_argument("-outname",
                        type=str,
                        # action='store',
                        required=True,
                        nargs="?",
                        help="The name you want the output to have")
    parser.add_argument("-shift",
                        type=int,
                        # action='store',
                        required=False,
                        default=12,
                        nargs="?",
                        help="The shift of the P-site calculation (default = 12)")
    parser.add_argument("-kmer",
                        type=int,
                        action='store',
                        default=['Total'],
                        required=False,
                        nargs="+",
                        help="The Kmer size to consider for counting")
    args = parser.parse_args()
    return args


def coverage_rate(coverage):
    """ Return the rate of position with at least a read in coverage (list).
    """
    return round(sum((1 for x in coverage if x > 0)) / len(coverage), 5)


def correct_positions_arrays(reads_p0,reads_p1,reads_p2):
    arrays = [reads_p0,reads_p1,reads_p2]
    arrays_lengths = [len(i) for i in arrays]
    for i,arr in enumerate(arrays):
        if arrays_lengths[i] < max(arrays_lengths):
            for j in range(max(arrays_lengths) - arrays_lengths[i]):
                arrays[i].append(0)
    return(arrays[0],arrays[1],arrays[2])


def count_percentage_reads_to_file(file_output, elements_in, elements_out, gff_iterator, bam, shift,kmers):
    print("BAM : %s" % bam)
    with open(file_output + "_reads.tab", "w") as wtab, \
            open(file_output + "_periodicity_start.tab", "w") as wpstart, \
            open(file_output + "_periodicity_stop.tab", "w") as wpstop, \
            open(file_output+ "_periodicity_all.tab","w") as wper:
        start = datetime.now()
        # Initialize counters
        features_found = 0
        features_count = 0
        # Write the title of the file
        wtab.write("{:<50s}\t{:>15s}\t{:>15s}\t{:>15s}\t{:>15s}\t{:>10s}\t{:>10s}\t{:>10s}\n".format("Seq_ID","Num_reads","Num_p0","Num_p1","Num_p2","Perc_p0","Perc_p1","Perc_p2"))
        for x, feature in enumerate(gff_iterator):
            # Filtering of the features
            if re.search(elements_out, feature.ftype) and not re.search(elements_in, feature.ftype):
                continue

            if re.search(elements_out, feature.ftype) and re.search(elements_in, feature.ftype):
                print("{} include and exclude at the same time!".format(feature.ftype))
                exit()

            if not re.search(elements_in, feature.ftype) and elements_in != "(all)":
                continue
            features_found += 1
            # Check if there is a problem with the P-site reduction
            try:
                # Coverage of the feature in the 3 frames
                coverage_by_frame = feature.frames_coverage(bam,shift=shift,kmers=kmers)
                reads_p0 = coverage_by_frame[0]
                reads_p1 = coverage_by_frame[1]
                reads_p2 = coverage_by_frame[2]

                # Number of reads by frame
                nb_reads_p0 = sum(coverage_by_frame[0])
                nb_reads_p1 = sum(coverage_by_frame[1])
                nb_reads_p2 = sum(coverage_by_frame[2])

                nb_reads_gene = nb_reads_p0 + nb_reads_p1 + nb_reads_p2
            except:
                # If there is a problem with the P-site reduction we do not consider this ORF
                continue

            # Percentage of reads by frame
            try:
                perc_reads_p0 = round(nb_reads_p0 / nb_reads_gene * 100, 2)
                perc_reads_p1 = round(nb_reads_p1 / nb_reads_gene * 100, 2)
                perc_reads_p2 = round(nb_reads_p2 / nb_reads_gene * 100, 2)
            except ZeroDivisionError:
                perc_reads_p0 = 0.0
                perc_reads_p1 = 0.0
                perc_reads_p2 = 0.0

            # Write on the reads count tab the total number of reads, the number and the percentage per phase for the feature
            wtab.write("{:<50s}\t{:>15d}\t{:>15d}\t{:>15d}\t{:>15d}\t{:>10.2f}\t{:>10.2f}\t{:>10.2f}\n".format(feature.ID, nb_reads_gene, nb_reads_p0,
                                                                           nb_reads_p1, nb_reads_p2, perc_reads_p0,
                                                                           perc_reads_p1, perc_reads_p2))
            features_count += 1


            # If the positions of the frame are not equals, we correct them
            if not len(reads_p0) == len(reads_p1) == len(reads_p2):
                reads_p0,reads_p1,reads_p2 = correct_positions_arrays(reads_p0,reads_p1,reads_p2)

            #Periodicity of the feature of all length
            for aa in range(len(reads_p0)):
                wper.write("{}\t{}\t{}\t{}\t{}\n".format(feature.ID, aa+1, reads_p0[aa], reads_p1[aa], reads_p2[aa]))

            # Periodicity of the feature of length superior to 50
            if len(reads_p0) > 50:
                # We write the periodicity of the first 50 AA positions
                wpstart.write('{}\tp0\t'.format(feature.ID))
                for i in range(0, 51):
                    wpstart.write('{}\t'.format(reads_p0[i]))
                wpstart.write('\n')
                wpstart.write('{}\tp1\t'.format(feature.ID))
                for i in range(0, 51):
                    wpstart.write('{}\t'.format(reads_p1[i]))
                wpstart.write('\n')
                wpstart.write('{}\tp2\t'.format(feature.ID))
                for i in range(0, 51):
                    wpstart.write('{}\t'.format(reads_p2[i]))
                wpstart.write('\n')
                # We write the periodicity of the last 50 AA positions
                wpstop.write('{}\tp0\t'.format(feature.ID))
                for i in range(len(reads_p0) - 50, len(reads_p0)):
                    wpstop.write('{}\t'.format(reads_p0[i]))
                wpstop.write('\n')
                wpstop.write('{}\tp1\t'.format(feature.ID))
                for i in range(len(reads_p0) - 50, len(reads_p0)):
                    wpstop.write('{}\t'.format(reads_p1[i]))
                wpstop.write('\n')
                wpstop.write('{}\tp2\t'.format(feature.ID))
                for i in range(len(reads_p0) - 50, len(reads_p0)):
                    wpstop.write('{}\t'.format(reads_p2[i]))
                wpstop.write('\n')
        end = datetime.now()

    print("Total ORFs\t:\t",features_found,"\t(100.00 %)")
    print("ORFs kept \t:\t",features_count,"\t(",str(round((features_count/features_found*100),2)),"%)")

    print("Duration : {}\n {} features corresponding the selection.\nThe mean time per selected feature is : {}\n".format(end-start, features_found,(end-start)/features_found))
 


def BAM2Reads(rname, gff_file, kmer, outname, elements_in=None, elements_out=None):
    if elements_in is None:
        elements_in = ["all"]
    if elements_out is None:
        elements_out = ["None"]

    elements_in = "(" + ")|(".join(elements_in) + ")"
    elements_out = "(" + ")|(".join(elements_out) + ")"
    print('Read the GFF file')
    start_time_gff = datetime.now()
    gff = Gff(gff_file, all_as_high=True)
    gff_iterator = sorted(gff.all_features(cast_into="Igorf"))
    end_time_gff = datetime.now()
    output_log = ""
    output_log = output_log+ 'Duration gff: {}'.format(end_time_gff - start_time_gff)
    print('GFF file read \t DONE')
    for size in kmer:
        if outname == "phasing":  # ORFphase
            file_input = "./kmer_{}/{}_kmer_{}_phasing".format(str(size), rname, str(size))
        else:  # ORFribomap
            file_input = "./kmer_{}/{}_kmer_{}_all".format(str(size), rname, str(size))
        bam_file = file_input + "_sorted_mapped.bam"
        file_output = "./kmer_{}/{}_kmer_{}_{}".format(str(size), rname, str(size), outname)
        output_log = output_log + "Kmer : {}\n".format(size)+count_percentage_reads_to_file(file_output, elements_in, elements_out, gff_iterator, bam_file)
    end_time_all = datetime.now()
    
    return output_log+ '\n\nDuration b2r: {}'.format(end_time_all - start_time_gff)


def main():
    parameters = get_args()
    print('Read the GFF file')
    GFF = Gff(parameters.gff, all_as_high=True)
    gff_iterator = sorted(GFF.all_features(cast_into="Igorf"))
    print('GFF file read \t DONE')
    elements_in = "(" + ")|(".join(parameters.features_include) + ")"
    elements_out = "(" + ")|(".join(parameters.features_exclude) + ")"
    file_output = parameters.outpath + "/" + parameters.outname

    count_percentage_reads_to_file(file_output, elements_in, elements_out, gff_iterator, bam = parameters.bam, shift = parameters.shift,kmers=parameters.kmer)


if __name__ == "__main__":
    main()
