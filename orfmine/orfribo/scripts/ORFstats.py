#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 16:53:28 2022

@author: christospapadopoulos
"""

import argparse
import statistics
from collections import Counter
from pathlib import Path



def get_args():
    """
    :return: parameters
    """
    parser = argparse.ArgumentParser(description='Bam2Read')
    parser.add_argument("-tab",
                        type=str,
                        required=True,
                        nargs="?",
                        help="Table of reads as generated from BAM2Reads")
    parser.add_argument("-N",
                        type=int,
                        required=False,
                        default=10,
                        nargs="?",
                        help="The minimum number of reads to calculte the phase"),
    parser.add_argument("-out",
                        type=str,
                        required=False,
                        default="./",
                        nargs="?",
                        help="Path to output directory (default='./')")

    args = parser.parse_args()
    return args



def my_mode(sample):
    c = Counter(sample)
    return [k for k, v in c.items() if v == c.most_common(1)[0][1]]


def main():
    parameters = get_args()
    
    output  = Path(parameters.out) / Path(parameters.tab).stem
    
    dico = {"P0": [], "P1": [], "P2": []}

    with open(parameters.tab, 'r') as f:
        next(f)
        for line in f:
            if int(line.split()[1]) >= parameters.N:
                dico["P0"].append(float(line.split()[5]))
                dico["P1"].append(float(line.split()[6]))
                dico["P2"].append(float(line.split()[7]))
            
            
    #  M E A N :
    try:
        mean_P0 = round(statistics.mean(dico["P0"]),2)
    except:
        mean_P0 = 0.0
    try:
        mean_P1 = round(statistics.mean(dico["P1"]),2)
    except:
        mean_P1 = 0.0
    try:
        mean_P2 = round(statistics.mean(dico["P2"]),2)
    except:
        mean_P2 = 0.0
    
    #  M E D I A N :
    try:
        median_P0 = round(statistics.median(dico["P0"]),2)
    except:
        median_P0 = 0.0
    try:
        median_P1 = round(statistics.median(dico["P1"]),2)
    except:
        median_P1 = 0.0
    try:
        median_P2 = round(statistics.median(dico["P2"]),2)
    except:
        median_P2 = 0.0
        
    #  M O D E :
    try:
        mode_P0 = round(list(my_mode(dico["P0"]))[0],2)
    except:
        mode_P0 = 0.0
    try:
        mode_P1 = round(list(my_mode(dico["P1"]))[0],2)
    except:
        mode_P1 = 0.0
    try:
        mode_P2 = round(list(my_mode(dico["P2"]))[0],2)
    except:
        mode_P2 = 0.0
    
    
    with open(str(output) + ".stats", "w") as fw:
        fw.write("{:10s}\t{:>10s}\t{:>10s}\t{:>10s}\n".format("Metric", "P0", "P1", "P2"))
        fw.write("{:10s}\t{:10.2f}\t{:10.2f}\t{:10.2f}\n".format("Mean", mean_P0, mean_P1, mean_P2))
        fw.write("{:10s}\t{:10.2f}\t{:10.2f}\t{:10.2f}\n".format("Median", median_P0, median_P1, median_P2))
        fw.write("{:10s}\t{:10.2f}\t{:10.2f}\t{:10.2f}\n".format("Mode", mode_P0, mode_P1, mode_P2))



            
if __name__ == "__main__":
    main()        
        
        
        
