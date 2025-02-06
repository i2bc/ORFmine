#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 15:38:13 2022

@author: christospapadopoulos
"""

import argparse
from pathlib import Path


def get_args():
    """
    :return: parameters
    """
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument("-tables",
                        type=str,
                        required=True,
                        nargs="*",
                        help="Tables of reads to be concatenated")
    parser.add_argument("-outpath",
                        type=str,
                        required=True,
                        nargs="?",
                        help="The path to save the dictionary output")
    parser.add_argument("-outname",
                        type=str,
                        required=True,
                        nargs="?",
                        help="The name you want the output to have")

    args = parser.parse_args()
    return args


def concatenate(tables):
    dico = {}
    for tab in tables:
        with open(tab,'r') as file:
            for line in file:
                ID = line.split()[0].strip()
                if ID == "Seq_ID":
                    continue
                if ID not in dico:
                    dico[ID] = {
                        "Num_reads": 0,
                        "Num_p0": 0,
                        "Num_p1": 0,
                        "Num_p2" :0,
                        "Perc_p0": 0.0,
                        "Perc_p1": 0.0,
                        "Perc_p2":0.0
                    }

                dico[ID]["Num_reads"] = dico[ID]["Num_reads"] + int(line.split()[1])
                dico[ID]["Num_p0"] = dico[ID]["Num_p0"] + int(line.split()[2])
                dico[ID]["Num_p1"] = dico[ID]["Num_p1"] + int(line.split()[3])
                dico[ID]["Num_p2"] = dico[ID]["Num_p2"]+ int(line.split()[4])

                try:
                    dico[ID]["Perc_p0"]  = float(dico[ID]["Num_p0"] / dico[ID]["Num_reads"]*100)
                except:
                    pass
                try:
                    dico[ID]["Perc_p1"]  = float(dico[ID]["Num_p1"] / dico[ID]["Num_reads"]*100)
                except:
                    pass
                try:
                    dico[ID]["Perc_p2"]  = float(dico[ID]["Num_p2"] / dico[ID]["Num_reads"]*100)
                except:
                    pass

    return dico
    
def write_output(concatenated, out_filename):
    with open(out_filename, "w") as fw:
        fw.write("{:<50s}\t{:>15s}\t{:>15s}\t{:>15s}\t{:>15s}\t{:>10s}\t{:>10s}\t{:>10s}\n".format("Seq_ID","Num_reads","Num_p0","Num_p1","Num_p2","Perc_p0","Perc_p1","Perc_p2"))
        for ID in concatenated:
            fw.write(
                "{:<50s}\t{:>15d}\t{:>15d}\t{:>15d}\t{:>15d}\t{:>10.2f}\t{:>10.2f}\t{:>10.2f}\n".format(
                    ID,concatenated[ID]["Num_reads"],
                    concatenated[ID]["Num_p0"],concatenated[ID]["Num_p1"],concatenated[ID]["Num_p2"],
                    concatenated[ID]["Perc_p0"],concatenated[ID]["Perc_p1"],concatenated[ID]["Perc_p2"]))
        


def main():
    parameters = get_args()

    Path(parameters.outpath).mkdir(parents=True, exist_ok=True)
    out_filename = Path(parameters.outpath) / (str(Path(parameters.outname)) + "_reads_concatenated.tab")

    concatenated = concatenate(tables = parameters.tables)
    write_output(concatenated, out_filename=out_filename)
    

if __name__ == "__main__":
    main()
