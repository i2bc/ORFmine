#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 01:21:05 2020

@author: christospapadopoulos


#@TODO To be integrated later!!!
if "H" in parameters.options and parameters.barcodes:
    # First we calculate ALL the Barcodes at ones:
    barcodes = calculate_HCA_barcodes(sequences=sequences)

    with open(out_path / str(fasta_basename + ".barcodes"), "w") as barcw:
        for i in barcodes:
            barcw.write(">{}\n{}\n".format(i,barcodes[i]))

"""
from datetime import datetime
from functools import partial
import importlib.util
import importlib.resources
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import ModuleType
from typing import Dict, Union
import signal
import warnings


# Temporarily redirect stderr | this is a hack done to remove the error print message from hca
from io import StringIO
original_stderr = sys.stderr
sys.stderr = StringIO()

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from pyHCA.core.annotateHCA import _annotation_aminoacids
    from pyHCA.core.classHCA import compute_disstat

# Restore stderr
sys.stderr = original_stderr

from orfmine import DOCKER_IMAGE
from orfmine.orfold.lib import utils as orfold_utils
from orfmine.orfold.lib import arguments
from orfmine.utilities.container import ContainerCLI


TANGO_EXEC = {
    "darwin": "tango2_3_1",
    "win32": "Tango.exe",
    "linux": "tango_x86_64_release",
}


SUFFIX_MAP = {
    "H": "_HCA.gff",
    "I": "_IUPRED.gff",
    "T": "_TANGO.gff"
}

def signal_handler(opened_files, sig, frame):
    print('You pressed Ctrl+C! Closing files...')
    for _, file_dict in opened_files.items():
        file_dict["file"].close()
    exit(1)




def calculate_HCA_barcodes(sequences):
    """
    We get the Barcode sequences of HCA for all the sequences
    """
    print("Calculating the HCA barcodes ... ")

    from pyHCA import HCA
    hca = HCA(seq=list(sequences.values()),querynames=list(sequences.keys()))
    # You can get the barcode sequence of ONE sequnce using the module get_hca_barcode
    # If you want to take the barcode of all the sequnces do a loop like:
    barcodes = {}
    for orf in list(sequences.keys()):
        barcodes[orf] = orfold_utils.get_hca_barcode(hca = hca, orf = orf)

    print("Barcodes calculation : DONE")
    print("\n")

    return barcodes

    
def calculate_tango_one_sequence(tango_path, seq, seqid, to_keep):
    tf = tempfile.NamedTemporaryFile(prefix="tango")
    tmp_name = tf.name.split("/")[-1]

    tango_command = tango_path + " " + tmp_name + " ct=\"N\" nt=\"N\" ph=\"7.4\" te=\"298\" io=\"0.1\" seq=\"" + seq + "\""
    
    process = subprocess.Popen(tango_command, stdout=subprocess.PIPE, stderr=None, shell=True)
    tango_output = process.communicate()    
       
    if str(tango_output[0]).replace("\\n", "").replace("\'", "").strip() == "b88, File not properly written, try writing it up again,": 
        b_aggregation = 'None'
        TANGO_portion = "NaN"
    else:
        try:
            aa_seq, b_aggregation, h_aggregation = orfold_utils.read_tango_seq(tmp_name + '.txt')
            TANGO_portion = orfold_utils.calculate_proportion_of_seq_aggregable(b_aggregation)
        except:
            # The txt file was never created :(
            TANGO_portion = "NaN"
       
    if os.path.exists(tmp_name + ".txt"):
        if not to_keep:
            os.system("rm " + tmp_name + ".txt")
        else:
            os.system("mv " + tmp_name + ".txt " + seqid + ".txt" )
            os.system("mv " + seqid + ".txt ./TANGO/")
    else:
        pass
    
    return TANGO_portion


def make_tmp_directories(out_path: Union[str, Path], to_keep: bool=False):
    if to_keep:
        tango_path = Path(out_path) / "TANGO"
        tango_path.mkdir(exist_ok=True, parents=True)
    

def import_optional_tools(options: str):
    iupred2a_lib = None
    tango_executable = None
    
    # get path to external softwares given in config.ini
    external_softwares = orfold_utils.read_config_file()
    
    # Check if the tools asked for the analysis are well Installed
    if "I" in options:
        is_valid, error_message = orfold_utils.check_path(path=external_softwares["iupred"], software="iupred")
        if not is_valid:
            print(orfold_utils.error_missing_softwares(software="iupred").format(error_message))
            exit()

        # import iupred2a_lib from its location
        spec = importlib.util.spec_from_file_location("iupred2a_lib", str(Path(external_softwares["iupred"]) / "iupred2a_lib.py"))
        iupred2a_lib = importlib.util.module_from_spec(spec)
        sys.modules["iupred2a_lib"] = iupred2a_lib
        spec.loader.exec_module(iupred2a_lib)

    if "T" in options:
        is_valid, error_message = orfold_utils.check_path(path=external_softwares["tango"], software="tango")
        if not is_valid:
            print(orfold_utils.error_missing_softwares(software="tango").format(error_message))
            exit()

        tango_path = Path(external_softwares["tango"])
        if tango_path.is_dir():
            tango_executable = str(tango_path / TANGO_EXEC[sys.platform])
        elif tango_path.is_file() and tango_path.name in TANGO_EXEC.values():
            tango_executable = str(tango_path)            

    return iupred2a_lib, tango_executable


def calculate_score(option, seq, iupred2a_lib=None, tango_path=None, to_keep: bool=False, seqid=None,):
    score, pvalue, mean = ["NaN"]*3

    if option == "H":
        score, pvalue = compute_disstat(0, len(seq), _annotation_aminoacids(seq=seq, method="domain", verbose=False)["cluster"])

    elif option == "I":
        pre_score  = iupred2a_lib.iupred(seq=seq, mode="short")[0]  # score
        score = orfold_utils.calculate_proportion_of_seq_disordered(pre_score)  # iupred_portion
        mean = round(sum(pre_score) / len(pre_score), 2)

    elif option == "T":
        score = calculate_tango_one_sequence(tango_path=tango_path, seq=seq, seqid=seqid, to_keep=to_keep)

    return round(score, 3), pvalue, mean


def write_gff_line(outfile, gff_dico, orf, value, gff_filename, nb_cols=20, minimum=0, maximum=1):
    try:
        new_gff_line = orfold_utils.change_color_in_gff_line(gff_dico=gff_dico, orf=orf, value=value, nb_cols=nb_cols, minimum=minimum, maximum=maximum)
        outfile.write(new_gff_line)
    except:
        print('An error occured at the writing of the {} file for the orf: {}'.format(gff_filename, orf))
        pass

def open_gff_file(gff_template: Union[str, Path], options: str, basename: str):
    opened_files = {}
    if gff_template:
        for option, suffix in SUFFIX_MAP.items():
            if option in options:
                opened_files[option] = {
                    "file": open(basename + suffix, 'w'),
                    "gff_filename": basename + suffix,
                    "gff_dict": orfold_utils.read_gff_file(gff_file=gff_template)
                }

    # Set the signal ensure file will be closed on abrupt preogram stop (ctrl-c)
    signal.signal(signal.SIGINT, partial(signal_handler, opened_files))

    return opened_files

def process_orf(orf: str, seq: str, options: str, iupred2a_lib: ModuleType=None, tango_path: Union[str, Path]="", opened_files: dict={}, to_keep: bool=False):
    scores = {}
    for option in ["H", "I", "T"]:
        scores[option] = "NaN"  # initialize scores to NaN for each property

        # skip score computation if option not asked
        if option not in options:
            continue

        # compute asked property score
        seqid = orf if option == "T" else None
        scores[option], _, _ = calculate_score(option=option, seq=seq, iupred2a_lib=iupred2a_lib, tango_path=tango_path, to_keep=to_keep, seqid=seqid)

        # write new gff with scores information
        if option in opened_files:
            write_gff_line(
                outfile=opened_files[option]["file"],
                gff_dico=opened_files[option]["gff_dict"],
                orf=orf,
                value=scores[option],
                gff_filename=opened_files[option]["gff_filename"],
                maximum=10 if option == "H" else 1,
                minimum=-10 if option == "H" else 0
            )

    return scores


def process_fasta_file(fasta_file: Union[str, Path], out_path: Union[str, Path], options: str="H", sample_size: Union[int, str]=None, gff_template: Union[str, Path]="", to_keep: bool=False):
    # import external optional softwares
    iupred2a_lib, tango_path = import_optional_tools(options=options)

    # get fasta file base name
    fasta_basename = Path(fasta_file).stem

    # get sequences in the fasta file
    sequences = orfold_utils.get_sequences(fasta_file=fasta_file, sample_size=sample_size)

    # get table output format of orfold
    out_format = orfold_utils.get_orfold_out_format(max_len_head=len(max(sequences.keys(), key=len)))

    with open(out_path / (fasta_basename + ".tab"), "w") as fw_output:
        fw_output.write(out_format.format("Seq_ID", "HCA", "Disord", "Aggreg"))

        opened_files = {}
        if gff_template:
            opened_files = open_gff_file(gff_template=gff_template, options=options, basename=fasta_basename)

        for orf, seq in sequences.items():
            scores = process_orf(orf=orf, seq=seq, options=options, iupred2a_lib=iupred2a_lib, tango_path=tango_path, opened_files=opened_files, to_keep=to_keep)

            # write scores in the table output
            scores = [orfold_utils.format_with_n_decimals(scores[opt]) if scores[opt] != "NaN" else "NaN" for opt in SUFFIX_MAP]
            fw_output.write(out_format.format(orf, *scores))

        # close opened files
        for _f in opened_files.values():
            _f["file"].close()


def run_orfold(fasta_file: Union[str, Path], out_path: Union[str, Path], options: str="H", sample_size: Union[int,str]=None, gff_template: Union[str, Path]="", to_keep: bool=False):
    out_path = Path(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    # create optional directories for tango if necessary
    make_tmp_directories(out_path=out_path, to_keep=to_keep)

    # process fasta file
    process_fasta_file(fasta_file=fasta_file, out_path=out_path, options=options, gff_template=gff_template, sample_size=sample_size, to_keep=to_keep)


def run_orfold_containerized(parameters: arguments.argparse.Namespace):

    # set list of flags related to input files
    input_args = ["--faa"]
    if parameters.gff:
        input_args += ["--gff"]

    # set external softwares paths
    software_bindings = orfold_utils.read_config_file()

    # set flag related to output path/file
    output_arg = "--out"

    # instantiate containerCLI handler
    cli = ContainerCLI(
            input_args=input_args,
            output_arg=output_arg,
            args=parameters,
            image_base=DOCKER_IMAGE,
            prog="orfold",
            container_type="docker" if parameters.docker else "singularity",
            software_bindings=software_bindings
        )

    cli.show()
    if not parameters.dry_run:
        cli.run()


def main():
    start_time = datetime.now()

    parameters = arguments.get_args()

    if parameters.docker or parameters.singularity:
        run_orfold_containerized(parameters=parameters)
    else:
        run_orfold(fasta_file=parameters.faa, out_path=parameters.out, options=parameters.options, gff_template=parameters.gff, sample_size=parameters.sample, to_keep=parameters.keep)

    end_time = datetime.now()
    print('\nDuration: {}'.format(end_time - start_time))


if __name__ == "__main__":
    main()





