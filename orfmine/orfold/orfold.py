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
from contextlib import ExitStack
from datetime import datetime
from functools import partial
import importlib.util
import importlib.resources
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
from types import ModuleType
from typing import Dict, List, Union
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
        print("Warning:", tango_output)
    else:
        try:
            aa_seq, b_aggregation, h_aggregation = orfold_utils.read_tango_seq(tmp_name + '.txt')
            TANGO_portion = orfold_utils.calculate_proportion_of_seq_aggregable(b_aggregation)
        except:
            # The txt file was never created :(
            TANGO_portion = "NaN"
            print("Error: tango tmp file might not have been created. Process aborted...")
            exit(0)
       
    if os.path.exists(tmp_name + ".txt"):
        if not to_keep:
            os.system("rm " + tmp_name + ".txt")
        else:
            os.system("mv " + tmp_name + ".txt " + seqid + ".txt" )
            os.system("mv " + seqid + ".txt ./TANGO/")
    
    return TANGO_portion


def make_tmp_directories(out_path: Union[str, Path], to_keep: bool=False):
    if to_keep:
        tango_path = Path(out_path) / "TANGO"
        tango_path.mkdir(exist_ok=True, parents=True)
    

def import_optional_tools(options: str, path_tango: str = None, path_iupred: str = None):
    """
    Imports optional tools based on the options provided.

    Args:
        options (str): The options specifying which tools to import ('I' for IUPred, 'T' for Tango).
        path_tango (str): Path to the Tango executable.
        path_iupred (str): Path to the IUPred directory.

    Returns:
        Tuple: The imported IUPred library and Tango executable path.
    """
    iupred2a_lib = None
    tango_executable = None

    if "I" in options:
        # Load IUPred library
        spec = importlib.util.spec_from_file_location("iupred2a_lib", str(Path(path_iupred) / "iupred2a_lib.py"))
        iupred2a_lib = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(iupred2a_lib)

    if "T" in options:
        # Validate Tango path
        if not Path(path_tango).exists():
            raise FileNotFoundError(f"Tango executable '{path_tango}' does not exist.")
        tango_executable = path_tango

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

    if score == "NaN":
        print("Error in computing the score for option {option}. Process aborted...")
        exit(0)

    return round(score, 3), pvalue, mean


def get_scores_from_tabline(line: str, options: str):
    opt_mapping = {opt:i for (i, opt) in enumerate(SUFFIX_MAP, start=1) if opt in options}

    return  {opt:line.split()[idx] for opt, idx in opt_mapping.items()}


def open_gff_file(gff_basename: str, options: str):
    gff_filenames = [f"{gff_basename}{suffix}" for (opt, suffix) in SUFFIX_MAP.items() if opt in options]

    opened_files = {opt: open(_file, 'w') for opt, _file in zip(options, gff_filenames)}

    # Set the signal ensure file will be closed on abrupt program stop (ctrl-c)
    signal.signal(signal.SIGINT, partial(signal_handler, opened_files))

    return opened_files


def annotate_gff(scores: Dict, outpath: Union[str, Path], out_basename: str, gff_template: Union[str, Path]="", options: str="H"):
    # set and open gff filenames to be written
    gff_basename = Path(outpath) / out_basename
    opened_files = open_gff_file(gff_basename=str(gff_basename), options=options)

    # iterate over each line of the gff template and process IDs found in scores
    with open(gff_template, "r") as template_file:
        for line in template_file:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            seqid = line.split('ID=')[-1].split(";")[0]
            if seqid in scores:
                for opt in options:
                    value = float(scores[seqid][opt])
                    maximum = 10 if opt == "H" else 1
                    minimum = -10 if opt == "H" else 0

                    new_gffline = orfold_utils.set_gffline_color(line_template=line, value=value, minimum=minimum, maximum=maximum)
                    opened_files[opt].write(new_gffline)
    
    # close opened files
    for _file in opened_files.values():
        _file.close()


def process_orf(orf: str, seq: str, options: str, iupred2a_lib: ModuleType=None, tango_path: Union[str, Path]="", to_keep: bool=False):
    scores = {}
    for option in ["H", "I", "T"]:
        scores[option] = "NaN"  # initialize scores to NaN for each property

        # skip score computation if option not asked
        if option not in options:
            continue

        # compute asked property score
        seqid = orf if option == "T" else None
        scores[option], _, _ = calculate_score(option=option, seq=seq, iupred2a_lib=iupred2a_lib, tango_path=tango_path, to_keep=to_keep, seqid=seqid)

    return scores


def process_fasta_file(fasta_file: Union[str, Path], out_path: Union[str, Path], options: str="H", sample_size: int=-1, to_keep: bool=False, path_tango: str = None, path_iupred: str = None):
    # import external optional softwares
    iupred2a_lib, tango_path = import_optional_tools(options=options, path_tango=path_tango, path_iupred=path_iupred)

    # get fasta file base name
    fasta_basename = Path(fasta_file).stem

    # set outfile
    scores_tab_file = out_path / f"{fasta_basename}.tab"

    # get table output format of orfold
    out_format = orfold_utils.get_orfold_out_format(max_len_head=60)

    # get sequences in the fasta file
    fasta_sequences = orfold_utils.fasta_generator(filename=fasta_file) if sample_size == -1 else orfold_utils.reservoir_sampling_fasta(filename=fasta_file, k=sample_size)

    all_scores = {}
    with open(scores_tab_file, "w") as fw_output:
        fw_output.write(out_format.format("Seq_ID", "HCA", "Disord", "Aggreg"))

        # iterate over each sequence
        for header, sequence in fasta_sequences:                
            scores = process_orf(orf=header, seq=sequence, options=options, iupred2a_lib=iupred2a_lib, tango_path=tango_path, to_keep=to_keep)

            # write scores in the table output
            fmt_scores = [orfold_utils.format_with_n_decimals(scores[opt]) if scores[opt] != "NaN" else "NaN" for opt in "HIT"]
            fw_output.write(out_format.format(header, *fmt_scores))
            all_scores[header] = scores

    return all_scores


def run_orfold(fasta_file: Union[str, Path], out_path: Union[str, Path], options: str = "H", sample_size: Union[int, str] = None, gff_template: Union[str, Path] = "",to_keep: bool = False,path_tango: str = None, path_iupred: str = None  ):
    out_path = Path(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    # create optional directories for tango if necessary
    make_tmp_directories(out_path=out_path, to_keep=to_keep)

    # process fasta file
    all_scores = process_fasta_file(fasta_file=fasta_file, out_path=out_path, options=options, sample_size=sample_size, to_keep=to_keep)

    # annotate gff if required
    if gff_template:
        annotate_gff(scores=all_scores, outpath=out_path, out_basename=Path(fasta_file).stem, options=options, gff_template=gff_template)


def run_orfold_containerized(parameters: arguments.argparse.Namespace):
    # Vérifier que les chemins sont fournis
    extra_bindings = {}
    if parameters.path_tango:
        extra_bindings[parameters.path_tango] = '/opt/tango'
    if parameters.path_iupred:
        extra_bindings[parameters.path_iupred] = '/opt/iupred'

    # Instantiate ContainerCLI
    cli = ContainerCLI(
        input_args=["--faa"],
        output_arg="--out",
        args=parameters,
        image_base=DOCKER_IMAGE,
        prog="orfold",
        container_type="docker" if parameters.docker else "singularity",
        extra_bindings=extra_bindings,
        dev_mode=parameters.dev,
        package_binding={"orfmine": "/home/orfuser/orfmine/orfmine"}
    )

    cli.show()
    if not parameters.dry_run:
        cli.run()

def main():
    parameters = arguments.get_args()

    if parameters.docker or parameters.singularity:
        run_orfold_containerized(parameters=parameters)
    else:
        start_time = datetime.now()

        run_orfold(
            fasta_file=parameters.faa,
            out_path=parameters.out,
            options=parameters.options,
            gff_template=parameters.gff,
            sample_size=parameters.sample,
            to_keep=parameters.keep,
        )

        end_time = datetime.now() 
        print('\nDuration: {}'.format(end_time - start_time))


if __name__ == "__main__":
    main()





