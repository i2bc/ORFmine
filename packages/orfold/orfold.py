#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 01:21:05 2020

@author: christospapadopoulos
"""

import os,sys,re,random
from datetime import datetime
import subprocess
import argparse
from pathlib import Path
import importlib.util
import tempfile
import sys
import warnings

import seaborn as sns

from packages.orfold.lib import utils as orfold_utils
from matplotlib.colors import to_hex


def get_args():
    """
    Returns:
        Parameters
    """

    parser = argparse.ArgumentParser(description='ORF Foldability Calculation')
    parser.add_argument("-faa",
                        type=str,
                        action='store',
                        required=True, 
                        nargs="*",
                        help="FASTA file containing the amino acid sequences to treat")
    
    parser.add_argument("-gff", 
                        required=False, 
                        type=str,
                        action='store',
                        nargs="*",
                        default=[],
                        help="GFF annotation file")
    
    parser.add_argument("-options",
                        type=list,
                        #action='store',
                        required=True, 
                        nargs="?",
                        default=["H"],
                        help=
                        '''Which properties are to be calculated. 
                             H for HCA (Default)
                             I for IUPred
                             T for Tango''')
    
    parser.add_argument("-plot", 
                        required=False, 
                        type=str,
                        action='store',
                        nargs="*",
                        default=False,
                        help=argparse.SUPPRESS)
    
    parser.add_argument("-barcodes", 
                        required=False, 
                        type=str,
                        #action='store',
                        nargs="?",
                        default=False,
                        help=argparse.SUPPRESS)
    
    parser.add_argument("-keep", 
                        required=False, 
                        type=list,
                        #action='store',
                        nargs="?",
                        default=[],
                        help="Option for keeping the Tango output files")
    
    parser.add_argument("-N", 
                        required=False, 
                        type=str,
                        action='store',
                        nargs="*",
                        default=["all"],
                        help="Size of sample(s) per FASTA file")
    
    args = parser.parse_args()
    return args


def get_root_name_of_files_list(files_list):
    """
    Removes the extentions and path of the files
    """
    files = {}
    for i in files_list:
        files[i.split("/")[-1].split(".")[0]] = i
    return files


def calculate_HCA_barcodes(sequences):
    """
    We get the Barcode sequences of HCA for all the sequences
    """
    from pyHCA import HCA
    hca = HCA(seq=list(sequences.values()),querynames=list(sequences.keys()))
    # You can get the barcode sequence of ONE sequnce using the module get_hca_barcode
    # If you want to take the barcode of all the sequnces do a loop like:
    barcodes = {}
    for orf in list(sequences.keys()):
        barcodes[orf] = orfold_utils.get_hca_barcode(hca = hca, orf = orf)  
    return barcodes

    
def calculate_tango_one_sequence(tango_path, name, to_keep):
    tf = tempfile.NamedTemporaryFile(prefix="tango")
    tmp_name = tf.name.split("/")[-1]

    tango_command = tango_path + " " + tmp_name + " ct=\"N\" nt=\"N\" ph=\"7.4\" te=\"298\" io=\"0.1\" seq=\"" + sequences[name] + "\""
    
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
        if 'T' not in to_keep:
            os.system("rm " + tmp_name + ".txt")
        else:
            os.system("mv " + tmp_name + ".txt " + name + ".txt" )
            os.system("mv " + name + ".txt ./TANGO/")
    else:
        pass
    
    return TANGO_portion


def make_tmp_directories(parameters):
    if "T" in parameters.keep:
        try:
            os.system("mkdir TANGO")
        except:
            pass
    

def make_files_associations(parameters):

    fastas = get_root_name_of_files_list(files_list=parameters.faa)
    # Check if gff files where given:
    if parameters.gff:
        gffs = get_root_name_of_files_list(files_list=parameters.gff)
    else:
        gffs = parameters.gff

    samples = parameters.N
    files_associations = {}

    # First I associate the number of sequences samples
    # It must be given in order!!! OBLIGATORY
    files_sampling = {}
    for n,name in enumerate(fastas):
        try:
            files_sampling[fastas[name]] = samples[n]
        except:
            files_sampling[fastas[name]] = "all"
    # ---------------------------------------- DONE
    result =  all(elem in fastas for elem in gffs)
    if result:
        # All the gff files found an associated fasta file
        print(
             '''
            These are the files associations I can make:''')
        print(
        ''' 
            {:^30s}\t{:^30s}\t{:^30s}
            {:^30s}\t{:^30s}\t{:^30s}'''.format("FASTA" , "GFF" , "Nb sequences","-----","---","------------")
                )
        for n,name in enumerate(fastas):
            try:
                print(
        ''' 
            {:^30s}\t{:^30s}\t{:^30s}'''.format(fastas[name].split("/")[-1],gffs[name].split("/")[-1],files_sampling[fastas[name]]))
                files_associations[fastas[name]] = gffs[name]
            except:
                print(
        ''' 
            {:^30s}\t{:^30s}\t{:^30s}'''.format(fastas[name].split("/")[-1],"",files_sampling[fastas[name]]))
                files_associations[fastas[name]] = ''

    else:
        # There is at least 1 GFF which could not be associated with FASTA
        print('''
             Oups! You provided GFF file(s) which has no correspodance to FASTA''')

        # BUT if we provide the same number of FASTA and GFFs then we can
        # associate them based on the order they passed in the terminal
        if len(fastas) == len(gffs):
            print(
              '''
             BUT i found {} FASTAs and {} GFFs 
             so I will associate them based in the order they are:'''.format(len(fastas),len(gffs)))
            print(
        ''' 
            {:^30s}\t{:^30s}\t{:^30s}
            {:^30s}\t{:^30s}\t{:^30s}'''.format("FASTA" , "GFF" , "Nb sequences","-----","---","------------")
                )
            for n,name in enumerate(fastas):
                print(
        ''' 
            {:^30s}\t{:^30s}\t{:^30s}''' \
              .format(fastas[list(fastas.keys())[n]].split("/")[-1],
                      gffs[list(gffs.keys())[n]].split("/")[-1],
                      files_sampling[fastas[name]]))
                files_associations[fastas[list(fastas.keys())[n]]] = gffs[list(gffs.keys())[n]]

        elif len(gffs) == 1:
            print(
              '''
             BUT i found {} FASTAs and {} unique GFF  
             so I will associate this GFF with ALL the FASTA:'''.format(len(fastas),len(gffs)))
            print(
        ''' 
            {:^30s}\t{:^30s}\t{:^30s}
            {:^30s}\t{:^30s}\t{:^30s}'''.format("FASTA" , "GFF" , "Nb sequences","-----","---","------------")
                )
            for n,name in enumerate(fastas):
                print(
        ''' 
            {:^30s}\t{:^30s}\t{:^30s}''' \
              .format(fastas[list(fastas.keys())[n]].split("/")[-1],
                      gffs[list(gffs.keys())[0]].split("/")[-1],
                      files_sampling[fastas[name]]))
                files_associations[fastas[list(fastas.keys())[n]]] = gffs[list(gffs.keys())[0]]


        else:
        # But if they do not have neither the same name with FASTA nor the same
        # number, sorry but I can do nothing for you! :)
            print('''
                  BEY
              ''')
            exit()
    return files_associations, files_sampling
    

def read_gff_file(gff_file):
    gff_dico = {}
    with open(gff_file,'r') as fi:
        for line in fi:
            if line.startswith("#") or line == "\n":
                continue
            try:
                gff_dico[line.split()[-1].split('ID=')[1].split(";")[0]] = line
            except:
                try:
                    gff_dico[line.split()[8].split('ID=')[1].split(";")[0]] = line
                except:
                    print("Check your last column of your GFF. There are spaces!!!")
    return gff_dico


def decide_which_color(value, nb_cols, minimum,maximum):
    step = (maximum - minimum) / nb_cols
    my_choice = round(abs(value - minimum) / step)
    my_rgb = sns.color_palette(palette="coolwarm", n_colors=nb_cols+1)[int(my_choice)]
    my_color = to_hex(my_rgb)
    return(my_color)


def change_color_in_gff_line(gff_dico, orf,value, nb_cols, minimum,maximum):
    my_line  = gff_dico[orf]
    my_color = decide_which_color(value=value,nb_cols=nb_cols, minimum=minimum,maximum=maximum)
    new_line = re.sub(pattern="color=.+", repl="color=" + my_color, string=my_line)
    new_line = new_line.strip() + ";element_value=" + str(value) + "\n"
    return new_line
    


TANGO_EXEC = {
    "darwin": "tango2_3_1",
    "win32": "Tango.exe",
    "linux": "tango_x86_64_release",
}


def main():
    start_time = datetime.now()
    parameters = get_args()

    # get path to external softwares given in config.ini
    external_softwares = orfold_utils.read_config_file()

    # Check if the tools asked for the analysis are well Installed
    if "H" in parameters.options:
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")

                from pyHCA import HCA
                from pyHCA.core.annotateHCA import _annotation_aminoacids
                from pyHCA.core.classHCA import compute_disstat
        except:
            print('''
                  Oups! pyHCA is not Installed. 
                  Please go to the link:    https://github.com/T-B-F/pyHCA 
                  and follow the installation instructions 
                  ''')
            exit()

    if "I" in parameters.options:
        is_valid, error_message = orfold_utils.check_path(path=external_softwares["iupred"], software="iupred")
        if not is_valid:
            print('''

{}

If you have installed IUPred2a on your computer, please edit the softwares.ini file
at the root of your ORFmine project and give the absolute root path of the 
parent directory where the IUPred2a library and data/ folder reside.

Otherwise, please go to this link:  https://iupred2a.elte.hu/download_new
and follow the install instructions. Then edit the config.ini file.
                  '''.format(error_message))
            exit()

        # import iupred2a_lib from its location
        spec = importlib.util.spec_from_file_location("iupred2a_lib", str(Path(external_softwares["iupred"]) / "iupred2a_lib.py"))
        iupred2a_lib = importlib.util.module_from_spec(spec)
        sys.modules["iupred2a_lib"] = iupred2a_lib
        spec.loader.exec_module(iupred2a_lib)

    if "T" in parameters.options:
        is_valid, error_message = orfold_utils.check_path(path=external_softwares["tango"], software="tango")
        if not is_valid:
            print('''

{}

If you have installed Tango on your computer, please edit the softwares.ini file
at the root of your ORFmine project and give the absolute root path of the 
parent directory where the Tango executable resides.

Otherwise, please go to this link: http://tango.crg.es/about.jsp
and follow the install instructions. Then edit the config.ini file.
                  '''.format(error_message))
            exit()


    #check_tools_installed(parameters=parameters)
    make_tmp_directories(parameters=parameters)
    files_associations , files_sampling = make_files_associations(parameters=parameters)

    for fasta_file in files_associations:
        name = os.path.basename(fasta_file)
        name = os.path.splitext(name)[0]
        size = files_sampling[fasta_file]

        # We read the MultiFasta file with the sequences:
        global sequences
        sequences = orfold_utils.read_multiFASTA(fasta_file)

        # Here we generate a sample of the initial fasta file
        if size != "all":
            indexes = random.sample(k=int(size),population=range(len(sequences)))
            indexes.sort()
            d = {list(sequences.keys())[i]:list(sequences.values())[i] for i in indexes}
            sequences = d.copy()
            del(d)


        #@TODO To be integrated later!!!
        if "H" in parameters.options and parameters.barcodes == 'True':
            # First we calculate ALL the Barcodes at ones:
            print("Calculating the HCA barcodes ... ")
            barcodes = calculate_HCA_barcodes(sequences = sequences)
            print("Barcodes calculation : DONE")
            print("\n")
            with open(name+".barcodes","w") as barcw:
                for i in barcodes:
                    barcw.write(">{}\n{}\n".format(i,barcodes[i]))
                    
        # Just some formating options for making beautiful the output file
        max_name = len(max(list(sequences.keys()), key=len))
        formating_a = "{:"+str(max_name+2)+"s}\t{:7s}\t{:7s}\t{:7s}\n"
        formating_b = "{:"+str(max_name+2)+"s}\t"

        with open(name+".tab","w") as fw_output:
            # We write the title in the output table
            fw_output.write(formating_a.format("Seq_ID","HCA","Disord","Aggreg"))


            if files_associations[fasta_file] != '':
                gff_dico = read_gff_file(gff_file = files_associations[fasta_file])
                if "H" in parameters.options:
                    fw_gff_H = open(name+'_HCA.gff','w')
                if "I" in parameters.options:
                    fw_gff_I = open(name+'_IUPRED.gff','w')
                if "T" in parameters.options:
                    fw_gff_T = open(name+'_TANGO.gff','w')

            print("\n\n")
            for i, orf in enumerate(sequences):
                seq = sequences[orf]

                if "H" in parameters.options:
                    score, pvalue = compute_disstat(0, len(seq), _annotation_aminoacids(seq=seq,method="domain",verbose=False)["cluster"])
                    score = round(score,2)

                    # If a gff file was given then we change the colour based on the HCA score
                    if files_associations[fasta_file] != '':
                        try:
                            new_gff_line = change_color_in_gff_line(gff_dico=gff_dico, orf=orf, value=score, nb_cols=20, minimum=-10, maximum=10)
                            fw_gff_H.write(new_gff_line)
                        except:
                            print('An error aoccured at the writing of the {}_HCA.gff file for the orf: {}'.format(name, orf))
                            pass
                else:
                    score = "NaN"


                if "I" in parameters.options:
                    try:
                        iupred_score  = iupred2a_lib.iupred(seq=seq ,mode="short")[0]
                        iupred_portion = orfold_utils.calculate_proportion_of_seq_disordered(iupred_score)
                        iupred_mean = round(sum(iupred_score) / len(iupred_score),2)
                    except:
                        iupred_mean, iupred_portion = "NaN","NaN"

                    # If a gff file was given then we change the colour based on the IUPRED score.
                    if files_associations[fasta_file] != '':
                        try:
                            new_gff_line = change_color_in_gff_line(gff_dico=gff_dico, orf=orf, value=iupred_portion, nb_cols=20, minimum=0, maximum=1)
                            fw_gff_I.write(new_gff_line)
                        except:
                            print('An error aoccured at the writing of the {}_IUPRED.gff file for the orf: {}'.format(name, orf))
                            pass
                else:
                    iupred_mean, iupred_portion = "NaN","NaN"


                if "T" in parameters.options:
                    tango_path = Path(external_softwares["tango"])
                    if tango_path.is_dir():
                        tango_executable = str(tango_path / TANGO_EXEC[sys.platform])
                    elif tango_path.is_file() and tango_path.name in TANGO_EXEC.values():
                        tango_executable = str(tango_path)

                    tango_portion = calculate_tango_one_sequence(tango_path=tango_executable, name=orf, to_keep=parameters.keep)

                    # If a gff file was given then we change the colour based on the Tango propensity.
                    if files_associations[fasta_file] != '':
                        try:
                            new_gff_line = change_color_in_gff_line(gff_dico=gff_dico, orf=orf, value=tango_portion, nb_cols=20, minimum=0, maximum=1)
                            fw_gff_T.write(new_gff_line)
                        except:
                            print('An error aoccured at the writing of the {}_TANGO.gff file for the orf: {}'.format(name, orf))
                            pass
                else:
                    tango_portion = "NaN"

                # --------------------------------------- #
                # We write line-by-line the table output  #
                # --------------------------------------- #

                fw_output.write(formating_b.format(orf))
                try:
                    fw_output.write("{:<7.3f}\t".format(float(score)))
                except:
                    fw_output.write("{:<7s}\t".format("NaN"))
                try:
                    fw_output.write("{:<7.3f}\t".format(float(iupred_portion)))
                except:
                    fw_output.write("{:<7s}\t".format("NaN"))
                try:
                    fw_output.write("{:<7.3f}\t".format(float(tango_portion)))
                except:
                    fw_output.write("{:<7s}\t".format("NaN"))
                #try:
                #    fw_output.write("{:s}\t".format(str(barcodes[orf])))
                #except:
                #    fw_output.write("{:<7.3f}\t".format(float("NaN")))
                fw_output.write("\n")


    try:
        fw_gff_H.close()
    except:
        pass
    # ------------------------------------------------------------------- #
    # If the plot option is True then we plot the HCA score distribution  #
    # ------------------------------------------------------------------- #
    if parameters.plot:

        tabs_to_plot = ""
        for file in files_associations:
            name = os.path.basename(file)
            name = os.path.splitext(name)[0]
            tabs_to_plot = tabs_to_plot + ' ' + name + ".tab"


        plot_command = "plot_orfold -tab {}".format(tabs_to_plot)
        subprocess.Popen(plot_command,stdout=subprocess.PIPE, stderr=None, shell=True)

    end_time = datetime.now()
    print("\n\n")
    print('Duration: {}'.format(end_time - start_time))
    exit()

    

    





