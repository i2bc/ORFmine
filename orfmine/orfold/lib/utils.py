#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:32:48 2020

@author: christospapadopoulos
"""
import configparser
from pathlib import Path
import random
import re
import sys
from typing import Union
import warnings

from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import seaborn as sns

from orfmine import ROOT_PATH


def error_missing_softwares(software: str):
    softwares_info = {
        "iupred": {
            "name": "IUPred2a",
            "link": "https://iupred2a.elte.hu/download_new"
        },
        "tango": {
            "name": "Tango",
            "link": "http://tango.crg.es/about.jsp"
        }
    }

    error_message_template = '''
{{}}

If you have installed {} on your computer, please edit the softwares.ini file
at the root of your ORFmine project and give the absolute root path of the 
parent directory where the IUPred2a library and data/ folder reside.

Otherwise, please go to this link:  {}
and follow the install instructions. Then edit the softwares.ini file.

'''.format("", softwares_info[software]["name"], softwares_info[software]["link"])
    
    return error_message_template



def parse_orfold_tab(tab):
    HCA = []
    IUPRED = []
    TANGO  = []
    with open(tab,"r") as f:
        for x,line in enumerate(f):
            if x!= 0:
                try:
                    HCA.append(float(line.split()[1]))
                except:
                    HCA.append(line.split()[1])
                try:
                    IUPRED.append(float(line.split()[2]))
                except:
                    IUPRED.append(line.split()[2])
                try:
                    TANGO.append(float(line.split()[3]))
                except:
                    TANGO.append(line.split()[3])
    return(HCA,IUPRED,TANGO)

    
def generate_the_plot(to_plot,colors,bins,labels):
    my_rgb = sns.color_palette(palette="coolwarm",n_colors=20,as_cmap=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        handle = []
        for i,data in enumerate(to_plot):
            if i < 1:
               sns_plot = sns.distplot(data,bins=bins[i],
                             hist=True,
                             kde = True,
                             kde_kws = {'linewidth': 0},
                             color=colors[i],label=None)
               handle.append(plt.Rectangle((0,0),1,1,color=colors[i],alpha=0.5,ec="k",linewidth=0))
            else:
               sns_plot = sns.distplot(data,bins=15,
                             hist=False,
                             kde = True,
                             kde_kws = {'linewidth': 3},
                             color=colors[i],label=labels[i-3])
               handle.append(plt.Line2D([0], [0],color=colors[i],lw=3))
        
        sns_plot = plt.xlim(-11, 11)
        sns_plot = plt.ylim(0, 0.4)
        sns_plot = plt.xlabel(xlabel="HCA score")
        #sns_plot = plt.text(x=-8,y=0.34,s="Disorder")
        #sns_plot = plt.text(x=0,y=0.34,s="Foldable")
        #sns_plot = plt.text(x=6,y=0.34,s="Aggregation")
        sns_plot = plt.legend(labels = labels,loc='upper left',handles=handle,fancybox=True, framealpha=0)
        sns_plot = plt.axvline(x = -3.319,ymin=0,ymax=5,ls = ':')
        sns_plot = plt.axvline(x =  5.214,ymin=0,ymax=5,ls = ':')

    return sns_plot
    

def read_multiFASTA(fasta_file):
    """_summary_

    @BUG: potential issue with `dico[name] + seq.replace("*","")`

    Args:
        fasta_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    dico = {}
    with open(fasta_file,'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                name = str(line.split()[0])[1:]
                dico[name] = ''
            elif line == '\n':
                continue
            else:
                seq = line.strip()
                dico[name] = dico[name] + seq.replace("*","")
    return dico

def get_seq_number(fasta_file):
    seq_number = 0
    with open(fasta_file,'r') as fasta:
        for line in fasta:
            if not line.startswith('>'):
                continue
            seq_number += 1

    return seq_number
    
    
def get_hca_barcode(hca,orf):
    '''
    This module generates the HCA barcode of the total sequence
    Clusters of <= 4 residues are neglected
    '''
    barcode = "." * len(hca.get_seqbin(orf))
    clusters = hca.get_clusters(orf)
    for x,cluster in enumerate(clusters):
        cluster_elements = str(cluster).split('\t')
        if len(cluster_elements[-1]) > 4:
            barcode = barcode[:int(cluster_elements[1])-1] + cluster_elements[-1] + barcode[int(cluster_elements[2]):]
    return barcode

  
def parse_iupred_output(iupred_output,check_seq,write_file = False):
    output = str(iupred_output[0]).split('\\n')
    
    sequence = []
    iupred_score = []
    anchor_score = []
    
    # If we gave an output file we will write the output
    if write_file:
        iupred_w = open(write_file,"w")
    
    for line in output:
        # Try to write the output of IUPRED in a file
        try:
            iupred_w.write(line.replace("\\t","\t").replace("b\'","")+"\n")
        except:
            pass
        # --------------------------------------------   
        if line.startswith("b") or line.startswith("#"):
            continue
        
        try:
            sequence.append(line.split("\\t")[1])
            iupred_score.append(float(line.split("\\t")[2]))
            anchor_score.append(float(line.split("\\t")[3]))
        except:
            continue
    # If the output file was opened here we close it
    try:
        iupred_w.close()
    except:
        pass
    # ----------------------------------------------
    if len(sequence) == len(iupred_score) and len(sequence) == len(anchor_score) and "".join(sequence) == check_seq:
        return(sequence , iupred_score , anchor_score)
    else:
        return "Something went wrong", "Something went wrong", "Something went wrong"
            
    
def calculate_proportion_of_seq_disordered(iupred_score):
    count_seg_tmp = 0
    count_agg_seg = 0
    for i,pos in enumerate(iupred_score):
        if pos > 0.5:
            count_seg_tmp += 1 
        elif pos <= 0.5 and count_seg_tmp >= 5:
            count_agg_seg = count_agg_seg + count_seg_tmp
            count_seg_tmp = 0
        elif pos <= 0.5 and count_seg_tmp < 5:
            count_seg_tmp = 0 
            continue
        
        if i == len(iupred_score) -1:
            if count_seg_tmp >= 5:
                count_agg_seg = count_agg_seg + count_seg_tmp
            else:
                continue
               
    return round(count_agg_seg/len(iupred_score), 3)

  
def read_tango_seq(file):
    aa_seq        = []
    b_aggregation   = []
    h_aggregation = []
    with open(file,'r') as f:
        for i,line in enumerate(f):
            if i == 0:
                continue
            aa_seq.append(line.split()[1])
            b_aggregation.append(float(line.split()[5]))
            h_aggregation.append(float(line.split()[6]))
            
    return(''.join(aa_seq),b_aggregation,h_aggregation)


def calculate_proportion_of_seq_aggregable(b_aggregation):
    count_seg_tmp = 0
    count_agg_seg = 0
    for i,pos in enumerate(b_aggregation):
        if pos > 5.0:
            count_seg_tmp += 1 
        elif pos <= 5.0 and count_seg_tmp >=5:
            count_agg_seg = count_agg_seg + count_seg_tmp
            count_seg_tmp = 0
        elif pos <= 5.0 and count_seg_tmp <5:
            count_seg_tmp = 0 
            continue
        
        if i == len(b_aggregation) -1:
            if count_seg_tmp >=5:
                count_agg_seg = count_agg_seg + count_seg_tmp
            else:
                continue
        
    return(round(count_agg_seg/len(b_aggregation),3))
    


def check_path(path: str="", software: str="iupred"):
    software_name = "IUPred2a" if software == "iupred" else "Tango"
    error_message = ""
    missing_files = []
    is_valid = True

    exec_sources = {
        "iupred": {
            "darwin": ["iupred2a_lib.py", "iupred2a.py", "data"],
            "win32": ["iupred2a_lib.py", "iupred2a.py", "data"],
            "linux": ["iupred2a_lib.py", "iupred2a.py", "data"],
        },
        "tango": {
            "darwin": ["tango2_3_1"],
            "win32": ["Tango.exe"],
            "linux": ["tango_x86_64_release"],
        }
    }

    # case where the path given by the user do not exist
    if not Path(path).exists():
        error_message = "Sorry, it looks like '{}' do not exist... Please ensure to give the correct absolute path where {} resides.".format(path, software_name)
        is_valid = False

        return is_valid, error_message

    if not path:
        error_message = """Sorry, no path has been given in {}...
        """.format(str(ROOT_PATH / 'softwares.ini'))
        is_valid = False

        return is_valid, error_message

    # get expected source files according to the software and platform
    if sys.platform not in exec_sources["iupred"]:
        source_files = exec_sources[software]["linux"]
    else:
        source_files = exec_sources[software][sys.platform]

    # cases where the given path exists but expected at least one expected source files is missing 
    for filename in source_files:
        source_file = Path(path) / filename
        if not source_file.exists():
            missing_files.append(str(source_file))
            is_valid = False

    if missing_files:        
        term = "files" if len(missing_files) > 1 else "file"        
        error_message = "Sorry, it looks like the following {} do not exist in {}: {}".format(term, path, ", ".join(missing_files))

    return is_valid, error_message


def get_root_name_of_files_list(files_list):
    """
    Removes the extentions and path of the files
    """
    return {Path(i).stem: i for i in files_list}


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
                gff_dico[line.split('ID=')[-1].split(";")[0]] = line
            except:
                try:
                    gff_dico[line.split()[8].split('ID=')[1].split(";")[0]] = line
                except:
                    print("Check your last column of your GFF. There are spaces!!!")
    return gff_dico


def decide_which_color(value, nb_cols, minimum, maximum):
    step = (maximum - minimum) / nb_cols
    my_choice = round(abs(value - minimum) / step)
    my_rgb = sns.color_palette(palette="coolwarm", n_colors=nb_cols+1)[int(my_choice)]
    my_color = to_hex(my_rgb)

    return my_color


def change_color_in_gff_line(my_line, value, minimum, maximum, nb_cols=20):
    my_color = decide_which_color(value=value, nb_cols=nb_cols, minimum=minimum, maximum=maximum)
    new_line = re.sub(pattern="color=.+", repl="color=" + my_color, string=my_line)
    new_line = new_line.strip() + ";element_value=" + str(value) + "\n"

    return new_line


def set_gffline_color(line_template, value, minimum, maximum, nb_cols=20):
    color = decide_which_color(value=value, nb_cols=nb_cols, minimum=minimum, maximum=maximum)

    # Check if 'color' tag exists in the line_template
    if "color=" in line_template:
        # new_line = re.sub(pattern="color=.+", repl="color=" + color, string=line_template)
        new_line = re.sub(pattern=r"color=[^;]*", repl="color=" + color, string=line_template)
    else:
        # If not, add the color tag. Add it before the potential last semicolon.
        if line_template.endswith(";"):
            new_line = line_template[:-1] + ";color=" + color
        else:
            new_line = line_template + ";color=" + color

    new_line += ";element_value=" + str(value) + "\n"

    return new_line



def get_orfold_out_format(max_len_head: int=12):
    return "{:"+str(max_len_head+2)+"s}\t{:7s}\t{:7s}\t{:7s}\n"


def get_sequences(fasta_file: str, sample_size: Union[str, int] = "all"):
    sequences = read_multiFASTA(fasta_file)

    if sample_size != -1:
        indexes = sorted(random.sample(k=int(sample_size), population=range(len(sequences))))
        sequences = {list(sequences.keys())[i]: list(sequences.values())[i] for i in indexes}

    if not sequences:
        raise ValueError(f"Failed to obtain sequences from the provided FASTA file: {Path(fasta_file.name)}")
        
    return sequences


def format_with_n_decimals(number, n=3):
    return f"{number:.{n}f}"


def fasta_generator(filename):
    header = None
    sequence = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:  # This indicates we're moving to a new sequence.
                    sequence[-1] = sequence[-1].rstrip("*")
                    yield header.lstrip(">"), ''.join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header:  # Yield the last sequence
            sequence[-1] = sequence[-1].rstrip("*")
            yield header.lstrip(">"), ''.join(sequence)


def iter_to_generator(iterable):
    """Make an iterable a generator"""
    for item in iterable:
        yield item


def reservoir_sampling_fasta(filename, k):
    """Sample k sequences from a FASTA file using reservoir sampling."""
    reservoir = []
    current_header = None
    current_sequence = []

    with open(filename, 'r') as file:
        for index, line in enumerate(file):
            line = line.strip()

            # If line is a header (i.e., a new sequence)
            if line.startswith('>'):
                # Check if there's a sequence stored and try to add it to the reservoir
                if current_header and current_sequence:
                    current_sequence[-1] = current_sequence[-1].rstrip('*')
                    sequence_entry = (current_header.lstrip(">"), ''.join(current_sequence))

                    if len(reservoir) < k:
                        reservoir.append(sequence_entry)
                    else:
                        rand_idx = random.randint(0, index)
                        if rand_idx < k:
                            reservoir[rand_idx] = sequence_entry

                # Reset for next sequence
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

        # Handle the last sequence in the file
        if current_header and current_sequence:
            current_sequence[-1] = current_sequence[-1].rstrip('*')
            sequence_entry = (current_header.lstrip(">"), ''.join(current_sequence))

            if len(reservoir) < k:
                reservoir.append(sequence_entry)
            else:
                rand_idx = random.randint(0, index)
                if rand_idx < k:
                    reservoir[rand_idx] = sequence_entry

    return iter_to_generator(reservoir)
