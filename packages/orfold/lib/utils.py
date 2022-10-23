#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:32:48 2020

@author: christospapadopoulos
"""

import matplotlib.pyplot as plt
import seaborn as sns
import random
import warnings
import configparser
from pathlib import Path
import sys

from packages.config import ROOT_PATH


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

    return(sns_plot)
    

def read_multiFASTA(fasta_file):
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
    return(dico) 
    
    
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
    return(barcode) 

  
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
        return("Something went wrong","Something went wrong","Something went wrong")
            
    
def calculate_proportion_of_seq_disordered(iupred_score):
    count_seg_tmp = 0
    count_agg_seg = 0
    for i,pos in enumerate(iupred_score):
        if pos > 0.5:
            count_seg_tmp += 1 
        elif pos <= 0.5 and count_seg_tmp >=5:
            count_agg_seg = count_agg_seg + count_seg_tmp
            count_seg_tmp = 0
        elif pos <= 0.5 and count_seg_tmp <5:
            count_seg_tmp = 0 
            continue
        
        if i == len(iupred_score) -1:
            if count_seg_tmp >=5:
                count_agg_seg = count_agg_seg + count_seg_tmp
            else:
                continue
               
    return(round(count_agg_seg/len(iupred_score),3))   

  

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
    

#Read config.ini file
def read_config_file():
    config_file = ROOT_PATH / 'softwares.ini'
    config = configparser.ConfigParser()
    config.read(str(config_file))

    for key, value in config["EXTERNAL_SOFTWARE"].items():
        config["EXTERNAL_SOFTWARE"][key] = value.strip('"')

    return config["EXTERNAL_SOFTWARE"]


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

    
