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
from orfold.lib.orfold_funcs import *
from matplotlib.colors import to_hex
import joblib

#from pyHCA import HCA
#from pyHCA.core.annotateHCA import _annotation_aminoacids
#from pyHCA.core.classHCA import compute_disstat


# -------------- # ========================================================= #
# Paths to set:  #
# -------------- # 
softwares_path="/Users/christospapadopoulos/Documents/de_novo/ORFmine/orfold_v1"
#  1. IUPRED
#  The python script of iupred
#iupred      = '/Users/christospapadopoulos/Documents/de_novo/iupred/iupred2a.py'
iupred_path  = softwares_path + '/orfold/softwares/iupred2a.py'
#  The path where is located the data file of iupred
#iupred_data = '/Users/christospapadopoulos/Documents/de_novo/iupred/'
iupred_data = softwares_path + '/orfold/softwares/'
#  2. TANGO
#  The script that excecutes Tango
#tango       = "/Users/christospapadopoulos/Documents/de_novo/TANGO/tango2_3_1"
#  If system is MacOS 
if sys.platform == "darwin":
    tango       = softwares_path + '/orfold/softwares/tango2_3_1'
#  If system is Windows
elif sys.platform == "win32":
    tango       = softwares_path + '/orfold/softwares/Tango.exe'
#  If system is Linux
elif sys.platform == "linux" or sys.platform == "linux2":
    tango       = softwares_path + '/softwares/tango_x86_64_release'
# ========================================================================== #


# ---------- # ============================================================= #
# FUNCTIONS  #
# ---------- #

def get_args():
    """
    Returns:
        Parameters
    """

    parser = argparse.ArgumentParser(description='ORF Foldability Calculation')
    parser.add_argument("-fna",
                        type=str,
                        action='store',
                        required=True, 
                        nargs="*",
                        help="Genomic fasta file(s) (.fna) ")
    
    parser.add_argument("-gff", 
                        required=False, 
                        type=str,
                        action='store',
                        nargs="*",
                        default=[],
                        help="GFF annotation file (.gff)")
    
    parser.add_argument("-options",
                        type=list,
                        #action='store',
                        required=True, 
                        nargs="?",
                        help="Which analyses you want to perform?")
    
    parser.add_argument("-plot", 
                        required=False, 
                        type=str,
                        action='store',
                        nargs="*",
                        default=False,
                        help="GFF annotation file (.gff)")
    parser.add_argument("-keep", 
                        required=False, 
                        type=list,
                        #action='store',
                        nargs="?",
                        default=[],
                        help="Option for keeping the IUPRED & TANGO output files")
    parser.add_argument("-N", 
                        required=False, 
                        type=str,
                        action='store',
                        nargs="*",
                        default=["all"],
                        help="Size of sample(s) per fasta file")
    
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
    hca = HCA(seq=list(sequences.values()),querynames=list(sequences.keys()))
    # You can get the barcode sequence of ONE sequnce using the module get_hca_barcode
    # If you want to take the barcode of all the sequnces do a loop like:
    barcodes = {}
    for orf in list(sequences.keys()):
        barcodes[orf] = get_hca_barcode(hca = hca, orf = orf)  
    return barcodes
    
def calculate_IUPRED_one_sequence(name):
    """
    Calculates the IUPRED for one sequence 
    """
    name = str(name)
    iupred_command = "python " + iupred + " -d " + iupred_data+ " -a" + ' ./FASTA_tmp/'+name+'.fasta_tmp ' + "short"
    # We create an output file for IUPRED if we have keep I option
    if "I" in parameters.keep:
        logfile = './IUPRED/'+name+'.iupred' 
    else:
        logfile = False
        
    process = subprocess.Popen(iupred_command, stderr=None, shell=True, stdout=subprocess.PIPE )   #stdout=subprocess.PIPE
    iupred_output = process.communicate()            
    #print(iupred_output)
    iupred_seq , iupred_score , anchor_score = parse_iupred_output(iupred_output = iupred_output , check_seq = sequences[name],write_file = logfile)
    
    if iupred_seq != 'Something went wrong':
        iupred_mean = round(sum(iupred_score) / len(iupred_score),2)
        iupred_portion = calculate_proportion_of_seq_disordered(iupred_score)
        #print(seq,round(sum(iupred_score) / len(iupred_score),2),calculate_proportion_of_seq_disordered(iupred_score))
        os.system("rm ./FASTA_tmp/"+name+'.fasta_tmp')
    else:
        iupred_mean    = "NaN"
        iupred_portion = "NaN"
        #print(seq , "----> PROBLEM")
    return iupred_mean , iupred_portion
    
def calculate_tango_one_sequence(tango_path  ,name):
    tf = tempfile.NamedTemporaryFile(prefix="tango")
    tmp_name = tf.name.split("/")[-1]

    tango_command = tango_path + " " + tmp_name + " ct=\"N\" nt=\"N\" ph=\"7.4\" te=\"298\" io=\"0.1\" seq=\"" + sequences[name] + "\""
    
    process = subprocess.Popen(tango_command, stdout=subprocess.PIPE, stderr=None, shell=True)
    tango_output = process.communicate()
    
       
    if str(tango_output[0]).replace("\\n","").replace("\'","").strip() == "b88, File not properly written, try writing it up again,": 
        b_aggregation = 'None'
        TANGO_portion = "NaN"
        #print(seq + "------> PROBLEM")
    else:
        try:
            aa_seq, b_aggregation, h_aggregation = read_tango_seq(tmp_name+'.txt')
            TANGO_portion = calculate_proportion_of_seq_aggregable(b_aggregation)
        except:
            # The txt file was never created :(
            TANGO_portion = "NaN"
       
    if os.path.exists(tmp_name+".txt"):
        if 'T' not in parameters.keep:
            os.system("rm "+tmp_name+".txt")
        else:
            os.system("mv "+tmp_name+".txt "+name+".txt" )
            os.system("mv "+name+".txt ./TANGO/")
    else:
        pass
    
    return TANGO_portion

# =========================================================================== #

# Start time of the program :     

#parameters = get_args()

# ------------------------------------------- #
# Create the tmp directory to save tmp files  #
# ------------------------------------------- #
def make_tmp_directories(parameters):
    if "I" in parameters.options:
        try:
            os.system("mkdir FASTA_tmp")
        except:
            pass

    if "I" in parameters.keep:
        try:
            os.system("mkdir IUPRED")
        except:
            pass

    if "T" in parameters.keep:
        try:
            os.system("mkdir TANGO")
        except:
            pass
    

# ---------------------------------------------- #
# We find the fasta and gff files associations : #
# ---------------------------------------------- #
def make_files_associations(parameters):

    fastas = get_root_name_of_files_list(files_list=parameters.fna)
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
#print(files_associations)
    
    
# -------------------------------- #   
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

def decide_which_color(value,nb_cols,minimum,maximum):
    step = (maximum - minimum) / nb_cols
    my_choice = round((value + (0-minimum)) / step) - 1
    my_rgb = sns.color_palette(palette="coolwarm",n_colors=nb_cols)[int(my_choice)]
    my_color = to_hex(my_rgb)
    return(my_color)

def change_color_in_gff_line(gff_dico,orf,value,nb_cols,minimum,maximum):
    my_line  = gff_dico[orf]
    my_color = decide_which_color(value=value,nb_cols=nb_cols,minimum=minimum,maximum=maximum)
    new_line = re.sub(pattern="color=.+",repl="color="+my_color,string=my_line)
    new_line = new_line.strip() + ";element_value=" + str(value) + "\n"
    return new_line
    

# ----------- #
# MAIN SCRIPT #
# ----------- #
def main():
    start_time = datetime.now()
    global parameters
    parameters = get_args()
# -------------------------------------------------------------- #
# Check if the tools asked for the analysis are well Installed:  #
# -------------------------------------------------------------- #
    if "H" in parameters.options:
        try:
            # import HCA functions:
            from pyHCA import HCA
            from pyHCA.core.annotateHCA import _annotation_aminoacids
            from pyHCA.core.classHCA import compute_disstat
            print('''
                  pyHCA\t:\tCHECK''')
        except:
            print('''
                  Oups! pyHCA is not Installed. 
                  Please go to the link:    https://github.com/T-B-F/pyHCA 
                  and follow the installation instructions 
                  ''')
            exit()

    if "I" in parameters.options:
        try:
            os.path.exists(iupred_path)
            global iupred_func
            import orfold.softwares.iupred_funcs as iupred_funcs
            #from orfold.softwares.iupred_funcs import *
            print('''
                  IUPRED\t:\tCHECK''')

            # from orfold.softwares.iupred_funcs import *

        except:
            print('''
                  Oups! IuPred is not Insatalled. 
                  Please go to this link:  https://iupred2a.elte.hu/download_new
                  and follow the installation instructions 
                  ''')
            exit()
        try:
            os.path.exists(iupred_data)
        except:
            print('''
                  Oups! The data file of IuPred cannot be located! 
                  ''')

    if "T" in parameters.options:
        try:
            os.path.exists(tango)
            global tempfile
            import tempfile
            print('''
                  TANGO\t:\tCHECK''')



        except:
            print('''
                  Oups! Tango is not Insatalled. 
                  Please go to this link:  http://tango.crg.es/about.jsp
                  and follow the installation instructions 
                  ''')
            exit()
# ---------------------------------------------------------------------------
    #check_tools_installed(parameters=parameters)
    make_tmp_directories(parameters=parameters)
    files_associations , files_sampling = make_files_associations(parameters=parameters)

    for fasta_file in files_associations:
        name = os.path.basename(fasta_file)
        name = os.path.splitext(name)[0]
        size = files_sampling[fasta_file]

        # We read the MultiFasta file with the sequences:
        global sequences

        sequences = read_multiFASTA(fasta_file)

        # Here we generate a sample of the initial fasta file
        if size != "all":
            indexes = random.sample(k=int(size),population=range(len(sequences)))
            indexes.sort()
            d = {list(sequences.keys())[i]:list(sequences.values())[i] for i in indexes}
            sequences = d.copy()
            del(d)

        # Just some formating options for making beautiful the output file
        max_name = len(max(list(sequences.keys()), key=len))
        formating_a = "{:"+str(max_name+2)+"s}\t{:7s}\t{:7s}\t{:7s}\n"
        formating_b = "{:"+str(max_name+2)+"s}\t"

        with open(name+".tab","w") as fw_output:
            # We write the title in the output table
            fw_output.write(formating_a.format("Seq_ID","HCA","Disord","Aggreg"))

            #@TODO To be integrated later!!!
            #if "H" in parameters.options:
                # First we calculate ALL the Barcodes at ones:
            #    barcodes = calculate_HCA_barcodes(sequences = sequences)


            if files_associations[fasta_file] != '':
                gff_dico = read_gff_file(gff_file = files_associations[fasta_file])
                if "H" in parameters.options:
                    fw_gff_H = open(name+'_HCA.gff','w')

            print("\n\n")
            for i,orf in enumerate(sequences):

                #print("\r\t0%|{:10s}|{:>3.0f}%".format(str("#"*round(round(i/len(sequences),2)*100/10)),round(i/len(sequences),2)*100),end ='')

                seq = sequences[orf]

                # ------------------------------------------------ #
                # If HCA option is on we calculate the HCA score:  #
                # ------------------------------------------------ #
                if "H" in parameters.options:
                    #score, pvalue = compute_disstat(0, len(list(sequences.values())[i]), _annotation_aminoacids(seq=list(sequences.values())[i],method="domain",verbose=False)["cluster"])
                    score, pvalue = compute_disstat(0, len(seq), _annotation_aminoacids(seq=seq,method="domain",verbose=False)["cluster"])
                    score = round(score,2)
                    # If a gff file was given then we change the colour based on
                    # the HCA score.
                    if files_associations[fasta_file] != '':
                        try:
                            new_gff_line = change_color_in_gff_line(gff_dico = gff_dico,orf=orf,value=score,nb_cols=20,minimum=-10,maximum=10)
                            fw_gff_H.write(new_gff_line)
                        except:
                            print('There is a probleme at the writing of the {}_HCA.gff file for the orf: {}'.format(name,orf))
                            pass
                else:
                    score = "NaN"

                # ----------------------------------------- #
                # If the IUPRED option is on we run IUPRED  #
                # ----------------------------------------- #
                if "I" in parameters.options:
                    #print(iupred(seq=seq,mode="short")[0])
                    # We we must write a tmp fasta file to run IUPRED
                    #with open('./FASTA_tmp/'+orf+'.fasta_tmp',"w") as fw:
                    #    fw.write(">{}\n{}\n".format(orf,sequences[orf]))

                    # We check that the fasta file was created
                    #if os.path.exists('./FASTA_tmp/'+orf+'.fasta_tmp'):
                        # And we launch the calculation of IUPRED
                    #    iupred_mean,iupred_portion = calculate_IUPRED_one_sequence(name = orf)
                    try:
                        iupred_score   = iupred_funcs.iupred(seq=seq,mode="short")[0]
                        iupred_portion = calculate_proportion_of_seq_disordered(iupred_score)
                        iupred_mean = round(sum(iupred_score) / len(iupred_score),2)
                    except:
                        iupred_mean,iupred_portion = "NaN","NaN"
                else:
                    iupred_mean,iupred_portion = "NaN","NaN"

                # ----------------------------------------- #
                # If the TANGO option is on we run TANGO    #
                # ----------------------------------------- #
                if "T" in parameters.options:
                    tango_portion = calculate_tango_one_sequence(tango_path = tango , name = orf)
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

    

    





