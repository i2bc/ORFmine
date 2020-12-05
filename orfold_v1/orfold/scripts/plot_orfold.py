#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 11:25:53 2020

@author: christospapadopoulos
"""

import argparse,os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import random
from orfold.lib.orfold_funcs import parse_orfold_tab,generate_the_plot


def get_args():
    """
    Returns:
        Parameters
    """
    parser = argparse.ArgumentParser(description='ORF Foldability Calculation')
    parser.add_argument("-tab",
                        type=str,
                        action='store',
                        required=True, 
                        nargs="*",
                        help="Tables of foldability calculated by orfold")
    parser.add_argument("-names",
                        type=str,
                        action='store',
                        required=False, 
                        nargs="*",
                        help="Names of datasets for the legend")
    
    args = parser.parse_args()
    return args



#dis_ref_hca , dis_ref_iupred , dis_ref_tango = parse_orfold_tab(tab = '/Users/christospapadopoulos/Documents/de_novo/Project_Nicolas/disprot_v7_protein_predictors_minsize30.tab')          
glo_ref_hca , glo_ref_iupred , glo_ref_tango = parse_orfold_tab(tab = '/Users/christospapadopoulos/Documents/de_novo/Project_Nicolas/globular.tab')          
#tra_ref_hca , tra_ref_iupred , tra_ref_tango = parse_orfold_tab(tab = '/Users/christospapadopoulos/Documents/de_novo/Project_Nicolas/Transmembrane_helices_20_nonreduntant.tab')          
to_plot = [glo_ref_hca]
#to_plot = [dis_ref_hca , glo_ref_hca ,tra_ref_hca]
colors  = ['grey']
#colors  = ['forestgreen','grey','mediumvioletred']
bins    = [17]
#bins     = [10,17,20]
labels  = ["Foldable"]
#labels  = ["Disorder","Foldable","Aggregation"]


def main():
    parameters = get_args()

    for x,file in enumerate(parameters.tab):
        name = os.path.basename(file)
        name = os.path.splitext(name)[0]
        hca , iupred , tango = parse_orfold_tab(tab = file)
        to_plot.append(hca)
        #colors.append("blue")
        colors.append([random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)])
        try:
            labels.append(parameters.names[x])
        except:
            labels.append(name)
        bins.append(15)



    my_plot = generate_the_plot(to_plot = to_plot , colors = colors , bins = bins , labels = labels)

    # Add the statistical KS-test
    try:
        from scipy.stats import ks_2samp
        KS = []
        for j in range(len(to_plot)):
            KS_tmp = []
            if j > 0:
                for t in range(1,1000):
                    my_sample_values = []
                    my_sample_idx    = random.sample(k=559,population=range(len(to_plot[j])))
                    for l in my_sample_idx :
                        my_sample_values.append(to_plot[j][l])

                    ks_test = ks_2samp(to_plot[0], my_sample_values)
                    KS_tmp.append(ks_test[1])

                mean_KS = sum(KS_tmp)/len(KS_tmp)
                if mean_KS > 0.05:
                    KS.append("")
                elif mean_KS <= 0.05 and mean_KS > 0.01:
                    KS.append("*")
                elif mean_KS <= 0.01 and mean_KS > 0.001:
                    KS.append("**")
                elif mean_KS <= 0.001:
                    KS.append("***")
                else:
                    KS.append("PROBLEM")

        my_plot = plt.text(x=-9.8,y=0.37,s=str('\n'.join(KS)),
                           verticalalignment='top',horizontalalignment= 'center',
                           fontsize = "large",linespacing = 1.35)
    except:
        pass



    fig = my_plot.get_figure()
    fig.savefig('output.png',transparent=False)
    fig.savefig('output_transparent.png',transparent=True)
    fig.savefig('output.pdf',transparent=True,dpi=300)

    fig2, ax = plt.subplots(figsize=(6, 1))
    fig2.subplots_adjust(bottom=0.5)
    norm = mpl.colors.Normalize(vmin=-10, vmax=10)
    cb1      = mpl.colorbar.ColorbarBase(ax,
                                         cmap=mpl.cm.coolwarm,
                                         orientation='horizontal',
                                          norm=norm)
                                         #boundaries=[-10,-5,0,5,10])
    fig2.savefig('Scale.pdf',transparent=True,dpi=300)

