#!/bin/miniconda3/envs/ORFmine_env/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 19:13:53 2020

@author: christospapadopoulos
"""


import os

iupred_files_path = "/ORFmine/orfold_v1/orfold/softwares/"

with open(iupred_files_path + "iupred2a.py","r") as f, open(iupred_files_path + "iupred_funcs.py","w") as fw:

    fw.write("#!/bin/miniconda3/envs/ORFmine_env/bin/python3")

    in_def = False
    count_import = 0
    for x,line in enumerate(f):
        line2 = line.strip()

        if line2.startswith("import"):
            count_import+= 1
            #print(x,in_def , line)
            fw.write(line)
            if count_import == 4:
                fw.write("\nPATH=\"{}\"".format(iupred_files_path))

        elif line2 == "":
            #print(in_def ,line)
            fw.write(line)

        elif line2.startswith("def"):
            #print(x,in_def , line)
            fw.write(line)
            in_def = True
            continue


        elif line2.startswith("return"):
            in_def = False
            #print(x,in_def , line)
            fw.write(line)
            continue

        elif in_def == True:
           # print(x,in_def , line)
            fw.write(line)

        else:
            continue
