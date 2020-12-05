#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 15:33:25 2020

@author: christospapadopoulos
"""

import argparse,os,random
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from Bio import SeqIO
import re
import linecache



def get_args():
    """
    Returns:
        Parameters
    """
    parser = argparse.ArgumentParser(description='ORF Foldability Calculation')
    parser.add_argument("-fna",
                        type=str,
                        #action='store',
                        required=True, 
                        nargs="?",
                        help="Genomic fasta file(s) (.fna) ")
    parser.add_argument("-gff",
                        type=str,
                        #action='store',
                        required=True, 
                        nargs="?",
                        help="Annotation file(s) (.gff) ")
    parser.add_argument("-o",
                        type=str,
                        #action='store',
                        required=True, 
                        nargs="?",
                        help="The output name of your fasta")
    parser.add_argument("-elements_include",
                        type=str,
                        action='store',
                        required=False, 
                        nargs="*",
                        default = ['all'],
                        help="Annotation elmements to be considered (By definition is all)")
    parser.add_argument("-elements_exclude",
                        type=str,
                        action='store',
                        required=False, 
                        nargs="*",
                        default = ["None"],
                        help="Annotation elmements not to be considered (By definition is None)")
    parser.add_argument("-chr_exclude",
                        type=str,
                        action='store',
                        required=False, 
                        nargs="?",
                        default = [],
                        help="Exclude chromosomes that do not interest you (mitochondrial for example)")
    parser.add_argument("-N",
                        type=int,
                        #action='store',
                        required=False, 
                        nargs="?",
                        default = False,
                        help="The size of the sample to be generated")

    args = parser.parse_args()
    return args


def read_genome(genome_fasta):
    with open(genome_fasta, 'r') as fasta_file:
        temp = ''.join([line.rstrip() if not line.startswith('>') else '\n' + line for line in fasta_file.readlines()]).split('\n')
        genome = {temp[i].replace('>', '').split()[0] : temp[i+1] for i in range(1, len(temp), 2)}
    return(genome)
    
    

def read_gff_info(gff,elements_in,elements_out):
    elements_in = "(" + ")|(".join(elements_in) + ")"
    elements_out = "(" + ")|(".join(elements_out) + ")"
    dico_info = {}
    with open(gff,'r') as f:
        print('\n')
        count = 0
        for x,line in enumerate(f):
            # If the line is not a comment
            if line.startswith('#') == False:
                # If we finished with all the sample we stop the loop
#                if N != False and len(N) == 0:
#                    break
                # If i dont want a sample but all the gff file
                # OR if i have a sample and the line is in this sample
#                if N == False or N != False and x in N:
                    # Remove the line number from the sample
#                    try:
#                        N.remove(x)
#                    except:
#                        pass
                    
                # And now read the line and extract the sequence
                element= line.split()[2].rstrip()
                
                if re.match(elements_out,element) and not re.match(elements_in,element):
                    continue
                
                if re.match(elements_out,element) and re.match(elements_in,element):
                    print("{} include and exclude at the same time!".format(element))
                    exit()
                
                if not re.match(elements_in,element) and elements_in != "(all)":
                    continue
                
                strand = line.split()[6].rstrip()
                start  = int(line.split()[3])
                stop   = int(line.split()[4])
                chrom  = line.split()[0].rstrip()
                if chrom in parameters.chr_exclude:
                    continue
                info   = line.split()[8]
                frame  = line.split()[7].rstrip()
                gene   = info.split(';')[0].split('=')[1].rstrip()
                
                if gene not in dico_info.keys():
                    count += 1
                    print('\r\t' + 'GFF' + '\t:\t' + str(count)+'\tsequences read', end = '')
                    dico_info[gene] = {}
                    dico_info[gene]['positions']=[]
                    dico_info[gene]['positions'].append(start)
                    dico_info[gene]['positions'].append(stop)
                    dico_info[gene]['DNA_seq'] = ''
                    dico_info[gene]["DNA_parts"] = []
                    gene_seen = False
                else:
                    gene_seen = True
                    
                dico_info[gene]['strand'] = strand
                dico_info[gene]['chrom']  = chrom
                dico_info[gene]['info']   = info
                dico_info[gene]['frame']  = frame
     
               
                if strand == '+':
                    
                    seq_dna = genome[chrom][start-1:stop]
                    #my_seq_dna = Seq(seq_dna.upper(), IUPAC.unambiguous_dna)
                    my_seq_dna = Seq(seq_dna.upper())
    
                    
                    if gene_seen == False:
                        dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                    elif gene_seen == True:
                        
                        if stop < min(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq'] 
    
                        elif start > max(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                        else:
                            print(gene)
                            continue

                    
                    dico_info[gene]['positions'].append(start)
                    dico_info[gene]['positions'].append(stop)
                    
                    try:
                        #dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'], IUPAC.unambiguous_dna)))#.replace('*','')
                        dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')
    
                    except:
                        dico_info[gene]['AA_seq']  = ''
                
                elif strand == '-':
                    seq_dna = genome[chrom][start-1:stop]
                    #my_seq_dna = Seq(seq_dna.upper(), IUPAC.unambiguous_dna)
                    my_seq_dna = Seq(seq_dna.upper())
                    my_seq_dna = my_seq_dna.reverse_complement()
                    #######dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                    dico_info[gene]["DNA_parts"].append(my_seq_dna)
                    
                    if gene_seen == False:
                        dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                    elif gene_seen == True:
                        if stop < min(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
    
                        elif start > max(dico_info[gene]['positions']):
                            dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq'] 
    
                        else:
                            print(gene)
                            continue

                    
                        dico_info[gene]['positions'].append(start)
                        dico_info[gene]['positions'].append(stop)
    
                    
                    try:
    #                        dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'], IUPAC.unambiguous_dna)))#.replace('*','')
                        dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')
    
                    except:
                        dico_info[gene]['AA_seq']  = ''
                            

                    
    return(dico_info)


def write_multifastas(dico_info,outname):
    names = []
    #isoforms = ["B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
#    type_of_data = "CDS"  # JUST TO NOT BUG - I WILL FIX IT
    with open(outname+'.nmultifasta','w') as fwn,open(outname+'.pmultifasta','w') as fwp:
        for gene in dico_info:
            # Name the sequences:
#            if type_of_data == 'CDS':
#                name = str(dico_info[gene]['chrom'])  +'_'+\
#                       str(dico_info[gene]['strand']) +'_'+\
#                       str(min(dico_info[gene]['positions']))+'-'+\
#                       str(max(dico_info[gene]['positions']))
#            elif type_of_data == 'IGORF':
#                if dico_info[gene]['strand'] == '+':
#                    name = str(dico_info[gene]['chrom'])  +'_'+\
#                       str(dico_info[gene]['strand']) +'_'+\
#                       str(min(dico_info[gene]['positions']))+'-'+\
#                       str(max(dico_info[gene]['positions']))+'_'+\
#                       str(dico_info[gene]['frame'])
#                elif dico_info[gene]['strand'] == '-':
#                    name = str(dico_info[gene]['chrom'])  +'_'+\
#                       str(dico_info[gene]['strand']) +'_'+\
#                       str(max(dico_info[gene]['positions']))+'-'+\
#                       str(min(dico_info[gene]['positions']))+'_'+\
#                       str(dico_info[gene]['frame'])
                       

            
#            if name not in names:
#                isoform_count = 0
#                names.append(name)
#            elif name in names:
#                isoform_count += 1
#                name = name + '_iso' +str(isoform_count)
                
            # -------------------------------------------------------- #
            
            # Here is a check point for sequences that have * in the middle
            aa_seq = dico_info[gene]['AA_seq']
            aa_seq = aa_seq.replace("*","X")
            if aa_seq[-1] == "X":
                aa_seq = aa_seq[:-1] + "*"
            else:
                # Here we do not keep genes that did not stop to STOP coddon
                # Probably the sequence stops abruptly
                continue
            
            if dico_info[gene]['DNA_seq'] != '' and dico_info[gene]['DNA_seq'][-3:] in ["TAG","TGA","TAA"]:
                fwn.write('>'+gene+'\n')
                #fwn.write('>'+gene+'\n')
                fwn.write(dico_info[gene]['DNA_seq']+'\n')
            
            if dico_info[gene]['AA_seq'] != '' and dico_info[gene]['AA_seq'][-1] == "*":
                fwp.write('>'+gene+'\n')
                #fwp.write(dico_info[gene]['AA_seq']+'\n')
                fwp.write(aa_seq+'\n')
                

def find_matches(gff , elements_in , elements_out):
    elements_in = "(" + ")|(".join(elements_in) + ")"
    elements_out = "(" + ")|(".join(elements_out) + ")"
    matches = []
    with open(gff,'r') as f:
        #count = 0
        for x,line in enumerate(f):
            if line.startswith('#') == False:
                element= line.split()[2].rstrip()
                
                if re.match(elements_out,element) and not re.match(elements_in,element):
                    continue
                    
                if re.match(elements_out,element) and re.match(elements_in,element):
                    print("{} include and exclude at the same time!".format(element))
                    exit()
                    
                if not re.match(elements_in,element) and elements_in != "(all)":
                    continue
                
                matches.append(x)
   
    return matches



def read_matched_gff_lines(matches , gff_file , genome):
    count = 0
    dico_info = {}
    for line in matches:
        my_line = linecache.getline(gff_file,line)
        my_line = my_line.split("\t")
        strand = my_line[6].rstrip()
        start  = int(my_line[3])
        stop   = int(my_line[4])
        chrom  = my_line[0].rstrip()
        info   = my_line[8]
        gene   = info.split(';')[0].split('=')[1].rstrip()
        
        if gene not in dico_info.keys():
            count += 1
            print('\r\t' + 'GFF' + '\t:\t' + str(count)+'\tsequences read', end = '')
            dico_info[gene] = {}
            dico_info[gene]['positions']=[]
            dico_info[gene]['positions'].append(start)
            dico_info[gene]['positions'].append(stop)
            dico_info[gene]['DNA_seq'] = ''
            dico_info[gene]["DNA_parts"] = []
            gene_seen = False
        else:
            gene_seen = True
            
        dico_info[gene]['strand'] = strand
        dico_info[gene]['chrom']  = chrom
        dico_info[gene]['info']   = info
        #dico_info[gene]['frame']  = frame
 
       
        if strand == '+':
            
            seq_dna = genome[chrom][start-1:stop]
            #my_seq_dna = Seq(seq_dna.upper(), IUPAC.unambiguous_dna)
            my_seq_dna = Seq(seq_dna.upper())

            
            if gene_seen == False:
                dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
            elif gene_seen == True:
                
                if stop < min(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq'] 

                elif start > max(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
                else:
                    print(gene)
                    continue

            
            dico_info[gene]['positions'].append(start)
            dico_info[gene]['positions'].append(stop)
            
            try:
                #dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'], IUPAC.unambiguous_dna)))#.replace('*','')
                dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')

            except:
                dico_info[gene]['AA_seq']  = ''
        
        elif strand == '-':
            seq_dna = genome[chrom][start-1:stop]
            #my_seq_dna = Seq(seq_dna.upper(), IUPAC.unambiguous_dna)
            my_seq_dna = Seq(seq_dna.upper())
            my_seq_dna = my_seq_dna.reverse_complement()
            #######dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
            dico_info[gene]["DNA_parts"].append(my_seq_dna)
            
            if gene_seen == False:
                dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)
            elif gene_seen == True:
                if stop < min(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = dico_info[gene]['DNA_seq'] + str(my_seq_dna)

                elif start > max(dico_info[gene]['positions']):
                    dico_info[gene]['DNA_seq'] = str(my_seq_dna) + dico_info[gene]['DNA_seq'] 

                else:
                    print(gene)
                    continue

            
                dico_info[gene]['positions'].append(start)
                dico_info[gene]['positions'].append(stop)

            
            try:
#                        dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'], IUPAC.unambiguous_dna)))#.replace('*','')
                dico_info[gene]['AA_seq']  = str(Seq.translate(Seq(dico_info[gene]['DNA_seq'])))#.replace('*','')

            except:
                dico_info[gene]['AA_seq']  = ''

    return(dico_info)
        
    

def main():
    global parameters
    parameters = get_args()

    print('''
          Genome file              : {}
          GFF file                 : {}
          Features to keep         : {}
          Features to exclude      : {}
          Chromosome to exclude    : {}
          Sample of population     : {}
          ''' . format(parameters.fna , parameters.gff , parameters.elements_include , parameters.elements_exclude , parameters.chr_exclude ,parameters.N))

    genome_file = parameters.fna
    global genome
    genome  = read_genome(genome_fasta = genome_file)

    elements_in = parameters.elements_include
    # The output type of the sequences name (With or without the frame in the end)
    #if "CDS" in elements:
    #    type_of_data = "CDS"
    #else:
    #    type_of_data = "IGORF"

    chomosomes_exclude = parameters.chr_exclude

    gff_file = parameters.gff

    if parameters.N != False:
        print("Reading sequences and generating sample of {}".format(parameters.N))
        matches = find_matches(gff = gff_file,
                               elements_in =parameters.elements_include,
                               elements_out=parameters.elements_exclude )

        matches = sorted(random.sample(k=parameters.N , population= matches))

        seqs = read_matched_gff_lines(matches , gff_file = gff_file , genome = genome)

    else:



        seqs = read_gff_info(gff=gff_file,
                             elements_in=parameters.elements_include,
                             elements_out=parameters.elements_exclude)

    #outname = os.path.basename(gff_file)
    #outname = os.path.splitext(outname)[0]
    #outname = outname + "_" + type_of_data
    #outname  = parameters.o


    write_multifastas(dico_info = seqs , outname=parameters.o)
