#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 00:40:27 2022

@author: christospapadopoulos
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import time


def get_args():
    """
    Returns:
        All the Parameters given at the terminal
    """
    parser = argparse.ArgumentParser(description='ORFs sequence extractor')
    parser.add_argument("-fna",
                        type=str,
                        #action='store',
                        required=True,
                        nargs="?",
                        help="Genomic fasta file (.fna) ")
    parser.add_argument("-gff",
                        type=str,
                        #action='store',
                        required=True,
                        nargs="?",
                        help="Annotation file (.gff) ")
    parser.add_argument("-o",
                        type=str,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = ['/database/amino_acid_sequence'],
                        help="Output name of fasta file(s)")
    parser.add_argument("-features_include",
                        type=str,
                        action='store',
                        required=False,
                        nargs="*",
                        default = ['all'],
                        help="Annotation features to be considered (By definition is all)")
    parser.add_argument("-features_exclude",
                        type=str,
                        action='store',
                        required=False,
                        nargs="*",
                        default = ["None"],
                        help="Annotation features not to be considered (By definition is None)")
    parser.add_argument("-chr_exclude",
                        type=str,
                        action='store',
                        required=False,
                        nargs="?",
                        default = [],
                        help="Chromosomes to be excluded (By definition is None)")
    parser.add_argument("-N",
                        type=int,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = False,
                        help="Size of the sample to be generated")
    parser.add_argument("-type",
                        type=str,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = "prot",
                        help="Sequences type of output [prot (default) / nucl / both];")
    parser.add_argument("-table",
                        type=str,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = "Standard",
                        help='Which codon table to use?  This can be either a name (string), an NCBI identifier (integer), or a CodonTable  object (useful for non-standard genetic codes). This defaults to the "Standard" table.')
    parser.add_argument("-check",
                        #type=bool,
                        action='store_true',
                        required=False,
                        #nargs="?",
                        default = False,
                        help='DO NOT write sequences that DO NOT finish with a STOP codon')
    parser.add_argument("-elongate",
                        type=int,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = False,
                        help="Number of nucletides to elongate vers 5 & 3 UTR")
    parser.add_argument("-name_attribute",
                        type=str,
                        #action='store',
                        required=False,
                        nargs="?",
                        default = "Name",
                        help="GFF 'Name' attribute")

    args = parser.parse_args()
    return args


class GFF_element:
    '''

    '''
    def __init__(self,gff_line,genome=None,name_attribute="Name",genetic_code=1):
        self.chromosome = gff_line.split("\t")[0]
        self.info       = self.__organize_info__(gff_line.split("\t")[-1])
        self.identity   = self.__get_feature_identity__(self.info,name_attribute)
        self.phase      = gff_line.split("\t")[-2]
        self.strand     = gff_line.split("\t")[-3]
        self.end        = int(gff_line.split("\t")[-5])
        self.start      = int(gff_line.split("\t")[-6])
        self.feature    = gff_line.split("\t")[-7]
        if genome:
            self.seq_nucl = self.__get_nucletide_seq__(self.chromosome,self.start,self.end,genome)
            self.seq_nucl_elongated = None
            self.UTR5_start = None
            self.UTR5_end   = None
            self.UTR3_start = None
            self.UTR3_end   = None


    def __organize_info__(self,info):
        dico = {}
        for i in info.strip().split(";"):
            try:
                dico[i.split("=")[0]] = i.split("=")[1]
            except:
                dico[i] = None
        #return({i.split("=")[0]:i.split("=")[1] for i in info.split(";")})
        return dico

    def __get_feature_identity__(self,info_organized,name_attribute):
        if "ID" in info_organized:
            return(info_organized["ID"])
        elif "Name" in info_organized:
            return(info_organized["Name"])
        elif "Parent" in info_organized:
            return(info_organized["Parent"])
        elif name_attribute in info_organized:
            return(info_organized[name_attribute])
        else:
            return(None)

    def __get_nucletide_seq__(self,chromosome,start,end,genome):
        return(genome[chromosome][start-1:end])

    def __update_feature__(self,gff_element):
        '''
        Updates the nucleotide sequence of the feature by concatenating it
        with the nucleotide sequence of another gff_element when they share the same identity.
        This function concatenates the exons of the CDS for example.
        '''
        if gff_element.end < self.start:
            self.seq_nucl = gff_element.seq_nucl + self.seq_nucl
            self.start    = gff_element.start
        elif gff_element.start > self.end:
            self.seq_nucl = self.seq_nucl + gff_element.seq_nucl
            self.end      = gff_element.end

    def __correct_feature_sequences__(self,genetic_code=1):
        '''
        Makes REV-COM for the sequences in strand - and translates them into
        amino acids based on the genetic code asked.
        All the sequences are given back 5` --> 3`
        '''
        if self.strand == "+":
            self.seq_prot = self.seq_nucl.translate(table=genetic_code)

        elif self.strand == "-":
            self.seq_nucl = self.seq_nucl.reverse_complement()
            self.seq_prot = self.seq_nucl.translate(table=genetic_code)

    def __elongate__(self,genome,elongate = 50):
        '''
        Elongates the nucleotide sequence towards the 5 and 3 UTRs
        '''
        # We get the nucl seq of the elongations
        elongate_5UTR_nucl = self.__get_nucletide_seq__(self.chromosome,(self.start - elongate),(self.start-1),genome)
        elongate_3UTR_nucl = self.__get_nucletide_seq__(self.chromosome,(self.end+1),(self.end + elongate),genome)
        # We construct the elongated nucl seq
        if self.strand == "-":
            self.seq_nucl_elongated = elongate_3UTR_nucl.reverse_complement() + self.seq_nucl + elongate_5UTR_nucl.reverse_complement()
        else:
            self.seq_nucl_elongated = elongate_5UTR_nucl + self.seq_nucl + elongate_3UTR_nucl
        ####self.seq_nucl_elongated = elongate_5UTR_nucl + self.seq_nucl + elongate_3UTR_nucl
        ###self.seq_nucl_elongated = self.__get_nucletide_seq__(self.chromosome,(self.start - elongate),(self.end + elongate),genome)
        # If is in the - strand we get the REV-COMP of the sequence

        # And we keep the indexes of the elongations
        self.UTR5_start = 1
        self.UTR5_end   = int(elongate)
        self.UTR3_start = int(len(self.seq_nucl_elongated) - elongate)
        self.UTR3_end   = int(len(self.seq_nucl_elongated))


class GFF_iterator:
    '''
    This class parses a GFF3 file into multiple features
    '''
    def __init__(self,gff_file,outname,output_type,chr_exclude=[None],gff_types=["all"],genome=None,genetic_code=1,check=True,elongate=False,name_attribute="Name"):
        self.file         = gff_file
        self.genetic_code = genetic_code
        self.chr_exclude  = chr_exclude
        self.gff_types    = gff_types
        self.outname      = outname
        self.output_type  = output_type
        self.name_attribute = name_attribute
        self.features     = self.__get_features__(gff_file,gff_types,genome,chr_exclude,genetic_code,outname,output_type,check,elongate,name_attribute)

    def __get_features__(self,gff_file,gff_types,genome,chr_exclude,genetic_code,outname,output_type,check,elongate,name_attribute):
        # Initialize the list of features
        features = []
        # Initialize the variable previous_feature to test if this feature.identity exists already and update it
        previous_feature = ""

        if output_type == "prot" or output_type == "both":
            fwp = open(outname + '.pfasta', 'w')

        if output_type == "nucl" or output_type == "both":
            fwn = open(outname + '.nfasta', 'w')

        if elongate != False:
            fasta_elongate = open(outname + '_elongated.nfasta','w')
            gff_elongate   = open(outname + '_elongated.gff','w')

        # Loop in every line of the GFF file:
        with open(gff_file,"r") as gff:
            for line in gff:
                if not line.startswith("#"):
                    # For every line we create a feature class
                    feature = GFF_element(line,genome,name_attribute)

                    # We check if the filter is in the list of features wanted
                    if feature.feature in gff_types and feature.chromosome not in chr_exclude or gff_types == ["all"]:
                        # If the features list is empty (case of the first feature)
                        # then we add the feature without making any other action
                        if len(features) == 0:
                            previous_feature = feature.identity
                            features.append(feature)
                            continue
                        # If the feature is not already in the list we add it (I check only the last feature)
                        if feature.identity != previous_feature:
                            previous_feature = feature.identity
                            # We append the new feature in the list
                            features.append(feature)

                            # We translate the previous feature into AA
                            features[0].__correct_feature_sequences__(genetic_code=genetic_code)

                            # If the check option is active then we check if the
                            # protein sequence ends to a STOP codon (*). If NO we
                            # skip this feature and do not write it in the output:
                            if not features[0].seq_prot.seq.endswith("*") and check == True:
                                del features[0]
                                continue


                            # If the elongation option is activate, we elongate the feture sequence:
                            if elongate != False:
                                # We elongate the sequence of the feature
                                features[0].__elongate__(genome=genome,elongate=elongate)

                                fasta_elongate.write(">{}\n{}\n".format(features[0].identity+"_mRNA",str(features[0].seq_nucl_elongated.seq)))

                                gene_end  = len(features[0].seq_nucl_elongated.seq)
                                cds_start = elongate+1
                                cds_end   = elongate+len(str(features[0].seq_nucl.seq))

                                # We write the GENE feature in the GFF : ID=identity
                                gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","gene","1",gene_end,".","+",".","ID=" + features[0].identity))
                                # We write the mRNA feature in the GFF : ID=identity_mRNA
                                gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","mRNA","1",gene_end,".","+",".","ID="+features[0].identity+"_mRNA;Parent="+features[0].identity))
                                # We write the CDS feature in the GFF  : ID=identity_CDS
                                gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","CDS",cds_start,cds_end,".","+",".","ID=" + features[0].identity + "_CDS;Parent="+features[0].identity+"_mRNA"))
                                # We write the 5UTR feature in the GFF  : ID=identity_5UTR
                                gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","five_prime_UTR",features[0].UTR5_start,features[0].UTR5_end,".","+",".","ID=" + features[0].identity + "_5UTR;Parent="+features[0].identity+"_mRNA"))
                                # We write the 3UTR feature in the GFF  : ID=identity_5UTR
                                gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","three_prime_UTR",cds_end+1,gene_end,".","+",".","ID=" + features[0].identity + "_3UTR;Parent="+features[0].identity+"_mRNA"))


                            # If the check option is False or the sequence just passed
                            # the check then we write it in the output:
                            try:
                                fwn.write(">{}\n{}\n".format(features[0].identity,str(features[0].seq_nucl.seq)))
                            except:
                                pass
                            try:
                                fwp.write(">{}\n{}\n".format(features[0].identity,str(features[0].seq_prot.seq)))
                            except:
                                pass
                            # We delete the previous feature
                            del features[0]

                        # Else we update the already existing feature
                        else:
                            features[-1].__update_feature__(feature)

        ## Final correction of the sequences:
        ##for feature in features:
        ##    feature.__correct_feature_sequences__(genetic_code=genetic_code)

        # Write the last feature left in the list
        features[0].__correct_feature_sequences__(genetic_code=genetic_code)
        try:
            fwn.write(">{}\n{}\n".format(features[0].identity,str(features[0].seq_nucl.seq)))
        except:
            pass
        try:
            fwp.write(">{}\n{}\n".format(features[0].identity,str(features[0].seq_prot.seq)))
        except:
            pass

        # And if elongate is not False, then write the last feature elongated
        if elongate != False:
            features[0].__elongate__(genome=genome,elongate=elongate)
            fasta_elongate.write(">{}\n{}\n".format(features[0].identity+"_mRNA",str(features[0].seq_nucl_elongated.seq)))
            gene_end  = len(features[0].seq_nucl_elongated.seq)
            cds_start = elongate+1
            cds_end   = elongate+len(str(features[0].seq_nucl.seq))
            # We write the GENE feature in the GFF : ID=identity
            gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","gene","1",gene_end,".","+",".","ID=" + features[0].identity))
            # We write the mRNA feature in the GFF : ID=identity_mRNA
            gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","mRNA","1",gene_end,".","+",".","ID="+features[0].identity+"_mRNA;Parent="+features[0].identity))
            # We write the CDS feature in the GFF  : ID=identity_CDS
            gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","CDS",cds_start,cds_end,".","+",".","ID=" + features[0].identity + "_CDS;Parent="+features[0].identity+"_mRNA"))
            # We write the 5UTR feature in the GFF  : ID=identity_5UTR
            gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","five_prime_UTR",features[0].UTR5_start,features[0].UTR5_end,".","+",".","ID=" + features[0].identity + "_5UTR;Parent="+features[0].identity+"_mRNA"))
            # We write the 3UTR feature in the GFF  : ID=identity_5UTR
            gff_elongate.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(features[0].identity+"_mRNA","elongated","three_prime_UTR",cds_end+1,gene_end,".","+",".","ID=" + features[0].identity + "_3UTR;Parent="+features[0].identity+"_mRNA"))


        # And close the opened fasta files
        try:
            fwn.close()
        except:
            pass
        try:
            fwp.close()
        except:
            pass
        try:
            fasta_elongate.close()
        except:
            pass
        try:
            gff_elongate.close()
        except:
            pass

        return(features)



def main():
    parameters     = get_args()
    genome_file    = parameters.fna

    #genome_file =  "/Users/christospapadopoulos/Documents/de_novo/Fungi_BLASTs/Scer/Scer.fna"
    print("Started\t:\t",time.ctime())
    my_fasta = SeqIO.to_dict(SeqIO.parse(open(genome_file),'fasta'))
    my_gff = GFF_iterator(gff_file      = parameters.gff,
                          genome        = my_fasta,
                          gff_types     = parameters.features_include,
                          outname       = parameters.o,
                          output_type   = parameters.type,
                          chr_exclude   = parameters.chr_exclude,
                          genetic_code  = parameters.table,
                          check         = parameters.check,
                          elongate      = parameters.elongate,
                          name_attribute= parameters.name_attribute)
    print("Ended \t:\t",time.ctime())
    return()


if __name__ == "__main__":
    main()
