#!/usr/bin/env python
""" The ioHCA module regroups the functions linked to file inputs and outputs
"""
from __future__ import print_function
import os, sys, gzip
from Bio import SeqIO
from collections import OrderedDict

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "t.bitard.feildel [you know what] uni-muenster.de"
__institute__ = "Institute for Evolution and Biodiversity, Muenster Germany"


def read_singlefasta(path, verbose=False):
    """ use Bio.SeqIO to read the fasta file and assert that it only contains 
    one protein
    """
    cnt = 0
    query, sequence = None, None
    for record, seq in read_multifasta_it(path, verbose=False):
        if cnt > 0:
            raise RuntimeError("Error, multiple sequences found in file {}".format(path))
        query = record
        sequence = str(seq.seq)
        cnt += 1
    return query, sequence

def read_multifasta(path, verbose=False):
    """ use Bio.SeqIO to read the fasta file and convert it into a dictionary
    
    Parameter
    ---------
    path : string
        path to the fasta file
        
    Return
    ------
    record_dict : dict
        a dictionary containing the Bio sequence object
    """
    if verbose:
        print("Read fasta inputfile ..")
    record_dict  = OrderedDict()
    if os.path.splitext(path)[1] in [".gz", ".gzip"]:
        with gzip.open(path, 'rt', encoding='utf-8') as handle: #Python3 fix
            for rec in SeqIO.parse(handle, "fasta"):
                record_dict[rec.id] = rec
            #record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        with open(path, "rU") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                record_dict[rec.id] = rec
            #record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return record_dict
    
def read_multifasta_it(path, verbose=False):
    """ use Bio.SeqIO to read the fasta file and convert it into a dictionary
    
    Parameter
    ---------
    path : string
        path to the fasta file
        
    Return
    ------
    record_dict : dict
        a dictionary containing the Bio sequence object
    """
    if verbose:
        print("Read fasta inputfile ..")
    if os.path.splitext(path)[1] in [".gz", ".gzip"]:
        with gzip.open(path, 'rt', encoding='utf-8') as handle: #Python3 fix
            for rec in SeqIO.parse(handle, "fasta"):
                yield rec.id, rec
    else:
        with open(path, "rU") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                yield rec.id, rec

def write_annotHCA(output, dannotate, sizes, verbose=False):
    """ write annotation output
    
    Parameters
    ----------
    outputf: string
        path to the output file
    dannotate: dict
        the annotation per proteins
    sizes: dict
        protein sizes
    verbose: bool
        print useful stuff
    """
    if verbose:
        print("Writting output of annotation to file {}".format(output))
            
    with open(output, "w") as outf:
        for prot in dannotate:
            outf.write(">{} {}\n".format(prot, sizes[prot]))
            for annotation in dannotate[prot]["domain"]:
                outf.write("{}\n".format(str(annotation)))
            for annotation in dannotate[prot]["cluster"]:
                outf.write("{}\n".format(str(annotation)))

def read_tremolo(path, fetch="domain"):
    """ read tremolo results
    """
    if fetch=="domain":
        return read_tremolo_domains(path)
    elif fetch=="all":
        return read_tremolo_all(path)
    else:
        print("Tremolo results fetch parameter {} is not supported".format(fetch), file=sys.stderr)
        sys.exit(1)
        
def read_tremolo_all(path):
    """ read table tremolo results
    """
    proteins = set()
    order = list()
    hits = dict()
    domains = dict()
    prot = None
    with open(path) as inf:
        for line in inf:
            if line[0] == "\n" and line[0] == "#":
                prot = None
                continue
            tmp = line.strip().split("\t")
            if line[0] == ">":
                prot, size = line[1:].split()
            elif prot != None and line.startswith("Qdom"):
                qdom = tmp[0].split()[1]
                domains.setdefault(prot, list())
                proteins.add(prot)
                if tmp[1] != "None":
                    dom = tmp[1]
                    start = int(tmp[2]) - 1
                    stop = int(tmp[3])
                    domain = (start, stop, dom)
                    domains.setdefault(prot, list()).append(domain)
            elif tmp[0] == "Hit":
                #"score", prot, dom, e_val, prob, score, ident, sim
                qdom = tmp[1]
                hitnum = int(tmp[2])
                order.append((prot, hitnum))
                hits.setdefault(prot, dict())
                hits[prot][hitnum] = {"E-value":float(tmp[3]), "Probab":float(tmp[4]), "Bitscore":float(tmp[5]), "Identities":float(tmp[6]), "Similarity":float(tmp[7])}
            elif tmp[0] == "HitQali":
                qdom = tmp[1]
                hitnum = int(tmp[2])
                hits[prot][hitnum]["Qstart"] = int(tmp[3])-1
                hits[prot][hitnum]["Qstop"] = int(tmp[4])
                hits[prot][hitnum]["Qali"] = tmp[5]
            elif tmp[0] == "HitTali":
                qdom = tmp[1]
                hitnum = int(tmp[2])
                hits[prot][hitnum]["Tstart"] = int(tmp[3])-1
                hits[prot][hitnum]["Tstop"] = int(tmp[4])
                hits[prot][hitnum]["Tali"] = tmp[5]
    
    return order, hits, domains, proteins
        
def read_tremolo_domains(path):
    """ read tremolo domain results
    """
    domains = dict()
    dsizes = dict()
    Tname = None
    protein_hit = False
    with open(path) as inf :
        for line in inf:
            #print(line)
            if line[0] == "\n" or line[0] == "#":
                Tname = None
                continue
            tmp = line.strip().split("\t")
            if line.startswith("Qdom") and len(tmp) == 4 and not protein_hit:
                domain = tmp[1]
                start, stop = int(tmp[2])-1, int(tmp[3])
                domains.setdefault(domain, dict())
                domains[domain]["QPos"] = (start, stop)
            elif line.startswith(">"):
                protein_hit = True
                Tname, Tsize = line[1:].strip().split()
                dsizes[Tname] = int(Tsize)
            elif line.startswith("Qdom") and Tname != None and protein_hit:
                domain = tmp[1]
                Tdomain = tmp[2]
                if Tdomain.startswith("IPR"):
                    Tdomain = Tdomain.split("/")[0]
                domains[domain].setdefault(Tname, dict())
                start, stop = int(tmp[3])-1, int(tmp[4])
                domains[domain][Tname].setdefault("Tpos", list()).append((start, stop, Tdomain))
            elif tmp[0] == "Hit" and Tname != None:
                domain = tmp[1]
                domains[domain].setdefault(Tname, dict())
                hitnb = tmp[2]
                evalue = float(tmp[3])
                domains[domain][Tname].setdefault("Hit", dict()).setdefault(hitnb, dict())
                domains[domain][Tname]["Hit"][hitnb]["evalue"] = evalue
            elif tmp[0] == "HitQali" and Tname != None:
                domain = tmp[1]
                hitnb = tmp[2]
                start, stop = int(tmp[3])-1, int(tmp[4])
                domains[domain][Tname]["Hit"][hitnb]["Qali"] = (start, stop)
            elif tmp[0] == "HitTali" and Tname != None:
                domain = tmp[1]
                hitnb = tmp[2]
                start, stop = int(tmp[3])-1, int(tmp[4])
                domains[domain][Tname]["Hit"][hitnb]["Tali"] = (start, stop)
            elif tmp[0] == "/":
                protein_hit = False
    return domains, dsizes

def read_annotation(inputfile, formatf):
    """ read domain annotation from pfam or HCA domain file

    Parameters
    ----------
    inputfile: string
        path to the annotation
    formatf: string
        file format either pfam or seghca

    Return
    ------
    annotaton: dict
        the domain annotation
    """
    annotation = dict()
    if formatf == "seghca":
        annotation = read_hcadomain(inputfile)
    elif formatf == "pfam":
        annotation = read_pfamdomain(inputfile)
    else:
        print("Error, no function defined to read domains in format {}".format(formatf), file=sys.stderr)
        sys.exit(1)
    return annotation

def read_hcadomain(inputfile):
    """ read hca domain

    Parameters    
    ---------- 
    inputfile: string
        path to the annotation

    Return
    ------
    annotaton: dict
        the domain annotation
    """
    annotation = dict()
    with open(inputfile) as inf:
        for line in inf:
            if line[0] == "#" or line[0] == " " or line[0] == "\n":
                continue
            elif line[0] == ">":
                prot, size, score= line[1:-1].split()
                annotation[prot] = []
            else:
                tmp = line.split()
                if tmp[0] == "domain":
                    start, stop = int(tmp[1])-1, int(tmp[2])
                    annotation[prot].append((start, stop, tmp[0], "!", None))
    return annotation

def read_pfamdomain(inputfile):
    """ read pfam domain

    Parameters 
    ---------- 
    inputfile: string
        path to the annotation

    Return   
    ------ 
    annotaton: dict 
        the domain annotation
    """  
    annotation = dict()
    with open(inputfile) as inf:
        for line in inf: 
            if line[0] == "#" or line[0] == "\n":
                continue
            tmp = line.split()
            prot = tmp[0]
            name = tmp[1]
            start, stop = int(tmp[2]), int(tmp[3])
            status = "!"
            annotation.setdefault(prot, []).append((start, stop, name, status, None))
    return annotation      

## reorganize results
def flatres(targets, proteins):
    """ flatten results for sorting from a list of proteins
    """
    # flatten dictionary
    #flattargets = list()
    flattargets = dict()
    prot_orders = dict()
    dsize = dict()
    for name in proteins:
        flattargets[name] = dict()
        for hitnum in targets[name]:
            res = targets[name][hitnum]
            if name in prot_orders:
                if res["E-value"] < prot_orders[name]:
                    prot_orders[name] = res["E-value"]
            else:
                prot_orders[name] = res["E-value"]
            flattargets[name][hitnum] = [res["E-value"], res["descr"], 
                res["Probab"], res["Score"], res["Identities"], res["Similarity"], res["Sum_probs"],
                res["Qstart"], res["Qstop"], res["Tstart"], res["Tstop"], 
                res["Qali"], res["Qcons"], res["Tali"], res["Tcons"]]
            if name not in dsize:
                dsize[name] = res["Tsize"]
            else:
                if dsize[name] != res["Tsize"]:
                    print("Warning different sizes ({} and {}) found for target protein {} in hhblits results".format(dsize[name], res["Tsize"], name))
            #flattargets.append([name, hitnum, res["E-value"], res["descr"], 
                #res["Probab"], res["Score"], res["Identities"], res["Similarity"], res["Sum_probs"],
                #res[d"Qstart"], res["Qstop"], res["Tstart"], res["Tstop"], 
                #res["Qali"], res["Qcons"], res["Tali"], res["Tcons"]])
    flatorders = sorted([(prot_orders[prot], prot) for prot in prot_orders])
    evalues, order = zip(*flatorders) 
    order_and_size = [(prot, dsize[prot]) for prot in order]
    #flattargets.sort()
    return order_and_size, flattargets

def orderda(arrangements, targets):
    """ order domain arrangements according to best evalue in the set of hit
    """
    keptda = list()
    for da in arrangements:
        #if da == None:
            #continue
        for prot in arrangements[da]:
            for hitnum in targets[prot]:
                keptda.append((targets[prot][hitnum]["E-value"], da))
    # sort 
    orderedda = list()
    keptda.sort()
    visited = dict()
    for evalue, da in keptda:
        if da not in visited:
            orderedda.append(da)
            visited[da] = 1
    return orderedda
            
def summary_tremolo(targets, groups, xbest=100):
    """ write a short summary table for the X best hits

    Parameters
    ----------
    targets: dict
        contains for each query domain, the protein and the hits from hhblits
    groups: dict
        proteins grouped by domain arrangement
    xbest: int
        number of hits to display

    Return
    ------
    summary: string
        the summary table in string, ready to be written
    """
    summary = ""
    for querydom in targets:
        prot2grp = dict()
        for da in groups[querydom]:
            for prot in groups[querydom][da]:
                prot2grp[prot] = da
        cnt = 0
        data = []
        for prot in targets[querydom]:
            for hitnum in targets[querydom][prot]:
                data.append((targets[querydom][prot][hitnum]["E-value"],
                             targets[querydom][prot][hitnum]["Probab"],
                             prot, hitnum, prot2grp[prot]))
        data.sort()
        data = data[:xbest]
        for evalue, probab, prot, hitnum, da in data:
            summary += "Qdom\t{}\t{}\t{}\t{}\t{}\t{}\n".format(querydom+1, prot, hitnum, da, evalue, probab)
    return summary

def write_tremolo_results(query, positions, targets, cddres, groups, output, xbest=100):
    """ write grouped results for domain res
    
    Parameters 
    ---------- 
    targets: dict
        contains for each query domain, the protein and the hits from hhblits
    cddres: dict
        contains the domain annotation
    groups: dict
        group the proteins per domain arrangement
    output: string
        path to the output file

    """
    with open(output, "w") as outf:
        # header
        outf.write("# Traveling into REmote hoMOLOgy with HCA (TREMOLO) \n")
        outf.write("# tremolo v1.0 was developped by Guillem Faure at Isabelle Callebaut's team\n")
        outf.write("# tremolo v2.0 is developped and maintained by Tristan Bitard-Feildel\n")
        outf.write("#\n")
        outf.write("# Please cite\n")
        outf.write("# Identification of hidden relationships from the coupling of hydrophobic cluster analysis and domain architecture information.\n")
        outf.write("# G. Faure & I. Callebaut.\n")
        outf.write("Bioinformatics. 2013 Jul 15;29(14):1726-33. doi: 10.1093/bioinformatics/btt271.")
        outf.write("\n\n")

        # write query information
        for i, name in enumerate(query.name):
            outf.write("Qname\t{}\n".format(name))
            outf.write("Qdesc\t{}\n".format(query.descr[i]))
            outf.write("Qseq\t{}\n\n".format(query.seq[i]))
        # write domain positions
        querydoms = []
        for ite, (start, stop) in enumerate(positions):
            outf.write("Qdom\t{}\t{}\t{}\n".format(ite+1, start+1, stop))
            querydoms.append(ite)
        outf.write("\n")

        # write summary
        outf.write(summary_tremolo(targets, groups, xbest))
        outf.write("\n")
        for querydom in querydoms:
            outf.write("# Query domain {}\n\n".format(querydom+1))
            # order da depending on best evalue
            orderedda = orderda(groups[querydom], targets[querydom])
            # write all domain arrangement at the beginning an the number of proteins
            outf.write("# Domain Domain_arrangement number_of_protein\n")
            for da in orderedda:
                outf.write("INFO\t{}\t{}\t{}\n".format(querydom+1, da, len(groups[querydom][da])))
            outf.write("\n")
            for da in orderedda:
                # sort prot by evalues
                proteins = groups[querydom][da]
                outf.write("##\n")
                order, flat = flatres(targets[querydom], proteins)
                for (prot, size) in order:
                    outf.write(">{} {}\n".format(prot, size))
                    if prot in cddres:
                        for start, stop, dom, d_e_val, bitscore, types in cddres[prot]:
                            outf.write("Qdom\t{}\t{}\t{}\t{}\n".format(querydom+1, dom, start+1, stop))
                    else:
                        outf.write("Qdom\t{}\t{}\n".format(querydom+1, "None"))
                    for hit in flat[prot]:
                        e_val, descr, prob, score, ident, sim, sprob, qstart, qstop, tstart, tstop, qali, qcons, tali, tcons = flat[prot][hit]
                        outf.write("Hit\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(querydom+1, hit, e_val, prob, score, ident, sim))
                        outf.write("HitQali\t{}\t{}\t{}\t{}\t{}\n".format(querydom+1, hit, qstart+1, qstop, qali))
                        outf.write("HitQcon\t{}\t{}\t{}\t{}\t{}\n".format(querydom+1, hit, qstart+1, qstop, qcons))
                        outf.write("HitTcon\t{}\t{}\t{}\t{}\t{}\n".format(querydom+1, hit, tstart+1, tstop, tcons))
                        outf.write("HitTali\t{}\t{}\t{}\t{}\t{}\n".format(querydom+1, hit, tstart+1, tstop, tali))
                        outf.write("//\n")
