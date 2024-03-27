#!/usr/bin/env python
""" add domains to tree
"""

from pyHCA.core.ioHCA import read_tremolo
import os, sys, argparse, ete3, json, gzip
import numpy as np
import random
import requests
from Bio import SeqIO
from ete3 import NCBITaxa
#try:
from ete3 import SeqMotifFace, TreeStyle, NodeStyle, add_face_to_node
#except ImportError:
#    print('An error occured while importing ete3, please check that ete3 is correctly installed allong PyQt4 (run python -c "from ete3 import SeqMotifFace"', file=sys.stderr)
#    sys.exit(1)
    

def retrieve_sp(proteins):
    """ retrieve uniprot data from API
    """
    prot2taxid = dict()
    taxid2sp = dict()
    to_remove = list()
    #cmd = "http://www.uniprot.org/uniprot/?query=id:{}&sort=score&columns=id,entry%20name,length,reviewed,organism-id,organism&format=json"
    cmd = "http://www.uniprot.org/uniprot/?query=id:{}&sort=score&columns=id,entry%20name,length,reviewed,organism-id,organism&format=tab"
    #Entry	Entry name	Length	Status	Organism ID	Organism
    #H0H275	H0H275_SACCK	245	unreviewed	1095631	Saccharomyces cerevisiae x Saccharomyces kudriavzevii (strain VIN7) (Yeast)
    for fullprot in proteins:
        prot = fullprot.split("|")[1]
        query = cmd.format(prot).strip()
        r = requests.get(query)
        if r.status_code == 200:
            #prot_data = r.json()[0]
            lines = r.content.decode().split("\n")
            if len(lines) >= 2:
                tmp = lines[1].strip().split("\t")
                if len(tmp) >= 4:
                    #print(tmp)
                    taxid = tmp[4]
                    sp = tmp[5]
                    #taxid = prot_data["organism-id"]
                    #sp = prot_data["organism"]
                    if taxid == "" or sp == "":
                        print("Unable to retrieve species information for protein {}".format(prot), file=sys.stderr)
                        to_remove.append(fullprot)
                    else:
                        if "(" in sp:
                            start_parenthesis = sp.index("(")
                            sp = sp[:start_parenthesis].strip()
                        prot2taxid[fullprot] = taxid
                        taxid2sp[taxid] = sp
                else:
                    print("Unable to retrieve uniprot data for protein {}".format(prot), file=sys.stderr)
                    to_remove.append(fullprot)
            else:
                print("Unable to retrieve uniprot data for protein {}".format(prot), file=sys.stderr)
                to_remove.append(fullprot)
        else:
            print("Unable to retrieve uniprot data for protein {}".format(prot), file=sys.stderr)
    if len(taxid2sp) == 0:
        print("Error, unable to retrieve the species of any of the protein provided, no tree can be constructed", file=sys.stderr)
        sys.exit(1)
    return prot2taxid, taxid2sp, to_remove

def read_taxid_uniprot(path):
    """ read taxid and uniprot informations
    the file is made of three columns separated by tabulations:
    <protein name>\t<ncbi taxonomic id>\t<species name>
    """
    prot2taxid = dict()
    taxid2sp = dict()
    with open(path) as inf:
        for line in inf:
            prot, taxid, sp = line.strip().split("\t")
            prot2taxid[prot] = taxid
            taxid2sp[taxid] = sp
    return prot2taxid, taxid2sp
            

def get_overlapping_domains(domains):
    """ get overlapping domains
    """
    to_remove = set()
    domains = sorted(list(domains))
    for i in range(len(domains)):
        starti, stopi, domi = domains[i]
        sizei = stopi - starti
        ki = (starti, stopi, domi)
        found = False
        if ki not in to_remove:
            for j in range(i+1, len(domains)):
                startj, stopj, domj = domains[j]
                sizej = stopj - startj
                kj = (startj, stopj, domj)
                if kj not in to_remove:
                    found = True
                    start = max(starti, startj)
                    stop = min(stopi, stopj)
                    diff = stop-start
                    if diff > 0: #overlapping 
                        if sizei > sizej:
                            to_remove.add(kj)
                        else:
                            to_remove.add(ki)
    to_keep = set()
    for i in range(len(domains)):
        starti, stopi, domi = domains[i]
        sizei = stopi - starti
        ki = (starti, stopi, domi)
        if ki not in to_remove:
            to_keep.add(ki)
    return to_remove, to_keep
    
def filter_domain_arch(domains):
    """ filter domain architecture based on domain overlaps
    """
    new_domains = list()
    for prot in domains:
        to_remove, new_domains = get_overlapping_domains(domains)
        while len(to_remove) > 0:
            to_remove, new_domains = get_overlapping_domains(new_domains)
    return new_domains
        

def random_color():
    """ generate a random RGB color
    """
    r = lambda: random.randint(0,255)
    #'#%02x%02x%02x' % (0, 128, 64)
    color = '#%02X%02X%02X' % (r(),r(),r())
    #color = np.random.rand(3,)
    return color

all_colors = set()
domain_color = dict()
def get_domain_color(dom):
    """ assign an uniq color to each domains
    """
    if dom not in domain_color:
        color = random_color()
        while color in all_colors:
            color = random_color()
        all_colors.add(color)
        domain_color[dom] = color
    return domain_color[dom]

def combine_features(data, dsizes, tree, taxid2sp, prot2taxid, taxa_to_merge):
    """ create a ete3 tree with domain features from jackhmmer
    """ 
    # motif example from http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#phylogenetic-trees-and-sequence-domains
    #simple_motifs = [
        ## seq.start, seq.end, shape, width, height, fgcolor, bgcolor
        #[10, 60, "[]", None, 10, "black", "rgradient:blue", "arial|8|white|long text clipped long text clipped"],
        #[120, 150, "o", None, 10, "blue", "pink", None],
        #[200, 300, "()", None, 10, "blue", "red", "arial|8|white|hello"],
    #]
    
    # add domain match
    # and ount number of sequences by taxid 
    # and get size
    #dsize = dict()
    motifs = dict()
    dcnt_sp = dict()
    for tremolo_dom in data:
        tremolo_dom_start, tremolo_dom_stop = data[tremolo_dom]["QPos"]
        del data[tremolo_dom]["QPos"]
        for prot in data[tremolo_dom]:
            taxid = prot2taxid[prot]
            sp = taxid2sp[taxid]
            sp = sp.replace("(", "").replace(")", "").replace(",", "").replace(";", "")
            taxid2sp[taxid] = sp
            
            full_domains = data[tremolo_dom][prot].get("Tpos", list())[:]
            
            target_domains = filter_domain_arch(full_domains)
            
            dom_arch = list()
            ordered_domains = sorted(target_domains)
            for start, stop, dom in ordered_domains:
                dom_arch.append(dom)
            dom_arch = ";".join(dom_arch)

            if sp not in dcnt_sp:
                dcnt_sp[sp] = dict()
            if dom_arch not in dcnt_sp[sp]:
                dcnt_sp[sp][dom_arch] = list()
            dcnt_sp[sp][dom_arch].append(prot)
            
            # add domains
            for start, stop, dom in ordered_domains:
                color = get_domain_color(dom)
                motifs.setdefault(sp, dict())\
                        .setdefault(prot, list())\
                        .append([start+1, stop, "[]", None, 10, "black", color, "arial|1|black|{}".format(dom)])
                #motifs[sp][prot].sort()
            for hitnb in data[tremolo_dom][prot]["Hit"]:
                start, stop = data[tremolo_dom][prot]["Hit"][hitnb]["Tali"]
                motifs.setdefault(sp, dict())\
                        .setdefault(prot, list())\
                        .append([start, stop, "o", None, 10, "black", "red", "arial|1|black|HCAdom {}-{}".format(tremolo_dom, hitnb)])
            
    #print(motifs)
    
    # merge taxonomic groups
    for taxid_taxa in taxa_to_merge:
        taxid, taxa = taxid_taxa.split(",")
        lnode = tree.search_nodes(name=taxid)
        if len(lnode) > 0:
            taxnode = lnode[0]
            children = [child.name for child in taxnode.children]
            leaves = list()
            for child in taxnode.traverse():
                if child.is_leaf():
                    leaves.append(child)
                    #print(taxid_taxa, child.name)
            taxid2sp[taxid] = taxa
            dcnt_sp[taxa] = dict()
            motifs[taxa] = dict()
            for leafnode in leaves:
                nodeid = leafnode.name
                nodesp = taxid2sp[nodeid]
                #print(taxid_taxa, nodeid, nodesp)
                # merge protein list
                for dom_arch in dcnt_sp[nodesp]:
                    for prot in dcnt_sp[nodesp][dom_arch]:
                        dcnt_sp[taxa].setdefault(dom_arch, list()).append(prot)
                        
                # merge motif
                for prot in motifs[nodesp]:
                    for m in motifs[nodesp][prot]:
                        motifs[taxa].setdefault(prot, list()).append(m)
                # delete obsoletes species
                del dcnt_sp[nodesp]
                del motifs[nodesp]
                del taxid2sp[nodeid]
            #for dom_arch in dcnt_sp[taxa]:
                #for prot in dcnt_sp[taxa][dom_arch]:
                    #print(prot, dom_arch)
            for child in children:
                node = tree.search_nodes(name=child)[0]
                node.detach()
            #print(taxnode.get_ascii(show_internal=True))
        else:
            print("Unable to find taxid {} in tree".format(taxid))
            print(lnode)

    for taxid in taxid2sp:
        #print(taxid)
        node = tree.search_nodes(name=taxid)
        if node != []:
            node[0].name = taxid2sp[taxid]
        else:
            print("Unable to find node for taxid {}".format(taxid))
        
    # expand taxonomic tree by the number of sequences in each taxa
    for node in tree:
        if node.is_leaf():
            if node.name in dcnt_sp:
                node_sp = node.name
                proteins = list()
                features = dict()
                for dom_arch in dcnt_sp[node_sp]:
                    sizes_and_proteins = list()
                    for prot in dcnt_sp[node_sp][dom_arch]:
                        new_name = node_sp+" | "+prot.split("|")[1]
                        if len(dcnt_sp[node_sp][dom_arch]) > 1:
                            new_name += " [+{}]".format(len(dcnt_sp[node_sp][dom_arch])-1)
                        sizes_and_proteins.append(
                            (dsizes[prot], 
                             new_name,
                             dom_arch, prot)
                            )
                    sizes_and_proteins.sort(reverse=True)
                    sizes, names, dom_archs, prots = zip(*sizes_and_proteins)
                    proteins.append(names[0].replace(":", " "))
                    features[names[0].replace(":", " ")] = (prots[0], dom_archs[0])
                subtree = ete3.PhyloTree("({});".format(", ".join(proteins)))
                node.add_child(subtree)
                for new_node in subtree:
                    prot, dom_arch = features[new_node.name]
                    seq = "G" * dsizes[prot]
                    m = motifs[node_sp][prot]
                    seqFace = SeqMotifFace(seq, seq_format="line", motifs=m)
                    new_node.add_face(seqFace, 0, "aligned")
    #for dom in domain_color:
        #print(dom, domain_color[dom])
    return tree

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="tremolores", 
                        help="tremolo results with domain matchs")
    parser.add_argument("-t", action="store", dest="treefile", 
                        help="phylogenetic tree with node as ncbi taxonomic ids (optional)")
    parser.add_argument("-s", action="store", dest="prot2species", 
                        help="file with prot to species informations (optional)")
    parser.add_argument("-n", action="store", dest="ncbitaxid", nargs="+", default=list(),
                        help="list of node for which leaves will be merged (internal node need to be in tree)")
    parser.add_argument("-o", action="store", dest="output", 
                        help="phylogenetic tree with tremolo hits")
    params = parser.parse_args()
    return params

def main():
    params = get_cmd()
    
    # read jackhmmer data
    tremolores, target_sizes = read_tremolo(params.tremolores, fetch="domain")
    
    to_remove = list()
    if params.prot2species:
        prot2taxid, taxid2sp = read_taxid_uniprot(params.prot2species)
    else:
        proteins = list(target_sizes.keys())
        prot2taxid, taxid2sp, to_remove = retrieve_sp(proteins)
    
    if to_remove:
        for prot in to_remove:
            for qdom in tremolores:
                if prot in tremolores[qdom]:
                    del tremolores[qdom][prot]
            if prot in target_sizes:
                del target_sizes[prot]
    
    # read tree
    if params.treefile:
        tree = ete3.Tree(params.treefile, format=1)
        for node in tree.traverse():
            #if node.name == "INT8457":
                #print(node.name)
            if node.name.startswith("INT"):
                node.name = node.name[3:]
    else:
        # download tree from ncbi taxonomic db
        ncbi = NCBITaxa()
        tree = ncbi.get_topology(list(taxid2sp.keys()), intermediate_nodes=True)
        
        
    # check that taxid in tree
    for taxid in taxid2sp:
        node = tree.search_nodes(name=taxid)
        if len(node) < 1:
            print("Unable to find taxonomic id {} in tree".format(taxid), file=sys.stderr)
            sys.exit(1)
    
    # create tree with domain
    new_tree = combine_features(tremolores, target_sizes, tree, taxid2sp, prot2taxid, params.ncbitaxid)
    
    # clean tree
    to_remove = list()
    for node in tree.traverse():
        if not node.is_root() and not node.is_leaf():
            if len(node.children) == 1:
                to_remove.append(node)
    while to_remove != []:
        for node in to_remove:
            node.delete()
        to_remove = list()
        for node in tree.traverse():
            if not node.is_root() and not node.is_leaf():
                if len(node.children) == 1:
                    to_remove.append(node)

    
    most_distant_leaf, tree_length = new_tree.get_farthest_leaf()
    current_dist = 0
    for postorder, node in new_tree.iter_prepostorder():
        if postorder:
            current_dist -= node.dist
        else:
            if node.is_leaf():
                node.dist += tree_length - (current_dist + node.dist)
            elif node.up: # node is internal
                current_dist += node.dist
    
    ts = TreeStyle()   
    ts.show_leaf_name = True
    #ts.force_topology = True
    ts.show_scale = False
    for node in new_tree.traverse():
        nodestyle = NodeStyle()
        nodestyle["fgcolor"] = "black"
        nodestyle["size"] = 3
        nodestyle["vt_line_color"] = "black"
        nodestyle["hz_line_color"] = "black"
        nodestyle["vt_line_width"] = 2
        nodestyle["hz_line_width"] = 2
        nodestyle["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
        nodestyle["hz_line_type"] = 0
        
    new_tree.render(params.output, dpi=600, w=8000)
    
    sys.exit(0)
    
if __name__ == "__main__":
    main()
