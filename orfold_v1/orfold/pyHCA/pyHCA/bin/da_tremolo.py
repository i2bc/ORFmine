#!/usr/bin/env python
""" report tremolo domain arrangement with unknown domain position
"""

import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_cmd():
    """ get command line argument
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="tremolores", help="tremolo result file")
    parser.add_argument("-c", action="store_true", dest="collapse", help="collapse domain arrangement")
    parser.add_argument("-o", action="store", dest="output", help="output plot with covered domains")
    parser.add_argument("--tcov", action="store", dest="min_tcov", type=float, default=0.5, help="minimal coverage of annotated domains to match")
    parser.add_argument("--evalue", action="store", dest="max_evalue", type=float, default=0.001, help="maximal evalue allowed")
    params = parser.parse_args()
    return params

def read_tremolo(path):
    """ read tremolo domain results
    """
    domains = dict()
    Tname = None
    with open(path) as inf :
        for line in inf:
            #print(line)
            if line[0] == "\n" or line[0] == "#":
                Tname = None
                continue
            tmp = line.strip().split("\t")
            if line.startswith("Qdom") and len(tmp) == 4:
                domain = tmp[1]
                start, stop = int(tmp[2])-1, int(tmp[3])
                domains.setdefault(domain, dict())
                domains[domain]["QPos"] = (start, stop)
            elif line.startswith(">"):
                Tname = line[1:].strip()
            elif line.startswith("Qdom") and Tname != None:
                domain = tmp[1]
                Tdomain = tmp[2]
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
    return domains

def collapse(da):
    da_collapsed = [da[0]]
    for i in range(1, len(da)):
        if da[i] != da_collapsed[-1]:
            da_collapsed.append(da[i])
    return da_collapsed
    

def filter_byda(tremolo_res, tcov, tevalue):
    """ filter results by domain arrangements
    """
    domain_arrangements = dict()
    collapsed_domain_arrangements = dict()
    for domain in tremolo_res:
        qdom_start, qdom_stop = tremolo_res[domain]["QPos"]
        del tremolo_res[domain]["QPos"]
        for prot in tremolo_res[domain]:
            target_domains = tremolo_res[domain][prot].get("Tpos", [])
            for hit in tremolo_res[domain][prot]["Hit"]:
                evalue = tremolo_res[domain][prot]["Hit"][hit]["evalue"]
                if evalue < tevalue:
                    thit_start, thit_stop = tremolo_res[domain][prot]["Hit"][hit]["Tali"]
                    qhit_start, qhit_stop = tremolo_res[domain][prot]["Hit"][hit]["Qali"]
                    target_domains.append((thit_start, thit_stop, "Hit_{}-{}_{}-{}".format(thit_start, thit_stop, qhit_start, qhit_stop)))
            target_domains.sort()
            da_gen = list()
            da = list()
            for tdom_start, tdom_stop, tdom in target_domains:
                if "/" in tdom and tdom.startswith("IPR"):
                    tdom = tdom.split("/")[0]
                da.append(tdom)
                if not tdom.startswith("Hit"):
                    da_gen.append(tdom)
            if da_gen != list():
                da_collapsed = collapse(da)
                da_gen_collapsed = [dom for dom in da_collapsed if not dom.startswith("Hit")]
                da_gen = ";".join(da_gen)
                domain_arrangements.setdefault(da_gen, list()).append(da)
                da_gen_collapsed = ";".join(da_gen_collapsed)
                collapsed_domain_arrangements.setdefault(da_gen_collapsed, list()).append(da_collapsed)
    return domain_arrangements, collapsed_domain_arrangements

def main():
    params = get_cmd()
    
    tremolo = read_tremolo(params.tremolores)
    
    domain_arrangements, domain_arrangements_co = filter_byda(tremolo, params.min_tcov, params.max_evalue)
    
    if params.collapse:
        for da_gen in domain_arrangements_co:
            print("{} {} {}".format("#"*6, da_gen, "#"*6))
            for da in domain_arrangements_co[da_gen]:
                print(" ".join(da))
            print()
    else:
        for da_gen in domain_arrangements:
            print("{} {} {}".format("#"*6, da_gen, "#"*6))
            for da in domain_arrangements[da_gen]:
                print(" ".join(da))
            print()
    
    sys.exit(0)


if __name__ == "__main__":
    main()
    