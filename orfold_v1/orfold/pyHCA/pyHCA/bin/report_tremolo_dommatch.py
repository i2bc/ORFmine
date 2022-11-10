#!/usr/bin/env python
""" report tremolo results directly matching an annotated protein domain
"""

import os, sys, argparse
from pyHCA.core.ioHCA import read_tremolo

def get_cmd():
    """ get command line argument
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="tremolores")
    parser.add_argument("--qcov", action="store", dest="min_qcov", type=float, default=0.8)
    parser.add_argument("--tcov", action="store", dest="min_tcov", type=float, default=0.5)
    parser.add_argument("--evalue", action="store", dest="evalue", type=float, default=0.001)
    params = parser.parse_args()
    return params

def find_dommatch(tremolo_res, qcov, tcov, cutoff_evalue):
    """ find positions matchin a domain
    """
    print("Domain Id\tTarget protein\tDomain\tDom beg\tDom end\t"
          "Thit beg\tThit end\tQhit beg\tQhit end\tevalue\ttarget cov.\tquery cov.")
    for domain in tremolo_res:
        qdom_start, qdom_stop = tremolo_res[domain]["QPos"]
        del tremolo_res[domain]["QPos"]
        for prot in tremolo_res[domain]:
            target_domains = tremolo_res[domain][prot].get("Tpos", [])
            for hit in tremolo_res[domain][prot]["Hit"]:
                evalue = tremolo_res[domain][prot]["Hit"][hit]["evalue"]
                if evalue < cutoff_evalue:
                    thit_start, thit_stop = tremolo_res[domain][prot]["Hit"][hit]["Tali"]
                    qhit_start, qhit_stop = tremolo_res[domain][prot]["Hit"][hit]["Qali"]
                    for tdom_start, tdom_stop, tdom in target_domains:
                        start = max(thit_start, tdom_start)
                        stop = min(thit_stop, tdom_stop)
                        diff = stop-start
                        if diff > 0:
                            c1 = diff / (tdom_stop-tdom_start)
                            c2 = diff / (qdom_stop-qdom_start)
                            if c1 >= tcov and c2 >= qcov:
                                print("Domain {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(domain, prot, tdom, tdom_start-1, tdom_stop, thit_start, thit_stop, qhit_start, qhit_stop, evalue, c1, c2))
    
def main():
    params = get_cmd()
    
    tremolo_res, dsizes = read_tremolo(params.tremolores, fetch="domain")
    find_match = find_dommatch(tremolo_res, params.min_qcov, params.min_tcov, params.evalue)
    
    sys.exit(0)
    
if __name__ == "__main__":
    main()