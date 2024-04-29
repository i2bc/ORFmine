#!/usr/bin/env python
""" get sequence from domain annotation
"""

import os, sys, argparse, string, re
from pyHCA.core.ioHCA import read_multifasta_it, read_annotation
from Bio import Seq

__author__ = "Tristan Bitard-Feildel, Guillem Faure, Isabelle Callebaut"
__licence__= "MIT"
__version__ = 0.1
__email__ = "t.bitard.feildel [you know what] uni-muenster.de"
__institute__ = "Institute for Evolution and Biodiversity, Muenster Germany"


def domain_sequence(inputf, domainf, outputf, verbose=False):
    """ get sequences of domain annotation
    """
    # read domain annotation
    annotation = read_annotation(domainf, "seghca")
    with open(outputf, "w") as outf:
        for prot, sequence in read_multifasta_it(inputf):
            seq = str(sequence.seq)
            for i, domain in enumerate(annotation.get(prot, [])):
                start, stop = domain[:2]
                outf.write(">{}-{} {}-{}\n{}\n".format(prot, i, start+1, stop, seq[start: stop]))

def _process_params():
    """ Process parameters when the script annotateHCA is directly called
    """
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "domainseq"))
    parser.add_argument("-i", action="store", dest="inputf", required=True,
        help="an amino-acid sequence files in fasta format")
    parser.add_argument("-d", action="store", dest="domainf", required=True,
        help="the hca domain annotation")
    parser.add_argument("-o", action="store", dest="outputf", required=True,
        help="the output file with annotation")
    parser.add_argument("-v", action="store_true", dest="verbose", default=False,
        help="keep temporary results")
    params = parser.parse_args()
    return params

    
def main_domainseq():
    """ the main function is called after direct invocation of the software
    """
    params = _process_params()
    
    # annotation
    domain_sequence(params.inputf, params.domainf, params.outputf, verbose=params.verbose)
    sys.exit(0)

    
if __name__ == "__main__":
    main_domainseq()
    
