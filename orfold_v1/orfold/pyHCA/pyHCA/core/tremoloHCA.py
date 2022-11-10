#!#/usr/bin/env python
""" this is an entirely revised version of tremoloHCA:

NEW features:
    - perform hca segmentation on fasta sequence before if not provided
    - option to select whole sequence or list of domain from hca segmentation
    - iterative search now use hhblits
    - display results by domain arrangement on a text file
    - textfile 2 html is a side utilitary script using tremolo text file result
"""

import os, sys, argparse
import configparser 
from pyHCA.core.annotateHCA import _annotation_aminoacids as segmentation
from pyHCA.core.ioHCA import read_multifasta, write_tremolo_results
from pyHCA.core.external import targets_hhblits, targets_jackhmmer_like, interpro_search
from pyHCA.core.classHCA import Seq

## domains
def read_domainpos(query, positions):
    """ read domain position in tuple format 1,10 20,30
    """
    domains = []
    if positions is None or positions == list():
        if isinstance(query.seq, str):
            # perform segmentation if it's only one sequence
            seg = segmentation(str(query.seq))
            for dom in seg["domain"]:
                domains.append((dom.start, dom.stop))
        else:
            #use orphhca on multiple sequence alignments
            # TODO
            print("Error, sequence ({}) not recognized for protein {}".format(query.name, query.seq), file=sys.stderr)
            sys.exit(1)
    elif positions[0] == "whole":
        # use the whole sequence
        domains = [(0, query.length)]
    else:
        # use user defined positions
        for val in positions:
            start, stop = val.split(",")
            start, stop = int(start)-1, int(stop)
            if stop<=start or start <0 or stop < 1:
                print("Error in start stop values {}".format(val), file=sys.stderr)
            domains.append((start, stop))
    return domains

def hhblits_search(query, domains, database, parameters, workdir):
    """ run hhblits
    """
    targets = dict()
    alltargetids = set()
    for i, (start, stop) in enumerate(domains):
        #print("domain {}".format(i))
        pathdom = os.path.join(workdir, "dom_{}".format(i))
        if not os.path.isdir(pathdom):
            os.makedirs(pathdom)
        # use sub part of sequence to search for targets
        pathquery = os.path.join(pathdom, "query_{}.fasta".format(i))
        with open(pathquery, "w") as outf:
            #for j, name in enumerate(query.name):
            subseq = str(query.seq)[start: stop]
            # IMPORTANT: the name of the sequence will be used as input for hhblits
            # a regular expression is set on "Q query_" to catch input name
            outf.write(">query_{} {} {}-{}\n{}\n".format(i, query.name, start+1, stop, subseq))
        # perform hhblits
        subtargets = targets_hhblits(pathquery, pathdom, database, parameters)
        alltargetids = alltargetids.union(set(subtargets.keys()))
        targets[i] = subtargets
    return targets, alltargetids

def jackhmmer_like_search(query, domains, database, parameters, workdir):
    """ run hhblits
    """
    targets = dict()
    alltargetids = set()
    for i, (start, stop) in enumerate(domains):
        #print("domain {}".format(i))
        pathdom = os.path.join(workdir, "dom_{}".format(i))
        if not os.path.isdir(pathdom):
            os.makedirs(pathdom)
        # use sub part of sequence to search for targets
        pathquery = os.path.join(pathdom, "query_{}.fasta".format(i))
        with open(pathquery, "w") as outf:
            #for j, name in enumerate(query.name):
            #subseq = str(query.seq[j])[start: stop]
            subseq = str(query.seq)[start: stop]
            # IMPORTANT: the name of the sequence will be used as input for hhblits
            # a regular expression is set on "Q query_" to catch input name
            outf.write(">query_{} {} {}-{}\n{}\n".format(i, query.name, start+1, stop, subseq))
        # perform hhblits
        subtargets = targets_jackhmmer_like(pathquery, pathdom, database, parameters)
        alltargetids = alltargetids.union(set(subtargets.keys()))
        targets[i] = subtargets
    return targets, alltargetids

## targets
def search_domains(query, domains, database, parameters, workdir):
    """ look for targets
    """
    #if method == "hhblits":
    targets, alltargetids = hhblits_search(query, domains, database, parameters, workdir)
    #elif method =="jackhmmer_like":
    #    targets, alltargetids = jackhmmer_like_search(query, domains, database, parameters, workdir)
    #else:
    #    raise ValueError("Unable to find method, check -m option [hhblits, or jackhmmer_like")
    #    sys.exit(0)
    return targets, list(alltargetids)

#### This is a useless comment on line 110 !

## group results by domain arrangements
def group_resda(targets, cddres):
    """ group results per domain arrangements
    """
    groups = dict()
    for querydom in targets:
        groups[querydom] = dict()
        for prot in targets[querydom]:
            if prot in cddres:
                cddres[prot].sort()
                da = ";".join(set([elmnt[2] for elmnt in cddres[prot]]))
                groups[querydom].setdefault(da, list()).append(prot)
            else:
                groups[querydom].setdefault("None", list()).append(prot)
    return groups

def get_cmd():
    """ get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputfasta", 
            help="input fasta file", required=True)
    parser.add_argument("-d", action="store", dest="domains", nargs="+", 
            help="list of domain positions (start and stop inclusive and "
            "separated by comma : -d 1,10 20,30 60,100. If not provided "
            "the search will be performed on each domain found after "
            "segmentation of the input sequence. "
            "To use the whole protein use -d whole.", default=list())
    parser.add_argument("-w", action="store", dest="workdir",
            help="working directory", required=True)
    #parser.add_argument("-m", action="store", dest="method", 
    #        choices=["hhblist", "jackhmmer_like"], default="hhblits",
    #        help="define sequence search method to use (default=%(default)s)")    
    parser.add_argument("--p2ipr", action="store", dest="p2ipr",
            help="path to the Interpro annotation of UniproKBt proteins, "
                 "gzip format supported.")
    parser.add_argument("--config", action="store", dest="configfile", 
            help="path to the configuration file for optional arguments")
    parser.add_argument("-o", action="store", dest="output", 
            help="output file")
    parser.add_argument("--target-db", action="store", dest="targetdb", 
            help="path to the target sequences database (fasta file for jackhmmer_like, hhm for hhblits", required=True)
    params = parser.parse_args()
        
    return params

def default_config(config):
    """ create a default configuration
    """
    config.optionxform=str
    config["path"] = dict()
    config["path"]["hhblits"] = "hhblits"
    config["path"]["phhmer"] = "phhmer"
    config["path"]["hmmbuild"] = "hmmbuild"
    config["path"]["hmmsearch"] = "hmmsearch"
    config["hhblits_options"] = dict()
    config["hhblits_options"]["E"] = "1.0"
    config["hhblits_options"]["cov"] = "10"
    config["hhblits_options"]["id"] = "80"
    config["jackhmmer_like_options"] = dict()
    config["jackhmmer_like_options"]["coverage"] = "10"
    config["jackhmmer_like_options"]["identity"] = "10"
    config["jackhmmer_like_options"]["kept_reg"] = ""
    config["jackhmmer_like_options"]["remove_reg"] = ""
    config["phmmer_options"] = dict()
    config["phmmer_options"]["E"] = "1.0"
    config["phmmer_options"]["domE"] = "1.0"
    config["hmmsearch_options"] = dict()
    config["hmmsearch_options"]["E"] = "1.0"
    config["hmmsearch_options"]["domE"] = "1.0"

def config_setup(path):
    """ process config file and user parameters
    """
    config = configparser.ConfigParser()
    if path is not None:
        config.read(path)
    else:
        default_config(config)
    return config

#### MAIN
def main():
    # main tremolo program
    params = get_cmd()
    configuration = config_setup(params.configfile)
    
    if not os.path.isdir(params.workdir):
        os.makedirs(params.workdir)
  
    # read input sequence
    inputquery = read_multifasta(params.inputfasta)
    names, seqs, descrs = list(), list(), list()
    for record in inputquery:
        names.append(inputquery[record].id)
        descrs.append(inputquery[record].description)
        seqs.append(str(inputquery[record].seq))
    
    if len(seqs) == 0:
        print("Error, no fasta sequences found in inputfile {}".format(
            params.inputfasta), sys.stder)
        sys.exit(1)
        
    is_multifasta = False
    if len(seqs) > 1:
        is_multifasta = True
    
    # output for multifasta 
    if is_multifasta :
        if os.path.isfile(params.output):
            print("Error, multifasta input file detected, output should be a directory", file=sys.stderr)
            sys.exit(1)
        else:
            if not os.path.isdir(params.output):
                os.makedirs(params.output)
    else:
        if os.path.isdir(params.output):
            print("Error, fasta file with a single entry detected, output should be a file", file=sys.stderr)
            sys.exit(1)
        else:
            output_file = params.output 
    
    # check domains positions
    if params.domains is None or params.domains == list():
        domains = [None] * len(seqs)
    else:
        if params.domains[0] == "whole":
            domains = ["whole"] * len(seqs)
        elif os.path.isfile(params.domains[0]):
            with open(params.domains[0]) as inf:
                prot2doms = dict()
                for line in inf:
                    tmp = line.split()
                    prot2doms.setdefault(tmp[0], list()).append("{},{}".format(tmp[1], tmp[2]))
            domains = []
            for prot in names:
                if prot in prot2doms:
                    domains.append(prot2doms[prot])
                else:
                    print("Warning, protein {} not found in domain list file {}".format(prot, params.domains[0]), end=". ")
                    print("Running HCA for this domain")
                    domains.append(None)
            assert len(domains) == len(seqs), "Error, the number of domain found is different than the number of sequence"
        else:
            if not is_multifasta:
                domains = params.domains[:]
            else:
                print("Error, a list of domain positions can only be specified "
                      "along a fasta file with a single protein (multi protein "
                      "sequences fasta file found)",  file=sys.stderr)
                sys.exit(1)
                
    for i in range(len(seqs)):
        query_workdir = os.path.join(params.workdir, "workdir_protein_{}".format(i+1))
        # only if input file has multiple fasta entries, then output_file is iteratively updated
        if len(seqs) > 1:
            output_file = os.path.join(params.output, "tremolo_output_{}.dat".format(i+1))

        if not os.path.isdir(query_workdir):
            os.makedirs(query_workdir)
            
        query = Seq(names[i], descrs[i], seqs[i], len(seqs[i]))
        # domains? whole sequence? segmentation?
        positions = read_domainpos(query, domains[i])
    
        # perform search method on each selected parts
        targets, alltargetids = search_domains(query, positions, params.targetdb, configuration, query_workdir)
        if alltargetids == []:
            print("Unable to find any targets with hhblits in database {}".format(params.targetdb), file=sys.stderr)
            #print("with parameters {}".format(params.hhblitsparams), file=sys.stderr)
            print("Please try less stringent parameters or a different database", file=sys.stderr)
            with open(output_file, "w") as outf: # output to output file
                outf.write("# Unable to find any targets with hhblits in database {}\n".format(params.targetdb))
                #outf.write("# with parameters {}\n".format(params.hhblitsparams))
                outf.write("# Please try less stringent parameters or a different database\n")
        
            sys.exit(0)
        
        # get domain from Interpro
        annotation = interpro_search(alltargetids, query_workdir, params.p2ipr)

        # group by domain arrangement
        groups = group_resda(targets, annotation)

        # write output
        write_tremolo_results(query, positions, targets, annotation, groups, output_file)

    sys.exit(0)

def foobar():
    print("pwet")

if __name__ == "__main__":
    main()


