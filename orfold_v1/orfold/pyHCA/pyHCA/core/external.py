#!/usr/bin/env python

""" The searchHCA module is an ensemble of functions to search a specified 
hydrophobic cluster or a list of  hydrophobic cluster (both borned by start and 
stop positions) on sequence databases.
"""

import os, sys, time, re
import subprocess, shlex
import sqlite3, gzip
from Bio import AlignIO

debug = True

### INTERPRO

def interpro_search(listofids, workdir, path_p2ipr):
    """ get Interpro annotation from a list of ids
    """

    if path_p2ipr:
        return local_interpro_search(listofids, path_p2ipr)
    else:
        #return web_interpro_search(listofids, workdir)
        print("Interpro web search not yet implemented, please use the local "
              "interpro file available at "
              "ftp://ftp.ebi.ac.uk/pub/databases/interpro/protein2ipr.dat.gz")
        sys.exit(1)

#def read_annotation_sqlite3(uniprotids, path_p2ipr):
    #""" read interpro annotation stored in a sqlite3 database
    #"""

    #annotation = dict()
    #t1 = time.time()
    #conn = apsw.Connection(path_p2ipr)
                
    #cur = conn.cursor()
                        
    #for prot in uniprotids:
        #print(prot)
        #name = uniprotids[prot]
        #cur.execute("select start, stop, interprodom, domain from interpro where protein=?", (prot,))
        #for row in cur:
            #annotation.setdefault(name, list()).append((row[0], row[1], row[2]))

    #cur.close()
    #conn.close()

    #t2 = time.time()
    #print("Done in {}".format(t2-t1))
    #return annotation

def read_annotation(uniprotids, path_p2ipr):
    """ read interpro annotation
    """
    annotation = dict()
    t1 = time.time()
    if path_p2ipr[-2:] == "gz":
        with gzip.open(path_p2ipr, "rt") as inf:
            for line in inf:
                tmp = line.strip().split("\t")
                if tmp[0] in uniprotids:
                    start, stop = int(tmp[4])-1, int(tmp[5])
                    prot = uniprotids[tmp[0]]
                    annotation.setdefault(prot, list()).append((start, stop, tmp[1]+"/"+tmp[3], -1, -1, -1))
                elif len(annotation) == len(uniprotids):
                    return annotation
    else:
        with open(path_p2ipr, "r") as inf:
            for line in inf:
                tmp = line.strip().split("\t")
                if tmp[0] in uniprotids:
                    start, stop = int(tmp[4])-1, int(tmp[5])
                    prot = uniprotids[tmp[0]]
                    annotation.setdefault(prot, list()).append((start, stop, tmp[1]+"/"+tmp[3], -1, -1, -1))
                elif len(annotation) == len(uniprotids):
                    return annotation
    t2 = time.time()
    print("Done in {}".format(t2-t1))
    return annotation

def local_interpro_search(listofids, path_p2ipr):
    """ read local interpro file to look for domain annotations
    """
    annotation = dict()
    uniprotids = dict([(prot.split("|")[1], prot) for prot in listofids])
    cnt = 0
    annotation = read_annotation(uniprotids, path_p2ipr)
    for prot in uniprotids:
        name = uniprotids[prot]
        if name in annotation:
            annotation[name].sort()
        else:
            annotation[name] = list()
    return annotation

### HHBLITS

def filter_hhblits(res, evalue): 
    """ filter targets by only keeping target with at least one hit under <evalue>
    """
    targets = dict()
    for k in res:
        hitnums = list()
        for hitnum in res[k]:
            if res[k][hitnum]["E-value"] < evalue:
                hitnums.append(hitnum)
        if hitnums != list():
            targets[k] = dict()
            for hitnum in hitnums:
                targets[k][hitnum] = dict()
                for param in res[k][hitnum]:
                    targets[k][hitnum][param] = res[k][hitnum][param]
    return targets    

def targets_hhblits(pathquery, workdir, pathdb, parameters):
    """ run and read results of hhblits
    """
    pathout = run_hhblits(pathquery, workdir, pathdb, parameters)
    targets = read_hhblits(pathout)
    try:
        cutoff_evalue = parameters["hhblits_options"].getfloat("E")
    except:
        raise ValueError("Unable to read value for 'E' keys of 'hhblits_options' section in the config file")
        
    new_targets = filter_hhblits(targets, cutoff_evalue)
    return new_targets

def read_hhblits(pathin):
    """ read hhblits results
    """
    alltargets = set()
    targets = dict()
    hitnumber = dict()
    with open(pathin) as inf:
        for line in inf:
            if line[0] == ">":
                tmp = line[1:].split()
                name = tmp[0]
                descr = " ".join(tmp[1:])
                hitnum = 0
                targets.setdefault(name, dict())
                if name in hitnumber:
                    hitnum = hitnumber[name] + 1
                hitnumber[name] = hitnum
                alltargets.add(name)
                targets[name][hitnum] = {"descr": descr, "Tstart":1e10, "Tstop":-1, "Tcons":"", "Tali":"", "Tsize": 0,
                                                "Qstart":1e10, "Qstop":-1, "Qcons":"", "Qali":"", "Qsize": 0,
                                "Probab":-1, "E-value":-1, "Score":-1, "Identities":-1, "Similarity":-1, "Sum_probs":-1,
                                }
            elif line.startswith("Probab="):
                #Probab=100.00  E-value=6.2e-40  Score=297.36  Aligned_cols=130  Identities=100%  Similarity=1.267  Sum_probs=129.4
                tmp = line.split()
                for keyval in tmp:
                    key, val = keyval.split("=")
                    if key == "Identities":
                        val = val[:-1]
                    targets[name][hitnum][key] = float(val)
            elif line.startswith("Q ") and line.split()[1] != "Consensus":
                m = re.match("\s*(\d+)\s+([\-\w]*)\s+(\d+)\s\((\d+)\)",line[16:])
                if m:
                    start, ali, stop, size = int(m.group(1))-1, m.group(2), int(m.group(3)), int(m.group(4))
                    targets[name][hitnum]["Qstart"] = min(targets[name][hitnum]["Qstart"], start)
                    targets[name][hitnum]["Qstop"] = max(targets[name][hitnum]["Qstop"], stop)
                    targets[name][hitnum]["Qali"] += ali
                    targets[name][hitnum]["Qsize"] = size
                else:
                    print("FAILED", line, pathin)
            elif line.startswith("Q Consensus"):
                m = re.match("\s*(\d+)\s+([\-\~\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = m.group(1), m.group(2), m.group(3)
                    targets[name][hitnum]["Qcons"] += ali
                else:
                    print("FAILED", line, pathin)
            elif line.startswith("T Consensus"):
                m = re.match("\s*(\d+)\s+([\-\~\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = m.group(1), m.group(2), m.group(3)
                    targets[name][hitnum]["Tcons"] += ali
                else:
                    print("FAILED", line, pathin)
            elif line.startswith("T "):
                m = re.match("\s*(\d+)\s+([\-\w]*)\s+(\d+)\s\((\d+)\)",line[16:])
                if m:
                    start, ali, stop, size = int(m.group(1)), m.group(2), int(m.group(3)), int(m.group(4))
                    targets[name][hitnum]["Tstart"] = min(targets[name][hitnum]["Tstart"], start)
                    targets[name][hitnum]["Tstop"] = max(targets[name][hitnum]["Tstop"], stop)
                    targets[name][hitnum]["Tali"] += ali
                    targets[name][hitnum]["Tsize"] = size
                elif debug:
                    print("FAILED", line, pathin)
    # check all hits have an entry, sometimes hhblits return empty matchs
    target_names = list(targets.keys())
    for name in target_names:
        hits_to_remove = list()
        for hitnum in targets[name]:
            if (targets[name][hitnum]["Tstart"] == 1e10 or targets[name][hitnum]["Qstart"] == 1e10 or
                targets[name][hitnum]["Tstop"] == -1 or targets[name][hitnum]["Qstop"] == -2 or
                targets[name][hitnum]["Tcons"] == "" or targets[name][hitnum]["Qcons"] == ""):
                hits_to_remove.append(hitnum)
        for hitnum in hits_to_remove:
            del targets[name][hitnum]
        if len(targets[name]) == 0:
            del targets[name]
    return targets

def run_hhblits(pathquery, workdir, pathdb, parameters):
    """ run hhblits
    """
    res = dict()
    parameters_cmd = " "
    for k in parameters["hhblits_options"].keys():
        try:
            parameters_cmd += " -{} {}".format(k.strip(), parameters["hhblits_options"].get(k))
        except:
            raise ValueError("Unable to parse value for key {} in section 'hhblits_options'".format(k))

    parameters_cmd  += " "
    pathout = os.path.join(workdir, "query_hhblits.hhr")
    #if not  os.path.isfile(pathout):
    pathsco = os.path.join(workdir, "query_hhblits.scores")
    pathlog = os.path.join(workdir, "query_hhblits.log")
    hhblits_path = parameters["path"].get("hhblits")
    #hhblits_path  = ast.literal_eval('\''+parameters["path"].get("hhblits")+'\'') # safe way to parse config arguments

    command = "{} -i {} -d {} -scores {} -o {} {}".format(hhblits_path, pathquery, pathdb, pathsco, pathout, parameters_cmd)
    print(command)
    with open(pathlog, "w") as logf:
        try:
            #print("Should no have gone here, hhblits")
            a = subprocess.check_call(shlex.split(command), stdout=logf, stderr=subprocess.STDOUT)
        except:
            print("Unable to run command: {}".format(command), file=sys.stderr)
            print("Please check path to hhblits, paths to the database and configuration options")
            sys.exit(1)
    return pathout


# jackhmmer_like

def targets_jackhmmer_like(pathquery, workdir, pathdb, parameters):
    """ run and read results of hhblits
    """
    pathout, proteins = run_jackhmmer_like(pathquery, workdir, pathdb, parameters)
    targets = read_jackhmmer_like(pathout, proteins)
    #new_targets = filter_hhblits(targets, cutoff_evalue)
    return targets


def read_jackhmmer_like(pathin, kept_proteins):
    """ read hhblits results
    """
    alltargets = set()
    targets = dict()
    hitnumber = dict()
    with open(pathin) as inf:
        search_ali = False
        next_to_domain_header = False
        for line in inf:
            if line[0] == ">" and line[1] == ">":
                tmp = line[2:].strip().split()
                name = tmp[0]
                descr = " ".join(tmp[1:])
                targets.setdefault(name, dict())
                alltargets.add(name)
                search_ali = False
                next_to_domain_header = False
            elif line.starswith("  ==") and not search_ali:
                #  == domain 2  score: 66.8 bits;  conditional E-value: 4e-20
                tmp = line.strip().split()
                hitnum = int(tmp[2])
                score = float(tmp[4])
                e_value = float(tmp[8])
                targets[name].setdefault(hitnum, {"descr": descr, "Tstart":1e10, "Tstop":-1, "Tcons":"", "Tali":"", "Tsize": 0,
                                                "Qstart":1e10, "Qstop":-1, "Qcons":"", "Qali":"", "Qsize": 0,
                                "Probab":-1, "E-value":e_value, "Score":score, "Identities":-1, "Similarity":-1, "Sum_probs":-1,
                                })
                search_ali = True
                next_to_domain_header = True
            elif search_ali and next_to_domain_header:
                tmp = line.strip().split()
                query = tmp[0]
                start = tmp[1]
                seq = tmp[2].replace(".", "-").upper()
                stop = tmp[2]
                targets[name][hitnum]["Qstart"] = int(start)-1
                targets[name][hitnum]["Qstop"] = int(stop)
                targets[name][hitnum]["Qali"] = seq
                next_to_domain_header = False
            elif search_ali and line.strip().startswith(name):
                tmp = line.split()
                query = tmp[0]
                start = tmp[1]
                seq = tmp[2].replace(".", "-").upper()
                stop = tmp[2]
                targets[name][hitnum]["Tstart"] = int(start)-1
                targets[name][hitnum]["Tstop"] = int(stop)
                targets[name][hitnum]["Tali"] = seq
                search_ali = False
            elif search_ali and not line.strip().startswith(name):
                tmp = line.strip()
                cons = line.replace(" ", "-")
                targets[name][hitnum]["Tcons"] = cons
                targets[name][hitnum]["Qcons"] = cons
                
    # check all hits have an entry, sometimes hhblits return empty matchs
    target_names = list(targets.keys())
    for name in target_names:
        hits_to_remove = list()
        for hitnum in targets[name]:
            if (targets[name][hitnum]["Tstart"] == 1e10 or targets[name][hitnum]["Qstart"] == 1e10 or
                targets[name][hitnum]["Tstop"] == -1 or targets[name][hitnum]["Qstop"] == -2 or
                targets[name][hitnum]["Tcons"] == "" or targets[name][hitnum]["Qcons"] == ""):
                hits_to_remove.append(hitnum)
        for hitnum in hits_to_remove:
            del targets[name][hitnum]
        if len(targets[name]) == 0:
            del targets[name]
    return targets


def run_phmmer(pathquery, workdir, pathdb, parameters):
    """ from an initial protein sequence query a database
    """
    outpath = os.path.join(workdir, os.path.splitext(os.path.basename(pathquery))[0])
    ali_out = outpath+"_msa.fasta"
    res_out = outpath+"_res.txt"
    log_file = outpath+"_plog.txt"
    
    phmmer_path = parameters["path"].get("phmmer")
    parameters_cmd = " "
    for k in parameters["phmmer_options"].keys():
        val = parameters["phmmer_options"].get(k)
        if len(k.strip()) == 1:
            try:
                parameters_cmd += " -{} {}".format(k.strip(), val)
            except:
                raise ValueError("Unable to parse value for key {} in section 'phmmer_options'".format(k))
        else:
            try:
                parameters_cmd += " --{} {}".format(k.strip(), val)
            except:
                raise ValueError("Unable to parse value for key {} in section 'phmmer_options'".format(k))
        
    command = "{} {} -A {} -o {} {} {}".format(phmmer_path, parameters_cmd, ali_out, res_out, pathquery, pathdb)
    cmd = shlex.split(command)
    with open(log_file, "w") as logf:
        try:
            subprocess.check_call(cmd, stdout=logf, stderr=subprocess.PIPE)
        except:
            print("Error, unable to run {}".format(command), file=sys.stderr)
            sys.exit(1)
    return ali_out

def run_hmmbuild(pathquery, workdir, parameters):
    """ build an hmm out of a fasta file
    """
    hmmout = outpath+".hmm"
    log_file = outpath+"_hmmlog.txt"
    hmmbuild_path = parameters["path"].get(hmmbuild)
    
    command = "{} {} {}".format(hmmbuild_path, hmmout, pathquery)
    cmd = shlex.split(command)
    with open(log_file, "w") as logf:
        try:
            subprocess.check_call(cmd, stdout=logf, stderr=subprocess.PIPE)
        except:
            print("Error, unable to run {}".format(command), file=sys.stderr)
            sys.exit(1)
    return hmmout

def run_hmmsearch(pathquery, workdir, pathdb, parameters):
    """ run hmmsearch and keep aligned hits
    """
    outpath = os.path.join(workdir, os.path.splitext(os.path.basename(pathquery))[0])
    ali_out = outpath+"_msa.fasta"
    res_out = outpath+"_res.txt"
    log_file = outpath+"_hlog.txt"
    
    hmmsearch_path = parameters["path"].get(hmmsearch)
    parameters_cmd = " "
    for k in parameters["hmmsearch_options"].keys():
        val = parameters["hmmsearch_options"].get(k)
        if len(val) == 1:
            try:
                parameters_cmd+= " -{} {}".format(k.strip(), val)
            except:
                raise ValueError("Unable to parse value for key {} in section 'hmmsearch_options'".format(k))
        else:
            try:
                parameters_cmd += " --{} {}".format(k.strip(), val)
            except:
                raise ValueError("Unable to parse value for key {} in section 'hmmsearch_options'".format(k))
        
        
    command = "{} {} --notextw -A {} -o {} {} {}".format(hmmsearch_path, parameters_cmd, ali_out, res_out, pathquery, pathdb)
    cmd = shlex.split(command)
    with open(log_file, "w") as logf:
        try:
            subprocess.check_call(cmd, stdout=logf, stderr=subprocess.PIPE)
        except:
            print("Error, unable to run {}".format(command), file=sys.stderr)
            sys.exit(1)
    return ali_out

def read_results(path):
    """ read aligned sequences of hmmsearch in stockholm format.
    """
    proteins = dict()
    with open(path) as inf:
        msa = AlignIO.read(inf, format="stockholm")
    for record in msa:
        proteins[record.id] = str(record.seq)
    return proteins

def read_input(path):
    """ read initial input
    """
    proteins = set()
    with open(path) as inf:
        for line in inf:
            if line[0] == ">":
                proteins.add(line[1:].strip().split()[0])
    return proteins

def filter_hmmsearch_alignment(proteins, name, workdir, init_seq, parameters):
    """ filter results based on hit coverage of initial input and on regexp presence/absence
    """
    init_length = len(init_seq)
    cutoff_identity = parameters["jackhmmer_like_options"].getfloat("identity")
    cutoff_coverage = parameters["jackhmmer_like_options"].getfloat("coverage")
    kept_regx = parameters["jackhmmer_like_options"].get("kept_reg")
    rm_regx = parameters["jackhmmer_like_options"].get("remove_reg")
    
    print(workdir, name)
    outpath = os.path.join(workdir, name+"_filteredmsa.fasta")
    query_msa = proteins[name]
    
    if rm_gx:
        kept_proteins = [name]
        r_remove = re.compile(rm_gx)
        for prot in proteins:
            seq_msa = proteins[prot]
            seq = seq_msa.replace("-", "").replace(".", "").upper()
            # regexp match
            match = r_remove.search(seq)
            if not match:
                kept_proteins.append(prot)
    else:
        kept_proteins = proteins[:]
    
    if kept_reg:
        r_kept = re.compile(kept_reg)
        kept = set()
        with open(outpath, "w") as outf:
            for prot in proteins:
                seq_msa = proteins[prot]
                seq = seq_msa.replace("-", "").replace(".", "").upper()
                # regexp match
                match = r_kept.search(seq)
                if match:
                    identity, coverage = compute_id_cov(query_msa, seq_msa, cutoff_coverage, cutoff_identity)
                    if identity and coverage:
                        outf.write(">{}\n{}\n".format(prot, seq_msa))
                        kept.add(prot.split()[0])
    return kept, outpath

def compute_id_cov(query_msa, seq_msa, cutoff_coverage, cutoff_identity):
    """ check if identity and coverage are sufficient
    """
    size, ident, cov = 0, 0, 0
    
    for i in range(len(query_msa)):
        if query_msa[i] != "-":
            if query_msa[i] == seq_msa[i]:
                ident += 1
            if seq_msa[i] != "-":
                cov += 1
            size += 1
    ident = (ident / size) * 100
    cov  = (cov / size) * 100
    coverage = cov > cutoff_coverage
    identity = ident > cutoff_identity
            
    return identity, coverage 

def read_results_and_filter(ali_results, name, n, workdir, init_seq, parameters):
    """ apply read and filter to results
    """
    # read results
    hit_proteins = read_results(ali_results)
    # filter
    res_proteins, filtered_ali = filter_hmmsearch_alignment(hit_proteins, name+"_iter_{}".format(n), workdir, init_seq, parameters)
    filtered_hmm = os.path.join(workdir, name+"_iter_{}_res.txt".format(n))
    return res_proteins, filtered_ali, filtered_hmm

def compute_init_values(path):
    """ compute initial mean sequence compute_mean_length
    """
    length = 0
    fasta = dict()
    with open(path) as inf:
        for line in inf:
            if line[0] == ">":
                prot = line[1:-1]
                fasta[prot] = ""
            else:
                fasta[prot] += line.strip().replace("-", "").replace(".", "")
    if len(fasta) > 1:
        raise ValueError("Multiple proteins found in initial fasta file {}".format(path))
    
    seq = fasta[prot]
    length = len(seq)
    return seq, length


def count_differences(query_proteins, targets, cov):
    notfound = 0
    overlapping = set()
    for prot in query_proteins:
        tmp = prot.split("/")
        name = tmp[0]
        start, stop = tmp[1].split("-")
        init_start, init_stop = int(start)-1, int(stop)
        if name in targets:
            found = False
            for new_start, new_stop in targets[name]:
                start = max(init_start, new_start)
                stop = min(init_stop, new_stop)
                diff = stop - start 
                if diff > 0:
                    c = diff / max(new_stop-new_start, init_stop-init_start)
                    if c > cov:
                        found = True
                        overlapping.add(prot)
            if not found:
                notfound += 1
        else:
            notfound += 1
    return notfound, overlapping

def check_set_differences(new_proteins, prev_proteins, cov=0.9):
    """ check hits, count number of new sequences and dropped sequences
    """
    init = dict()
    overlapping = set()
    for prot in prev_proteins:
        tmp = prot.split("/")
        start, stop = tmp[1].split("-")
        init.setdefault(tmp[0], list()).append((int(start)-1, int(stop)))
    
    new = dict()
    for prot in new_proteins:
        tmp = prot.split("/")
        start, stop = tmp[1].split("-")
        new.setdefault(tmp[0], list()).append((int(start)-1, int(stop)))
        
    nbnew, overlapping = count_differences(new_proteins, init, cov)
    dropped, _= count_differences(prev_proteins, new, cov)
            
    print(nbnew, dropped, len(new_proteins), len(prev_proteins), len(overlapping))
    inboth = len(overlapping)
    if len(overlapping) == len(new_proteins):
        return True
    return False


def run_jackhmmer_like(pathquery, workdir, pathdb, parameters):
    """ run jackhmmer_like : phmmer, hmmbuild, hmmsearch
    """
   
    name = os.path.splitext(os.path.basename(pathquery))[0]
    init_proteins = [name]
    stop = False
    init_size = 1
    init_seq, init_length = compute_init_values(pathquery)
    #print("Coverage threshold to add new protein: [{}; {}]".format(init_length-init_length*(1-params.coverage),
                                                                   #init_length+init_length*(1-params.coverage)))
    inpath = pathquery
    if init_size == 0:
        print("Error, no protein found in input fasta file, {}".format(pathquery), file=sys.stderr)
        sys.exit(1)

    coverage = parameters["jackhmmer_like_options"].getfloat("coverage")
    # run phmmer
    ali_results = run_phmmer(inpath, workdir, pathdb, parameters)
    res_proteins, filtered_ali_results, filtered_hmm_results = \
                            read_results_and_filter(ali_results, name, 0, workdir, init_seq, parameters)
    
    if len(res_proteins.intersection(init_proteins)) == len(res_proteins):
        # no new proteins
        stop = True
    else:
        init_proteins = res_proteins
        inpath = filtered_ali_results

    # run hmmsearch
    if not stop:
        niter = parameters["jackhmmer_like_options"].getint("nb_iterations") + 1
        for n in range(1, niter):
            outpath = os.path.join(params.workdir, name+"_iter_{}".format(n))

            # convert fasta inputfile to hmm file
            hmmpath = run_hmmbuild(inpath, outpath)

            # run hmmsearch
            ali_results, hmm_results = run_hmmsearch(hmmpath, outpath, db_target, parameters)
            res_proteins, filtered_ali_results, filtered_hmm_results =\
                            read_results_and_filter(ali_results, name, params.workdir, n, 
                                                    init_seq, parameters)

            if check_set_differences(res_proteins, init_proteins):
                # no new proteins
                break
            if len(res_proteins) == 0:
                print("Error, no protein found", file=sys.stderr)
                break
            else:
                # new iteration
                init_proteins = res_proteins
                inpath = filtered_ali_results

    return filtered_hmm_results, res_proteins
