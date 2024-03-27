#!/usr/bin/env python
""" compute predicted disorder per residue based on HCA profile
"""

import os, sys, argparse
import numpy as np
import pandas as pd
try:
    from sklearn.externals import joblib
except:
    import joblib
import lightgbm
from pyHCA.core.disorderHCA_data import model_descriptors, key_order, thresholds
from pyHCA import HCA 
from pyHCA.core.seq_util import transform_seq, check_if_msa

def prepare_sequence(seq):
    """ compute score for a given sequence
    """
    size = len(seq)
    hcaprot = HCA(seq=seq)
    hcaclusters = hcaprot.get_clusters()
    hcadomains  = hcaprot.get_domains()
    return  hcadomains, hcaclusters

def get_not_tagged(positions):
    not_tagged = {}
    idx = np.where(positions == 0)[0]
    not_dom_idx = 0
    if len(idx) > 0:
        start = idx[0]
        cur = start
        if len(idx) == 1:
            not_tagged[start] = (not_dom_idx, 1, start, start+1)
        else:
            for i in range(1, len(idx)):
                if cur+1 == idx[i]:
                    cur = idx[i]
                else:
                    stop = cur
                    for c in range(start, stop+1):
                        not_tagged[c] = (not_dom_idx, stop - start + 1, start, stop+1)
                    start = idx[i]
                    cur = idx[i]
                    not_dom_idx += 1
                if i+1 == len(idx):
                    stop = idx[i]
                    for c in range(start, stop+1):
                        not_tagged[c] = (not_dom_idx, stop - start + 1, start, stop+1)
    return not_tagged

def compute_features(seq, 
                     domains, clusters, 
                     feature_range_to_name, descriptors,
                     verbose=False):
    output = {}
    domain_area = list()
    domain_size = dict()
    domain_pos = np.zeros(len(seq))
    viewed = set()
    for dom_id, dom in enumerate(domains):
        domain_pos[dom.start: dom.stop] = 1
        for i in range(dom.start, dom.stop):
            domain_size[i] = (dom_id, dom.stop - dom.start + 1)
            if dom_id not in viewed:
                domain_area.append((dom.start, dom.stop, dom_id, True))
                viewed.add(dom_id)
    
    # get inter domain regions inside list of domains
    viewed = set()
    not_domains = get_not_tagged(domain_pos)
        
    for res in not_domains:
        not_dom_idx, size, start, stop =  not_domains[res]
        if not_dom_idx not in viewed:
            domain_area.append((start, stop, not_dom_idx, False))
            viewed.add(not_dom_idx)
            
    domain_area.sort()
    
    cluster_area = list()
    cluster_pos = np.zeros(len(seq))
    cluster_layout = np.zeros(len(seq))
    viewed = set()
    for cls_id, cls in enumerate(clusters):
        #cluster_pos[cls.start: cls.stop] = 1
        cluster_area.append((cls.start, cls.stop, cls_id, True))
        k = 0
        for i in range(cls.start, cls.stop):
            cluster_layout[cls.start: cls.stop] = 1
            if cls.hydro_cluster[k] == 1:
                cluster_pos[i] = 1
            k+=1
    
    # get inter cluster regions inside list of domains
    viewed = set()
    not_clusters = get_not_tagged(cluster_layout)
    for res in not_clusters:
        not_cls_idx, size, start, stop =  not_clusters[res]
        if not_cls_idx not in viewed:
            cluster_area.append((start, stop, not_cls_idx, False))
            viewed.add(not_cls_idx)
    
    cluster_area.sort()
    
    # get physico-chemical properties           
    seq_descriptors = {} #np.zeros((len(descriptors), len(seq)))
    for i, d in enumerate(descriptors):
        seq_descriptors[d] = np.zeros(len(seq))
        for j, aa in enumerate(seq):
            if aa in descriptors[d]:
                seq_descriptors[d][j] = descriptors[d][aa]
                
    # reverse, get seg idx for cluster and domain
    aa_to_dom = {}
    for idx, (start, stop, dom_idx, is_dom) in enumerate(domain_area):
        for i in range(start, stop):
            aa_to_dom[i] = idx
    
    aa_to_cls = {}
    for idx, (start, stop, cls_idx, is_cls) in enumerate(cluster_area):
        for i in range(start, stop):
            aa_to_cls[i] = idx
    
    # get domain properties
    for i in range(len(seq)):
        
        # disordered or ordered
        cls_idx_in_list = aa_to_cls[i]
        dom_idx_in_list = aa_to_dom[i]
        start_cls, stop_cls, cls_idx, is_cls = cluster_area[cls_idx_in_list]
        start_dom, stop_dom, dom_idx, is_dom = domain_area[dom_idx_in_list]
        
        dom_size = stop_dom - start_dom
        output.setdefault("domain_size", []).append(np.log(dom_size))
        output.setdefault("perc_hydro_domain", []).append(cluster_pos[start_dom: stop_dom].sum() / dom_size)

        # distance to other domains
        if dom_idx_in_list > 1:
            prev_dom_idx_in_list  = dom_idx_in_list - 2
            assert prev_dom_idx_in_list >= 0
            prev_start_dom, prev_stop_dom, prev_dom_idx, prev_is_dom = domain_area[prev_dom_idx_in_list]
            assert prev_is_dom == is_dom
            d_prev = start_dom - prev_start_dom
        else:
            d_prev = start_dom # distance to start
            
        if dom_idx_in_list < len(domain_area) - 2:
            next_dom_idx_in_list = dom_idx_in_list + 2
            next_start_dom, next_stop_dom, next_dom_idx, next_is_dom = domain_area[next_dom_idx_in_list]
            assert next_is_dom == is_dom
            d_next = start_dom - next_start_dom
        else:
            d_next = len(seq) - stop_dom + 1 # distance to start
        output.setdefault("d_prev_domain", []).append(d_prev)
        output.setdefault("d_next_domain", []).append(d_next)
        
        # aa index properties
        for r in feature_range_to_name:
            for d in feature_range_to_name[r]:
                val = descriptors[d].get(aa, 0)
                if r == "pm16":
                    w_start, w_stop = -16, 16
                    s = sum(seq_descriptors[d][max(0, i+w_start): min(i+1+w_stop, len(seq)+1)])
                    v = s / (min(i+w_stop+1, len(seq)+1) - max(0, i+w_start))
                    output.setdefault("{}_pm16".format(d), []).append(v)
                elif r == "cluster":
                    s = sum(seq_descriptors[d][max(0, i-start_cls): min(i+1+stop_cls, len(seq)+1)])
                    v = s / (min(i+stop_cls+1, len(seq)+1) - max(0, i-start_cls))
                    output.setdefault("{}_cluster".format(d), []).append(v)
                elif r == "domain":
                    s = sum(seq_descriptors[d][max(0, i-start_dom): min(i+1+stop_dom, len(seq)+1)])
                    v = s / (min(i+stop_dom+1, len(seq)+1) - max(0, i-start_dom))
                    output.setdefault("{}_domain".format(d), []).append(v)
                
            
    return output

def get_params():
    """ get command line ArgumentParser
    """
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "disorder"))
    parser.add_argument("-i", action="store", dest="fastafile", help="the fasta file", required=True)
    parser.add_argument("-o", action="store", dest="outputfile", help="output file in svg format", required=True)
    parser.add_argument("-m", action="store", dest="model", help="model to use", required=True)
    parser.add_argument("--verbose", action="store_true", dest="verbose", help="print information")
    params = parser.parse_args()
    
    return params

allowed_models = ["allowDSoverlap.h5", "noDSoverlap.h5"]

model_to_features = {
    "allowDSoverlap.h5":
            ['RACS820101_cluster',
             'FAUJ880111_domain',
             'YUTK870104_domain',
             'QIAN880128_cluster',
             'BAEK050101_pm16',
             'YUTK870104_cluster',
             'RACS820107_cluster',
             'BASU050102_pm16',
             'SUYM030101_pm16',
             'LEWP710101_domain',
             'DAYM780201_cluster',
             'SNEP660103_cluster',
             'GEOR030102_cluster',
             'WILM950103_cluster',
             'FINA910103_cluster'],
    "noDSoverlap.h5":
             ['LEWP710101_domain',
             'YUTK870104_cluster',
             'BASU050102_pm16',
             'WILM950103_domain',
             'SUYM030101_domain',
             'WILM950103_cluster',
             'QIAN880128_cluster',
             'DAYM780201_cluster',
             'GEOR030102_cluster',
             'RACS820107_cluster',
             'SUYM030101_pm16'],
}
# commont to both: d_next_domain, d_prev_domain, domain_size, perc_hydro_domain


def main():
    params = get_params()
    
    from pyHCA.core.ioHCA import read_multifasta_it

    model_name = os.path.basename(params.model)
    if model_name not in allowed_models:
        print("Error, unrecognized models {}".format(model_name), file=sys.stderr)
        print("Models can be downloaded from the github repository: https://github.com/T-B-F/pyHCA/data/", file=sys.stderr)
        sys.exit(1)
    # load model
    lgb_classifier = lightgbm.Booster(model_file=params.model)


    # get feature to compute
    feature_names_and_range = model_to_features[model_name]
    feature_range_to_name = {}
    for name in feature_names_and_range:
        name, r = name.split("_")
        feature_range_to_name.setdefault(r, []).append(name)

    feature_names = [name.split("_")[0] for name in feature_names_and_range]

    # get descriptors values
    descriptors = model_descriptors[model_name]

    keys = key_order[model_name]

    # prepare features and run prediction 
    if params.verbose:
        print("Processing")
    with open(params.outputfile, "w") as outf:
        for prot, sequence in read_multifasta_it(params.fastafile):
            if params.verbose:
                print(prot)
                print(sequence.seq)
            seq = str(sequence.seq).upper()
            domains, clusters = prepare_sequence(seq)
            features = compute_features(seq, domains, clusters, feature_range_to_name, descriptors, params.verbose)
            features = pd.DataFrame(features)[keys]
            probas = 1 - lgb_classifier.predict(features)
            classes = (probas > thresholds[model_name]).astype(int)
            outf.write(">{}\n".format(prot))
            cnt = 0
            for i in range(len(seq)):
                if seq[i] == "X":
                    outf.write("{} {} nan nan\n".format(i+1, seq[i]))
                else:
                    outf.write("{} {} {:.3f} {}\n".format(i+1, seq[i],  probas[cnt], classes[cnt]))
                    cnt += 1

    sys.exit(0)
    
if __name__ == "__main__":
    main()

