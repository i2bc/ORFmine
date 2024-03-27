#!/usr/bin/env python
""" package grouping utility functions around the HCA methodology
"""
import os, gzip, sys

def read_hclustdict(path, cutoff_coil=0):
    """ read the dictionary of hydrocluster
    """
    amas = set()
    with open(path) as inf:
        header = inf.readline()
        for line in inf:
            tmp = line.strip().split(";")
            val = float(tmp[4])
            if val > cutoff_coil:
                amas.add(tmp[1])
    return list(amas)

def read_hcaresults(path, domains=True, clusters=True):
    """ read hca domains and associated clusters
    """
    if os.path.splitext(path)[1] in [".gz", ".gzip"]:
        return _read_hcaresults_gzip(path, domains=domains, clusters=clusters)
    else:
        return _read_hcaresults(path, domains=domains, clusters=clusters)
    
def _read_hcaresults(path, domains=True, clusters=True):
    data_domain = dict()
    data_cluster = dict()
    data = dict()
    with open(path) as inf:
        for line in inf:
            tmp = line.split()
            if line[0] == ">":
                prot, size = tmp[0][1:], tmp[1]
                size = int(size)
                if prot not in data_domain:
                    data_cluster[prot] = \
                        dict([("size", size), ("data", list())])
                    data_domain[prot] = \
                        dict([("size", size), ("data", list())])
                else:
                    print("Error, the protein {} has already been parsed"\
                        .format(prot), file=sys.stderr)                
            elif domains and tmp[0] == "domain":
                tmp = line.split()
                data_domain[prot]["data"]\
                    .append((int(tmp[1])-1, int(tmp[2])))
            elif clusters and tmp[0] == "cluster":
                data_cluster[prot]["data"]\
                    .append((int(tmp[1])-1, int(tmp[2]), tmp[3]))
    if domains and clusters:
        for prot in data_domain:
            data[prot] = dict([("size", data_domain[prot]["size"]), 
                               ("domains", dict()), ("clusters", list())])
            for i, (clustart, clustop, clust) in enumerate(data_cluster[prot]["data"]):
                data[prot]["clusters"].append((clustart, clustop, clust))
                for domstart, domstop in data_domain[prot]["data"]:
                    if domstart <= clustart and clustop <= domstop:
                        data[prot]["domains"]\
                            .setdefault((domstart, domstop), list())\
                            .append(i)
    elif domains:
        data = data_domain
    else:
        data = data_cluster
    return data
                
def _read_hcaresults_gzip(path, domains=True, clusters=True):
    data_domain = dict()
    data_cluster = dict()
    data = dict()
    with gzip.open(path, "rt", encoding="UTF-8") as inf:
        for line in inf:
            tmp = line.split()
            if line[0] == ">":
                prot, size = tmp[0][1:], tmp[1]
                size = int(size)
                if prot not in data_domain:
                    data_cluster[prot] = \
                        dict([("size", size), ("data", list())])
                    data_domain[prot] = \
                        dict([("size", size), ("data", list())])
                else:
                    print("Error, the protein {} has already been parsed"\
                        .format(prot), file=sys.stderr)
                
            elif domains and tmp[0] == "domain":
                tmp = line.split()
                data_domain[prot]["data"]\
                    .append((int(tmp[1])-1, int(tmp[2])))
            elif clusters and tmp[0] == "cluster":
                data_cluster[prot]["data"]\
                    .append((int(tmp[1])-1, int(tmp[2]), tmp[3]))
    if domains and clusters:
        for prot in data_domain:
            data[prot] = dict([("size", data_domain[prot]["size"]), 
                               ("domains", dict()), ("clusters", list())])
            for i, (clustart, clustop, clust) in enumerate(data_cluster[prot]["data"]):
                data[prot]["clusters"].append((clustart, clustop, clust))
                for domstart, domstop in data_domain[prot]["data"]:
                    if domstart <= clustart and clustop <= domstop:
                        data[prot]["domains"]\
                            .setdefault((domstart, domstop), list())\
                            .append(i)
    elif domains:
        data = data_domain
    else:
        data = data_cluster
    return data
                
