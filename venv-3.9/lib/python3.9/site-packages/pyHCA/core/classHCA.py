#!/usr/bin/env python
""" The classHCA is an ensemble of class to describe hydrophobic clusters and
domain created through hydrophobic clusters.
"""

import numpy as np
import scipy.stats as st

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "tristan.bitard-feildel [you know what] impmc.upmc.fr"
__institute__ = "IMPMC, UPMC"



class Seq(object):
    """ A Sequence object """
    def __init__(self, name,  description, seq, length):
        """ hold the name, the description and the amino acid sequence of a fasta entry
        """
        self.length = length
        self.name = name
        self.descr = description
        self.seq = seq

PDB_SEQRES_FREQ = {
"A":       0.0859929227395,
"C":       0.0206946259759,
"E":       0.0648133592624,
"D":       0.0549383608817,
"G":       0.0830137574641,
"F":       0.0378506807406,
"I":       0.0547752297123,
"H":       0.0260798136952,
"K":       0.057416225372,
"M":       0.0230113512204,
"L":       0.0871762434274,
"N":       0.0408935875728,
"Q":       0.0364996548365,
"P":       0.045826092199,
"S":       0.0602061283907,
"R":       0.0504294572634,
"T":       0.05511062539,
"W":       0.0131711008412,
"V":       0.068839479608,
"Y":       0.0332613034068,
}

class HydroCluster(object):
    """ the class definning the amas of hydrophobic clusters
    """
    #__slots__  = ["__start", "__stop", "__hydro_cluster"]
    def __init__(self, start, stop, hydro_cluster):
        self.__start = start
        self.__stop = stop
        #self.__hydro_cluster = hydro_cluster[:]
        self.__str_hydro_cluster = hydro_cluster[:]
        self.__hydro_cluster = [int(val) for val in hydro_cluster]

    @property
    def start(self):
        return self.__start
    
    @property
    def stop(self):
        return self.__stop
    
    @property
    def hydro_cluster(self):
        return self.__hydro_cluster
    
    def __str__(self):
        """ write position, indexes are inclusive and start from 1
        """
        return "cluster\t{}\t{}\t{}".format(self.__start+1, self.__stop, 
                                                           self.__str_hydro_cluster)
    def add_offset(self, offset):
        """ add an offset to the hydrophobic clusters
        """
        self.__start += offset
        self.__stop += offset
        
    def get(self, method):
        if method == "all":
            return self.__start, self.__stop, self.__hydro_cluster
        elif hasattr(self, method):
            return getattr(self, method)
        else:
            raise ValueError("Unknown argument for function "
                             "HydroCluster.get({}".format(method))

def compute_disstat(start, stop, clusters):
    """ compute the domain pvalue based on length, and clusters
    """
    a, b, c = -10, 9, 10
    size = stop - start
    N = size # 2 * size
    nb_cluster = len(clusters)
    cov = np.zeros(size)
    cov.fill(a)
    for clust in clusters:
        if len(clust.hydro_cluster) > 2:
            cov[clust.start - start: clust.stop - start] = \
                [b if clust.hydro_cluster[i] == 1 else c for i in range(len(clust.hydro_cluster))]
    score = sum(cov) / size
    score = score

    # inverse gaussian parameters fitted from disprot v7 sequence scores
    # PDB fitted paramters (0.0028441615833769331, -36.33441401336944, 0.10425345380333628)
    # DisProt fitted paramters (0.14091823798877751, -11.362093616044803, 0.54142167247756956)
    
    # return st.recipinvgauss.sf(score, *(0.14091823798877751, -11.362093616044803, 0.54142167247756956))
    # PDB (0.0023527450356064881, -40.277311463079315, 0.095643136597401535) 
    # DisProt (0.17521210606656379, -11.060272542998465, 0.6129194848871018)        
    params = (0.17521210606656379, -11.060272542998465, 0.6129194848871018)  # hca2
    #params = (0.10859453489976689, -14.651939851385755, 1.1206990484737327) # hca1
    return score, st.recipinvgauss.sf(score, *params)

class DomHCA(object):
    """ the class definning the domain delineated by the HCA segmentation
    """
    __slots__ = ["__start", "__stop", "__pvalue", "__score"]
 
    def __init__(self, start, stop, clusters): #, list_of_hcclusters):
        self.__start = start
        self.__stop = stop
        self.__score = np.nan 
        self.__pvalue = self._compute_pvalue(clusters)
        if stop - start < 30:
            self.__pvalue = np.nan
        
        #self.__clusters = list_of_hcclusters[:]
    @property
    def start(self):
        return self.__start
    
    @property
    def pvalue(self):
        return self.__pvalue
    
    @property
    def score(self):
        return self.__score
    
    @property
    def stop(self):
        return self.__stop
    
    #@property
    #def clusters(self):
        #return self.__clusters

    def add_offset(self, offset):
        """ add an offset to start and stop and all hydrophobic clusters inside
        the domain
        """
        self.__start += offset
        self.__stop += offset
        #for hclust in self.__clusters:
            #hclust.add_offset(offset)

    #def __iter__(self):
        #return iter(self.__clusters)
            
    def __str__(self):
        domain = "domain\t{}\t{}\t{}\t{}".format(self.__start+1, self.__stop, self.__pvalue, self.__score)
        return domain

    #def _compute_pvalue2(self, clusters):
    #    """ compute the domain pvalue based on length, and clusters
    #    """
    #    coef =  [14.76376808, -34.28491775,   6.63229591] #[ 14.53120892, -34.21923594,   6.83311063]
    #    intercept = -4.67632801

    #    size = self.stop - self.start
    #    Hinside = 0
    #    Pinside = 0
    #    for clust in clusters:
    #        if len(clust.hydro_cluster) > 2:
    #            for i in range(len(clust.hydro_cluster)):
    #                if clust.hydro_cluster[i] == 1:
    #                    Hinside += 1
    #                else:
    #                    Pinside += 1
    #    outside = size - Hinside - Pinside

    #    x = np.asarray([outside, Hinside, Pinside], dtype=float)
    #    x /= size
    #    
    #    V = x * coef
    #    U = V.sum() + intercept
    #    U *= -1
    #    A = np.exp(U)+1
    #    P = np.reciprocal(A)
    #    self.__score = U

    #    return P 

    def _compute_pvalue(self, clusters):
        score, pval = compute_disstat(self.start, self.stop, clusters)
        self.__score = score
        return pval


