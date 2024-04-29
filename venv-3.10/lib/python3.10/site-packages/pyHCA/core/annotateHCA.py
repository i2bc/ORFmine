#!/usr/bin/env python
""" The annotateHCA file is an ensemble of functions to annotate amino acid 
sequences with domains defined by SegHCA or hydrophobic clusters
"""

import os, sys, argparse, string, re
from pyHCA.core.ioHCA import read_multifasta_it, write_annotHCA
from pyHCA.core.classHCA import HydroCluster, DomHCA
from pyHCA.core.classHCA import compute_disstat
from pyHCA.core.seq_util import six_frames
from Bio import Seq
import numpy as np

__author__ = "Tristan Bitard-Feildel, Guillem Faure, Isabelle Callebaut"
__licence__= "MIT"
__version__ = 0.1
__email__ = "t.bitard.feildel [you know what] uni-muenster.de"
__institute__ = "Institute for Evolution and Biodiversity, Muenster Germany"


def _transformSequence(seq, low_complexity=[], hydrophobe="YIMLFWVC"):
    """ transforme amino acid string chain into kind of HCA code as 1-> YIMLFWV  P-> P  0-> other 
    
    Parameters
    ----------
    seq ; string
        amino acid sequence
    low_complexity ; list
        list of positions of low complexity
    hydrophobe: string
        define hydrophobic residue used by the HCA method
    Return
    ------
    setrans : string
        code of the sequence as a string chain such as 1-> YIMLFWV  P-> P  0-> other 
    """
    #===========================================================================
    # get the binary code 
    # 1-> YIMLFWV 
    # P-> P
    # * -> * 
    # 0-> other 
    #===========================================================================
    seqtrans = ""
    for ite, aa in enumerate(seq):
        AA = aa.upper()
        if low_complexity and ite+1 in low_complexity:
            seqtrans += "0"
        elif AA  in hydrophobe:
            seqtrans += "1"
        elif AA == "P":
            seqtrans += "P"
        elif AA == "*": #stop 
            seqtrans += "*"
        else:
            seqtrans += "0"
    return seqtrans

def _getPosition(match, start):
    """
    """
    lite = []
    for ite, ii in enumerate(match):
        if ii == "1":
            lite.append(ite+start)
    return lite

def _manageReplacement(seqtrans):
    """
    brief this step delete small cluster 1 and 11 by connectivity distance or
    Proline/stop breakers
    param seqtrans is a string of sequence code 1-> hydrophobe P-> Proline *->* 
    and 0-> other
    return seqtrans is a string without small cluster there are replaced by 0 or 00
    """
    
    lsmall_cluster = [] # contains the position of small cluster 1 or 11
    
    #===========================================================================
    # 1/ small cluster 1 and 11
    #===========================================================================
    # #===========================================================================
    #===========================================================================
        
    #===========================================================================
    # we want to delete 000010000 and 0000110000
    # we proceed in 2 steps because of overlapping sometime like 000001000001100000000
    #===========================================================================
    for i in re.finditer(r"(?=(0{4}1{1,2}0{4}))", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        
        lsmall_cluster += _getPosition(match, start)
        
        #=======================================================================
        # replacement
        #=======================================================================
        seqtrans = seqtrans[:start] + "0"*lmatch + seqtrans[stop:]
    
    #===========================================================================
    # 2/ Proline surrounded
    #===========================================================================
    # #===========================================================================
    #===========================================================================
    
    #===========================================================================
    # we want to delete every 1 or 11 surounded by 2 prolines or 2 *
    #===========================================================================
    # no working
    #for i in re.finditer(r"(?=([P0,\*0]{0,3}1{1,2}0{0,3}[P,\*]))", seqtrans):
    for i in re.finditer(r"(?=(P0{0,3}1{1,2}0{0,3}[P,\*]))", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        
        lsmall_cluster += _getPosition(match, start)
        
        #=======================================================================
        # replacement
        #=======================================================================
        seqtrans = seqtrans[:start+1] + "0"*(lmatch-2) + seqtrans[stop-1:]
    
    #===========================================================================
    # we want to delete every 1 or 11 surounded by 1 proline and 0000
    #===========================================================================
    # PX0000 or *XOOO
    #for i in re.finditer(r"(?=([P0,\*]{0,3}1{1,2}0{4}))", seqtrans): #not working
    for i in re.finditer(r"(?=(P0{0,3}1{1,2}0{4}))", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        
        lsmall_cluster += _getPosition(match, start)
        
        #=======================================================================
        # replacement
        #=======================================================================
        seqtrans = seqtrans[:start+1] + "0"*(lmatch-5) + seqtrans[stop-4:]
        
    # 0000XP 
    for i in re.finditer(r"(?=(0{4}1{1,2}0{0,3}P))", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        
        lsmall_cluster += _getPosition(match, start)
        
        #=======================================================================
        # replacement
        #=======================================================================
        seqtrans = seqtrans[:start+4] + "0"*(lmatch-5) + seqtrans[stop-1:]
    
    #===========================================================================
    # remove begin + ruptor    
    #===========================================================================
    # P ruptor
    for i in re.finditer(r"(?=^(0{0,3}1{1,2}0{0,3}P))", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        seqtrans = "0"*(lmatch-1) + seqtrans[stop-1:]
    # 0000 ruptor
    for i in re.finditer(r"(?=^(0{0,3}1{1,2}0{4}))", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        seqtrans = "0"*(lmatch-1) + seqtrans[stop-4:]
    
    #===========================================================================
    # remove ruptor+end    
    #===========================================================================
    # 0000 ruptor
    for i in re.finditer(r"(?=(0{4}1{1,2}0{0,3})$)", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        seqtrans = seqtrans[:start]+"0"*(lmatch) 
    
    # P ruptor
    for i in re.finditer(r"(?=(P0{0,3}1{1,2}0{0,3})$)", seqtrans):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        seqtrans = seqtrans[:start+1]+"0"*(lmatch-1) 
    
    # print lsmall_cluster
        
    return seqtrans, lsmall_cluster
    
def _getAmas(seqtrans):
    """ the old _getAmas seems to be buggy, this new implementation is faster
    and safer
    """
    new_seq = seqtrans[:]
    seq_buffer, amas = [], []
    list_amas = []
    start = -1
    for i, aa in enumerate(seqtrans):
        # if proline or stop check if existing amas and add it to the list
        if aa == "P" or aa == "*":
            # replace hydrophilic buffer by Ruptor in the new sequence
            if seq_buffer:
                new_seq = new_seq[:i-len(seq_buffer)] + "R" * len(seq_buffer) + new_seq[i:]
                seq_buffer = []
            new_seq = new_seq[:i] + "R"  + new_seq[i+1:]
            if amas:
                #list_amas.append([start, "".join(amas)])
                hyclust = HydroCluster(start, start+len(amas), "".join(amas))
                list_amas.append(hyclust)
                amas = []
        elif aa == "0":
            # add hydrophilic to buffer
            seq_buffer.append(aa)
        elif aa == "1" and amas == []:
            # replace hydrophilic buffer by Ruptor in the new sequence
            if seq_buffer:
                new_seq = new_seq[:i-len(seq_buffer)] + "R" * len(seq_buffer) + new_seq[i:]
                seq_buffer = []
            start = i
            amas.append(aa)
        elif aa == "1" and amas != []:
            # should we create a new amas or not ? check the length of the buffer
            if len(seq_buffer) > 3:
                #list_amas.append([start, "".join(amas)])
                hyclust = HydroCluster(start, start+len(amas), "".join(amas))
                list_amas.append(hyclust)
                new_seq = new_seq[:i-len(seq_buffer)] + "R" * len(seq_buffer) + new_seq[i:]
                start = i
                amas, seq_buffer = [], []
            amas.extend(seq_buffer)
            amas.append(aa)
            seq_buffer = []
        else:
            if amas != []:
                hyclust = HydroCluster(start, start+len(amas), "".join(amas))
                list_amas.append(hyclust)
                #list_amas.append([start, "".join(amas)])
            new_seq = new_seq[:i-len(seq_buffer)] + "R" * len(seq_buffer) + new_seq[i:]
            seq_buffer, amas = [], []
    # don't forget last one
    if amas:
        #list_amas.append([start, "".join(amas)])
        hyclust = HydroCluster(start, start+len(amas), "".join(amas))
        list_amas.append(hyclust)
    if seq_buffer:
        new_seq = new_seq[:i-len(seq_buffer)+1] + "R" * len(seq_buffer)
        
    #list_amas_sorted = sorted(list_amas)
    list_amas_sorted = sorted(list_amas, key=lambda x: x.start)
    return list_amas_sorted, new_seq
    

def _removeSmallAmas(list_amas, seq):
    """ this step remove the small cluster 1 or 11 if any big cluster are close 
    before or after 7 amino acids. seq is modified this cluster is replace by R 
    or RR. A new list is generated containing each amas that we kept.
    
    Parameters
    ----------
    list_amas: list
        is a list [position, amas]
    seq: string 
        R -> ruptor(Proline, 0{4 or more}) 1 -> hydrophobe 0 -> the other
    
    Return
    -------
    keepCluster : list
        list of kept clusters
    newseq: string
        a list R -> ruptor(Proline, stop, 0{4 or more}, 1 and 11 close a big cluster) 
        1 -> hydrophobe 0 -> the other
    """
    
    keepCluster = []
    newseq = seq[:]
    for ite, hyclust in enumerate(list_amas):
        #pos, amas = pos_amas
        pos = hyclust.start
        amas = hyclust.hydro_cluster
        #hyclust = HydroCluster(pos, pos+len(amas), amas)
        add = False
        # only amas 1 or 11
        if len(amas) > 2:
            keepCluster.append(hyclust)
            continue
        
        # # is there a big amas close?
        # can we search after?
        #if ite != len(list_amas)+1:
            ## after?  ite+1 we search 1 amas to n after the current
            #for next_pos_amas in list_amas[ite+1:]:                
                #next_pos = next_pos_amas.start 
                #next_amas = next_pos_amas.hydro_cluster
                ## next_pos is the begin of the next amas pos+len(amas) is 
                ## the end of the current amas
                #distance = next_pos - pos + len(amas) 
                #if distance > 7:
                    #break
                #if len(next_amas) > 2:
                    #keepCluster.append(hyclust)
                    #add = True
                    #break
            
        # can we search before ?
        # we do not add if we already find a big cluster after
        #if ite != 0 and add == False:
            ## before? [::-1] -> reverse the list_amas len(list_amas)-ite begin 
            ## at the first previous amas
            #for prev_pos_amas in list_amas[::-1][len(list_amas)-ite:]:  
                #prev_pos = prev_pos_amas.start
                #prev_amas = prev_pos_amas.hydro_cluster
                ## pos is the begining of the current amas 
                ## prev_pos+len(prev_amas) is the end of the previous amas
                #distance = pos - (prev_pos+len(prev_amas)) 
                #if distance > 7:
                    #break
                #if len(prev_amas) > 2:
                    #keepCluster.append(hyclust)
                    #add = True
                    #break
        # small amas is transform in RUPTOR
        if add == False:
            newseq = newseq[:pos]+len(amas)*"R"+newseq[pos+len(amas):]
    
    keepCluster_sorted = sorted(keepCluster, key=lambda x: x.start)
    return keepCluster_sorted, newseq


def _codeSequence(keepCluster, seq, smooth = False):
    """ for each amas, we put only 1 as a code, example 10001 -> 11111
    
    Parameters
    ----------
    keepCluster : list
        a list [HydroCluster1, ] contains each amas and their position
    seq : string
    
    return seq is the sequence R -> ruptor and 1-> amas with transformation R -> 0 and 1-> 1 in integer
    """
    
    for hyclust in keepCluster:
        start, stop, amas = hyclust.get("all")
        if smooth:
            seq = seq[:start]+len(amas)*"1"+seq[stop:]
        #pass
    #list(map(int, list(seq.replace("R", "0"))))
    return np.array([int(val) for val in seq.replace("R", "0")], dtype=np.uint8)

def _removeProline(seqtransclean):
    """ this function replace Proline by 0 to obtain a perfect binary string, 
    after we redelete 1 or 11 surounded by 0000
    
    Parameter
    ---------
    seqtransclean: string
        a string 1-> hydrophobe P->Proline and 0->other
    
    Return
    ------
    seqpro : string
        a binary code 1->hydrophobe 0-> other
    """
    
    for i in re.finditer(r"(?=(P))", seqtransclean):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()        
        
    seqpro = seqtransclean.replace("P", "0")
    # we remove 1 or 11 surounded by 0000 because P becomes now 0
    # for instance 0P001000011001000 -> becomes 00001000011001000 
    # so we can remove the first 1 like  00000000011001000
    for i in re.finditer(r"(?=(0{4}1{1,2}0{4}))", seqpro):
        match  = i.group(1)
        lmatch = len(match)
        start  = i.start()
        stop   = start + lmatch
        #print (start, stop, seqpro[start:stop])
        
        # replacement
        seqpro = seqpro[:start] + "0"*lmatch + seqpro[stop:]
    
    return seqpro

  
def _getDensity(seqbin, t=0.1):
    """ computes the number of hydrophobe in a window
    
    Parameter
    ---------
    seqbin : string
        a string code of the sequence 0 or 1
    
    Return
    ------
    hydrotable : list
        a table containing for each window and position the number of hydrophobe
    """
    
    # each size of word
    begin = True
    
    # we reduce the windows if sequence are too small
    #if len(seqbin) <17:
    #    lwindows = [len(seqbin)]
    #else:
    #    lwindows = [17]#range(10,51)#[10, 20, 30, 40, 50]

    size_word = len(seqbin) if len(seqbin) < 17 else 17

    #for size_word in lwindows:
    lscore= []
    lthres=[]
    
    # improvement lscore += [-1] * (size_word/2)
    begin = True
    for iter_seq in range(0, len(seqbin)-size_word+1):
        score = sum(seqbin[iter_seq: iter_seq + size_word])/size_word
        
        if begin:
            # for the begining of the window we put the first score
            lscore += [score] * int(size_word/2)
            lthres += [score > t]* int(size_word/2)
            begin = False
        else:
            lscore.append(score)
            lthres.append(score > t)
    
    diff = len(seqbin) - len(lscore)
    # for the end of the position in the window we put the same score
    lscore += [score] * diff 
    lthres += [score > t] * diff 
    return lscore, lthres
    
    
def _sumDensityWindow(hydrotable):
    """sym density for each position in each window
    
    Parameter
    ---------
    hydrotable : list 
        a table containing [density for each position for window1, 
                            density for each position for window2, ...]
                            
    Return
    ------
    list_score_position : list
        for each position of sequence sum of hydrophobe computed by all windows
    """
    
    if len(hydrotable[0]) <17:
        lwindows = [len(hydrotable[0])]
    else:    
        lwindows = [17]#range(10,51)
        
    # SUM HYDROPHOBE WINDOWS
    list_score_position = []
    for position in range(len(hydrotable[0])):
        scorebyposition = 0 # sum the score for the current position
        scorebypositionlinker = 0 # sum the score for the current position
        lscorebyposition = [] # to count the position score of windows match
        for ite, window in enumerate(hydrotable):
            score = window[position]
            
            scorebyposition += window[position]
            lscorebyposition.append(window[position])
        list_score_position.append(scorebyposition/float(len(lscorebyposition)))
    return list_score_position
    

def _findMinima(list_score_position, thresholds, size, t=0.2, nb_under_threshold=4):
    """ alternative to the old _findMinima function
    list_score_position : list
        for each position, sum by windows / size windows of hydrophobe after smooth
    size : int
        length of the amino acid sequence
    """
    #print(thresholds)
    dlimit_domain = dict()
    pres_start, prev_stop = None, None
    start, stop = None, None
    for i in range(size):
        if thresholds[i]:
            if start == None:
                if prev_stop != None:
                    if i - prev_stop >= nb_under_threshold:
                        start = i
                    else:
                        start = prev_start
                else:
                    start = i
                    dlimit_domain[start] = start
            else:
                stop = i
                #print(start, stop)
                dlimit_domain[start] = stop
        else:
            #print(start, stop)
            if start != None or stop != None:
                prev_start = start
                prev_stop = stop
                start, stop = None, None

    #print(dlimit_domain)
    
    limit_domain = []
    for k in dlimit_domain:
        #if dlimit_domain[k]:
        limit_domain.append(k)
        limit_domain.append(dlimit_domain[k])
        
    return sorted(limit_domain)
    
    
def _findMinima_old(list_score_position, thresholds, size, t=0.2):
    """ find the position where the curbe throught the threishold. 
    And after find the local minima
    
    Parameters
    ----------
    list_score_position : list
        for each position, sum by windows / size windows of hydrophobe after smooth
    size : int
        length of the amino acid sequence
    t : float
        the threshold
    
    Return
    ------
    limit_domain : list
        is the position of minima [beg,end, beg, end, beg, end...]
    """
    size_under_threshold = 4
    
    # SEARCH MINIMA TO FIND LIMITS
    limit_domain = []
    state = ""
    for ite, i in enumerate(list_score_position):
        if ite == len(list_score_position)-1:
            continue
        
        first = i
        second = list_score_position[ite+1]
        # 0.0  .1 [.2 .4] .5 we throught the limit
        if first < t and second >= t:
            findminima = False
            state = "UP"
            # find the local minimal to go this limit
            prev = first
            cnt = 0
            start = None
            for j in range(ite-1, -1,-1):
                # check print j, list_score_position[j]
                if prev <= list_score_position[j]:
                    #if not start:
                        #start = j
                    #cnt += 1
                #else:
                    #cnt = 0
                    #start = None
                #if cnt == size_under_threshold:
                #if prev <= list_score_position[j]:
                    #print(">>", ite, start)
                    #print "minima", j
                    limit_domain.append(j)
                    findminima = True
                    break   
                prev = list_score_position[j]
            if not findminima:
                limit_domain.append(0)
        
        # .5 [.4 .2] .1 we down the limit
        if first >= t and second < t:
            state = "DOWN"
            findminima = False
            
            # find the local minimal after the limit
            prev = second
            cnt = 0
            start = None
            below = True
            for j in range(ite+2, len(list_score_position)):
                #print("#", j, list_score_position[j])
                #if list_score_position[j] >= t:
                    #below = False
                #if not below :
                    #break
                
                #if prev <= list_score_position[j]:
                #if list_score_position[j] < t :
                    #if not start:
                        #start = j
                    #cnt += 1
                #else:
                    #cnt = 0
                    #start = None
                #if cnt == size_under_threshold:
                if prev <= list_score_position[j]:
                    #print ("<<", ite, start, cnt, size_under_threshold)
                    findminima = True
                    # we add the first position because we down directly, it means there are hydrophobe in Nter
                    if len(limit_domain) == 0:
                        limit_domain.append(0)   
                    limit_domain.append(j)
                    #print(limit_domain)
                    break
                prev = list_score_position[j]
            
            # These is the end of the sequence no minima, so we put the last limit
            if findminima==False and below:
                limit_domain.append(size-1)
            
    # we add the last position (size of the sequence) because the Cter are hydrophobe
    if state == "UP":
        limit_domain.append(size-1)        
    
    
    if len(limit_domain) == 0:
        limit_domain.append(0)
        limit_domain.append(size-1)
    return limit_domain


def _findAccurateLimit(limit_domain, keepCluster, seq, final_only = False):
    """
    """
    #print(limit_domain)
    #print(list_of_group)
    #list_of_group.sort()
    list_final_limit = []
    #for i in range(0,len(limit_domain), 2):
        #if i >= len(limit_domain)-1:
            #break
        #beg, end = limit_domain[i], limit_domain[i+1]
        #potential_domain = []
        #for jb, je in list_of_group:
            #gap = abs(beg-jb) + abs(end-je)
            ##print(beg, end, jb, je)
            #potential_domain.append([gap, (jb, je)])
        ##print(">>", potential_domain, min(potential_domain))
        #if len(potential_domain) == 0:
            #if not final_only:
                #list_final_limit.append((beg, end))
        #else:
            #size, lim = min(potential_domain)
            #list_final_limit.append(lim)
            #tmp_list_of_group = list()
            #for jb, je in list_of_group:
                #if jb != lim[0] and je != lim[1]:
                    #tmp_list_of_group.append([jb, je])
            #list_of_group = tmp_list_of_group[:]
            
    for i in range(0,len(limit_domain), 2):
        if i >= len(limit_domain)-1:
            break
        beg, end = limit_domain[i], limit_domain[i+1]
        potential_domain = []
        for clust in keepCluster:
            start = max(clust.start, beg)
            stop = min(clust.stop, end)
            if stop-start > 0:
                potential_domain.append((clust.start, clust.stop))
        if len(potential_domain) == 0:
            if not final_only:
                list_final_limit.append((beg, end))
        else:
            potential_domain.sort()
            list_final_limit.append([potential_domain[0][0], potential_domain[-1][1]])
    
    final_domains = []
    for start, stop in list_final_limit:
        #print(start, stop)
        hydro_dom = []
        for hydroclust in keepCluster:
            hstart, hstop = hydroclust.get("start"), hydroclust.get("stop")
            #print(start, stop, hstart, hstop)
            if start <= hstart and hstop <= stop:
                hydro_dom.append(hydroclust)
        final_domains.append(DomHCA(start, stop, hydro_dom))
    return final_domains
    
        
def _clusterizeCluster(seqbin):
    """ clusterize the hydrophobic clusters based on their distances
    
    Parameters
    ----------
    seqbin : string
        positions of hydrophobic cluster 0011111100001111110
    
    Return
    ------
    selected: list
        the selected list of hydrophobic clusters
    """
    # get cluster and central position
    dcluster = {}
    cluster = False
    for ite, is_amas in enumerate(seqbin):
        if not cluster and is_amas == 1:
            start = ite
            cluster = True
        elif cluster and is_amas == 0:
            stop = ite
            ocluster = (start, stop)
            dcluster[start] = ocluster 
            cluster = False
    else:
        # is the last position a cluster?
        if cluster:
            stop = ite+1
            ocluster = (start, stop)
            dcluster[start] = ocluster 
            
    # no or one cluster so no need to cluster
    if len(dcluster) <=1:
        l = []
        for ite, ocluster in dcluster.items():
            l.append(ocluster)
        return l
    
    selected = set()
    nbgroup_previous = len(dcluster)
    nbgroup, loop = 0,  0
    dcluster_indexes = sorted(list(dcluster.keys()))
    while nbgroup_previous != nbgroup:
        # update nbgroup_previous
        nbgroup_previous = nbgroup #len(dcluster)
       
        if nbgroup == 1:
            break
        # compute distance between cluster
        name_select, min_new_name, dmin = [[],[]], None, 1e6 
        indexes = []
        for i in range(len(dcluster_indexes)-1):
            ite_i = dcluster_indexes[i]
            oclusteri = dcluster[ite_i]
            ite_j = dcluster_indexes[i+1]
            oclusterj = dcluster[ite_j]
        
            # search for the way
            if oclusterj[0] > oclusteri[1]:
                dij = oclusterj[0] - oclusteri[1]
                flag = 0
            else:
                dij = oclusteri[0] - oclusterj[1]
                flag = 1
            if dmin > dij:
                dmin = dij
                indexes = [ite_i, ite_j] if flag == 0 else [ite_j, ite_i]
            
        if indexes:
            new_clust = (dcluster[indexes[0]][0], 
                         dcluster[indexes[1]][1])
            # graph
            selectedi = dcluster[indexes[0]]
            selectedj = dcluster[indexes[1]]
            selected.add(selectedi)
            selected.add(selectedj)
            selected.add(new_clust)
            #print(selectedi, selectedj, new_clust)
            #new_index = min(indexes)
            dcluster_indexes.remove(max(indexes))
            dcluster[min(indexes)] = new_clust
            # delete old composante
            del dcluster[max(indexes)]
           
            # update nbgroup
            nbgroup = nbgroup_previous - 1
    
    return list(selected)

def compute_loc_score(seq, clusters, dist=16):
    """ compute local score
    """
    models ={"C": 0.38888888888888884, "P": 0.8333333333333333, "S": 0.8518518518518517, 
             "M": 0.27777777777777773, "Q": 0.8518518518518517, "V": 0.27777777777777773, 
             "F": 0.0185185185185186, "A": 0.6851851851851851, "Y": 0.4259259259259259, 
             "D": 0.9074074074074073, "I": 0.11111111111111116, "R": 0.8148148148148148, 
             "E": 0.9074074074074073, "W": 0.2962962962962963, "K": 1.0, "H": 0.6851851851851851, 
             "G": 0.7962962962962962, "L": 0.0, "T": 0.7592592592592592, "N": 0.8703703703703702}
    #models = {"K": 0.17, "A": 0.0, "C": -0.16, "S": 0.09, "F": -0.36,
              #"G": 0.06, "M": -0.22, "H": 0.0, "E": 0.12, "I": -0.31,
              #"R": 0.07, "Y": -0.14, "L": -0.37, "P": 0.08, "V": -0.22,
              #"D": 0.12, "T": 0.04, "Q": 0.09, "W": -0.21, "N": 0.1}
    scores = list()
    pos_lng_clusters = np.zeros(len(seq), dtype=np.uint8)
    pos_clusters = np.zeros(len(seq), dtype=np.uint8)
    part_clusters = np.zeros(len(seq), dtype=np.uint8)
    part_lng_clusters = np.zeros(len(seq), dtype=np.uint8)
    for clust in clusters:
        hclust = np.array([int(val) for val in clust.hydro_cluster], dtype=np.uint8)
        part_clusters[clust.start: clust.stop] = 1
        pos_clusters[clust.start: clust.stop] = hclust
        if len(hclust) > 2:
            part_lng_clusters[clust.start: clust.stop] = 1
            pos_lng_clusters[clust.start: clust.stop] = hclust
            
    for i in range(len(seq)):
        subseq = seq[max(0,i-dist): min(i+dist+1, len(seq))]
        hscore = 0
        for aa in models:
            cnt = subseq.count(aa)
            hscore += cnt * models[aa] #* -1
        hscore /= len(subseq)
        # -range(i, i+1+d), i, range(1, i+1+d)
        subh_exp = part_clusters[max(0,i-dist): min(i+dist+1, len(seq))]
        cnth_exp = sum(subh_exp)
        exp = cnth_exp - 1
        subh_obs = part_lng_clusters[max(0,i-dist): min(i+dist+1, len(seq))]
        obs = 0
        for j in range(0, len(subh_obs)-1):
            if subh_obs[j] == 1:
                if j+1 < len(subh_obs) and subh_obs[j+1] == 1:
                    obs += 1
                #elif j+2 < len(subh_obs) and subh_obs[j+2] == 1:
                    #obs += 1
                #elif j+3 < len(subh_obs) and subh_obs[j+3] == 1:
                    #obs += 1
                #elif j+4 < len(subh_obs) and subh_obs[j+4] == 1:
                    #obs += 1
        if cnth_exp >= 2:
            score2 = obs / exp
        else:
            score2 = 0
        #print(i, subh_exp, cnth_exp, subh_obs, obs, score2,  cnth_exp/len(subseq))
        scores.append(score2)
      
    #smooth_scores = list()
    ## smooth averaging
    #for i in range(len(seq)):
        #sub_scores = scores[max(0,i-dist): min(i+dist+1, len(seq))]
        #smooth_scores.append(sum(sub_scores)/len(sub_scores))
    return np.array(scores)

def _annotation_aminoacids(seq, t=0.1, method="domain", return_seqbin=False, verbose=False):
    """ The amino acids annotation function. Two methods are avaliable: 'domain'
    and 'cluster' 
    
    Parameters
    ----------
    seq : string
        the amino acid protein sequence
    t : float
        parameter controlling the domain creation based on cluster density
    method : string
        the method used, 
        domain: will return a list of domain positions
        cluster: will return a list of cluster positions
    verbose: bool
        print interesting stuff
    
    Return:
    -------
    annotat : list
        the annotation results
    """
    annotat = {"cluster": [], "domain": [], "scores": []}
    # obsolete identify low complexity segment by SEG algorithm
    low_complexity = []#getSEGLW(fseq)
    seqori = seq[:]
    # transforme sequence
    if verbose:
        print("Transform amino acid sequence into binary sequence")
    seqtrans2 = _transformSequence(seqori, low_complexity)
    #print(seqori)
    #print("".join(str(val) for val in seqtrans2))
    
    # get amas
    # list_amas -> [ [position of the amas, amas], [...], ... ]
    if verbose:
        print("Create hydrophobic clusters")
        
    #print(seqtrans2)
    if method == "cluster_removeP":
        seqtrans2 = seqtrans2.replace("P", "0")#_removeProline(seqtrans2)
    
    #print(seqtrans2)
    list_amas, seqamas = _getAmas(seqtrans2)
    annotat["cluster"] = list_amas[:]
    #print(seqamas)
    
    keepCluster, newseq = _removeSmallAmas(list_amas, seqamas)
    seqbin = _codeSequence(keepCluster, newseq, smooth = True)
    #print(seqbin)
    #print(keepCluster)
    
    if method == "domain":
        # get density of hydrophobe mean by windows
        hydrotable, thresolds = _getDensity(seqbin, t)
        #print(hydrotable)
        # mean by position of all windows
        #list_score_position = _sumDensityWindow(hydrotable) # TODO Useless, the size of the window is fixe at 17
        limit_domain = _findMinima(hydrotable, thresolds, len(seqbin), t)
        #print(limit_domain)
        if verbose:
            print("Group HC to delineate domains")
        #list_of_group = _clusterizeCluster(seqbin)
        # get limits
        #print(list_of_group)
        domains = _findAccurateLimit(limit_domain, keepCluster, seqamas, final_only=True)
        #compute_pvalues_domains(domains)
        annotat["domain"] = domains
    
    if return_seqbin:
        return annotat, seqbin
    else:
        return annotat
    


def _annotation(output, inputf, seq_type="aminoacid", t=0.1, method="domain", verbose=False):
    """ The main annotation function. Two methods are avaliable: 'domain' and 
    'cluster' 
    
    Parameters
    ----------
    inputf: string
        path of the input file
    seq_type: string, ["aminoacids", "nucleotides"]
        the type of biological sequence
    t : float
        parameter controlling the domain creation based on cluster density
    method : string
        the method used, 
        domain: will return a list of domain positions
        cluster: will return a list of cluster positions
    verbose: bool
        print interesting stuff
    
    Return:
    -------
    danno : dictionarry
        the annotation for each protein or each frame of each nucleotide
        sequences
    """
    with open(output, "w") as outf:
        outf.write("""# pyHCA v0.1 segmentation results
# 
# Format:
# 
# >'protein_id' 'protein_length' 'hca_score computed on the whole sequence'
# domain 'domain_start' 'domain_stop' 'hca_score' 'hca_pvalue' (if -m domain is used)
# cluster 'cluster_start' 'cluster_stop' 'cluster_pattern'
# 
# The hca_score and associated p-value provide a way to measure the foldability
# of a protein, i.e how similar is the score compared to scores from disordered
# sequences.
# Low p-values correspond to scores at the tail of the distribution of scores 
# for disordered protein sequences.
# 
# /!\ Warning /!\
# 1- The score computed at the whole protein level (in the line with '>') is for 
# information only as some people found it useful.
# No p-value is associated to this score as the scores used in the distributions
# don't come from full protein sequences but domain or "disordered regions" of 
# comparable lengths.
#
# 2- similarly, scores are displayed even for HCA domain shorted than 30 amino 
# acids. 
# As the sequences of length lower than 30 amino acids where filtered out to
# compute distributions of scores, no p-values are given.
# 
# In these two cases, the scores provided must be analyzed carefully, keeping
# in mind their origin and initial purpose
# /!\ Warning /!\
#
#

""")
        if seq_type == "aminoacid":
            for prot, sequence in read_multifasta_it(inputf, verbose):
                #for prot in dseq:
                #sequence = str(dseq[prot].seq)
                annotations = _annotation_aminoacids(sequence, t=t, method=method, verbose=verbose)
                score, pvalue = compute_disstat(0, len(sequence), annotations["cluster"])
                outf.write(">{} {} {:.3f} {:.3f}\n".format(prot, len(sequence), pvalue, score))
                for domannot in annotations["domain"]:
                    outf.write("{}\n".format(str(domannot)))
                for clustannot in annotations["cluster"]:
                    outf.write("{}\n".format(str(clustannot)))
        else:
            cnt, nb_dot = 0, 0
            #for name in dseq:
            for name, sequence in read_multifasta_it(path, verbose):
                #for strand, frame, start, protseq in six_frames(dseq[name]):
                for strand, frame, start, protseq in six_frames(sequence):
                    cnt += 1
                    if cnt == 1000:
                        cnt = 0
                        sys.stdout.write(".")
                        sys.stdout.flush()
                        nb_dot += 1
                    if nb_dot == 80:
                        nb_dot = 0
                        sys.stdout.write("\n")
                    if strand > 0:
                        new_name = "{}_5'3'_Frame_{}_start_{}".format(name, frame+1, start+1)
                    else:
                        new_name = "{}_3'5'_Frame_{}_start_{}".format(name, frame+1, start+1)
                    
                    annotations = {"cluster": [], "domain": []}
                    cur_annotation = _annotation_aminoacids(protseq, t=t, method=method, verbose=verbose)
                    for domannot in cur_annotation["domain"]:
                        annotations["domain"].append(domannot)
                    for clustannot in cur_annotation["cluster"]:
                        annotations["cluster"].append(clustannot)
                        
                    score, pvalue = compute_disstat(0, len(protseq), annotations["cluster"])
                    if annotations:
                        outf.write(">{} {} {:.3f}\n".format(new_name, len(protseq), score))
                        for domannot in annotations["domain"]:
                            outf.write("{}\n".format(str(domannot)))
                        for clustannot in annotations["cluster"]:
                            outf.write("{}\n".format(str(clustannot)))
            sys.stdout.write("\n")



def _scores(output, dseq, seq_type="aminoacid", t=0.1, method="domain", verbose=False, dist=16):
    """ The main annotation function. Two methods are avaliable: 'domain' and 
    'cluster' 
    
    Parameters
    ----------
    dseq : dictionary
        the biological sequences, keys are string, values are biopython Sequence
        object from SeqIO
    seq_type: string, ["aminoacids", "nucleotides"]
        the type of biological sequence
    t : float
        parameter controlling the domain creation based on cluster density
    method : string
        the method used, 
        domain: will return a list of domain positions
        cluster: will return a list of cluster positions
    verbose: bool
        print interesting stuff
    
    Return:
    -------
    danno : dictionarry
        the annotation for each protein or each frame of each nucleotide
        sequences
    """
    with open(output, "w") as outf:
        if seq_type == "aminoacid":
            #for prot in dseq:
            for prot, sequence in read_multifasta_it(inputfile):
                #sequence = str(dseq[prot].seq)
                annotations = _annotation_aminoacids(sequence, t=t, method=method, 
                                                    verbose=verbose, dist=dist)
                outf.write(">{} {}\n".format(prot, len(sequence)))
                if method =="domain":
                    posdomains = np.zeros(len(sequence), dtype=np.uint8)
                    for domannot in annotations["domain"]:
                        posdomains[domannot.start: domannot.stop] = 1
                    for i in range(len(sequence)):
                        outf.write("{:.5f}\t{}\n".format(annotations["scores"][i], posdomains[i]))
                else:
                    for i in range(len(sequence)):
                        outf.write("{:.5f}\tNaN\n".format(annotations["scores"][i]))
        else:
            cnt, nb_dot = 0, 0
            #for name in dseq:
            for name, sequence in read_multifasta_it(inputfile):
                for strand, frame, start, protseq in six_frames(sequence):
                    cnt += 1
                    if cnt == 1000:
                        cnt = 0
                        sys.stdout.write(".")
                        sys.stdout.flush()
                        nb_dot += 1
                    if nb_dot == 80:
                        nb_dot = 0
                        sys.stdout.write("\n")
                    if strand > 0:
                        new_name = "{}_5'3'_Frame_{}_start_{}".format(name, frame+1, start+1)
                    else:
                        new_name = "{}_3'5'_Frame_{}_start_{}".format(name, frame+1, start+1)
                    
                    annotations = _annotation_aminoacids(protseq, t=t, method=method, verbose=verbose, dist=dist)                            
                    if annotations:
                        outf.write(">{} {}\n".format(new_name, len(protseq)))
                        if method =="domain":
                            posdomains = np.zeros(len(protseq), dtype=np.uint8)
                            for domannot in annotations["domain"]:
                                posdomains[domannot.start: domannot.stop] = 1
                            for i in range(len(protseq)):
                                outf.write("{:.5f}\t{}\n".format(annotations["scores"][i], pos_domains[i]))
                        else:
                            for i in range(len(protseq)):
                                outf.write("{:.5f}\tNaN\n".format(annotations["scores"][i]))
            sys.stdout.write("\n")

def _process_params():
    """ Process parameters when the script annotateHCA is directly called
    """
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "annotate"))
    parser.add_argument("-i", action="store", dest="inputf", required=True,
        help="an amino-acid sequence files in fasta format")
    parser.add_argument("-o", action="store", dest="outputf", required=True,
        help="the output file with annotation")
    parser.add_argument("-v", action="store_true", dest="verbose", default=False,
        help="keep temporary results")
    parser.add_argument("-m", action="store", dest="method", default="cluster1", 
        choices=["cluster","domain"], help=("method to use, cluster: will "
        "report *the hydrophobic clusters found in the sequence, domain: will "
        "delineate domains based on the hydrophobic cluster profile of the "
        "sequence"))
    parser.add_argument("-t", action="store", dest="seqtype", 
        default="aminoacid", choices=["aminoacid","nucleotide"], help=("the "
        "type of the biological sequences passed in the input file"))
    params = parser.parse_args()
    return params

    
def main_segment():
    """ the main function is called after direct invocation of the software
    """
    params = _process_params()
    
    # annotation
    _annotation(params.outputf, params.inputf, seq_type=params.seqtype, t=0.1, 
                              method=params.method, verbose=params.verbose)
    sys.exit(0)

    
if __name__ == "__main__":
    main_segment()
    
