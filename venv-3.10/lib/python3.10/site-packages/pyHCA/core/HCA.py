#!/usr/bin/env python
""" The classHCA is an ensemble of class to describe hydrophobic clusters and
domain created through hydrophobic clusters.
"""

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "tristan.bitard-feildel [you know what] upmc.fr"
__institute__ = "UPMC"

import time, os, tempfile, sys
import Bio
from Bio import Seq 
from Bio import SeqIO
try:
    from Bio.Alphabet import IUPAC
except:
    from Bio.Data import IUPACData as IUPAC

import numpy as np
import warnings
import functools

from .classHCA import Seq as SeqHCA
from .annotateHCA import _annotation_aminoacids
from .drawHCA import make_svg, getSVGheader, colorize_positions
from .seq_util import compute_conserved_positions
from .tremoloHCA import search_domains, group_resda, write_tremolo_results
from .external import interpro_search

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used."""

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning) #turn off filter 
        warnings.warn("Call to deprecated function {}.".format(func.__name__), category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning) #reset filter
        return func(*args, **kwargs)

    return new_func

class HCA(object):
    """ HCA class provides an API interface to all the standalone programs
    """
    def __init__(self, querynames=None, seq=None, seqfile=None, file_format="fasta"):
        """ create an HCA instance, accept only one sequence at a time
        
        Parameters:
        -----------
        seq: list of instances or instance
            instance can either be a string, a Bio.Seq object or a Bio.SeqRecord object
        querynames: list or string
            if seq is instance of string or Bio.Seq, provide a name to the protein sequences, if not used "query_<num_idx>" is used
        seqfile: string
            path to the file containing protein sequences
        file_format: string
            BioPython supported file format for seqfile
            
            
        Usage:
        ------
        
        >>> # instanciation of a single sequence
        >>> seq_str = "ATGYHVVLIVQEAGFHILLV"
        >>> hca = HCA(seq=seq_str, querynames="my_query")
        >>>        
        >>> from Bio import Seq
        >>> seq_bio = Seq.Seq("ATGYHVVLIVQEAGFHILLV")
        >>> hca = HCA(seq=seq_bio) # without specifying querynames, automatically set it up to query_1
        >>>
        >>> from Bio import SeqRecord
        >>> seq_rec = SeqRecord.SeqRecord(id="protein_1", seq="ATGYHVVLIVQEAGFHILLV")
        >>> hca = HCA(seq=seq_rec) # SeqRecord's id attribute is used 
        >>>

        >>> # instanciation of a list of sequences
        >>> seq_str_lst = ["ATGYHVVLIVQEAGFHILLV", "AGVVLATGYHHILLVFHILLV"]
        >>> hca = HCA(seq=seq_str_lst, querynames=["my_query_1", "my_query_2"])
        >>>
        """
        
        if isinstance(seq, str) and seqfile is None and os.path.isfile(seq):
            raise ValueError("Error, please use seqfile= argument to pass a file as an argument")
        
        if seq is None and seqfile is None:
            raise ValueError("Error, you need to provide a sequence (or a list of sequences), or a valid path to a file to instantiate a HCA object")
        
        if seq is not None and seqfile is not None:
            raise ValueError("Error, you need to provide either a sequence (or a list of sequences), or a valid path to a file to instantiate a HCA object")
        
        self.sequences = list()
        self.msa_seq = list()
        self.is_msa = False
        self.querynames = list()
        self.__seqbin = dict()
        self.__domains = dict()
        self.__clusters = dict()
        self.__scores = dict()
        self._number_of_sequences = 0
        
        # attributes to keep in memory if computation were previously done
        self._segments_done = False
        self._segments_done_with_t = -1
        self._tremolo_done = False
        
        if seqfile is None and seq is not None:
            # list of sequences or single sequence
            qnames = list()
            if isinstance(seq, list):
                for s in seq:
                    new_seq, qname = check_seq_type(s)
                    self.sequences.append(new_seq)
                    qnames.append(qname)
                    self._number_of_sequences += 1
            else:
                seq, qname = check_seq_type(seq)
                self.sequences.append(seq)
                qnames = [qname]
                self._number_of_sequences += 1
            
            # adapt queryname to list type
            if querynames is None:
                all_qnames = list(set(qnames))
                if len(all_qnames) == 1 and all_qnames[0] == "query":
                    for idx in range(len(self.sequences)):
                        name = "query_{}".format(idx+1)
                        self.querynames.append(name)
                        # initialize containers for domains and clusters
                        if name not in self.__domains:
                            self.__seqbin[name] = ""
                            self.__domains[name] = list()
                            self.__clusters[name] = list()
                        else:
                            raise ValueError("Error, multiple proteins found with the same name {}".format(name))
                            
                else:
                    for name in qnames:
                        self.querynames.append(name)
                        # initialize containers for domains and clusters
                        if name not in self.__domains:
                            self.__seqbin[name] = ""
                            self.__domains[name] = list()
                            self.__clusters[name] = list()
                        else:
                            raise ValueError("Error, multiple proteins found with the same name {}".format(name))
            elif isinstance(querynames, list):
                if len(querynames) != len(self.sequences):
                    raise ValueError("Error, queryname argument must be a list of same size as seq argument (one name per seq")
                else:
                    for name in querynames:
                        self.querynames.append(name)
                        # initialize containers for domains and clusters
                        if name not in self.__domains:
                            self.__seqbin[name] = ""
                            self.__domains[name] = list()
                            self.__clusters[name] = list()
                        else:
                            raise ValueError("Error, multiple proteins found with the same name {}".format(name))
            elif isinstance(querynames, str):
                self.querynames.append(querynames)
            else:
                raise ValueError("Error, querynames should be a list of names or a string")
        else:
            if not os.path.isfile(seqfile):
                raise ValueError("Error, argument path ({}) of seqfile must be a file".format(seqfile))
            with open(seqfile) as inf:
                for record in SeqIO.parse(seqfile, file_format):
                    name = record.id
                    self.querynames.append(name)
                    self.sequences.append(str(record.seq))
                    self.__seqbin[name] = ""
                    self.__domains[name] = list()
                    self.__clusters[name] = list()
                    self._number_of_sequences += 1
                    
        # no sequences?
        if self._number_of_sequences == 0:
            raise ValueError("No sequences found")
        
        # check if sequence is a MSA sequence
        is_msa, msa_seq, sequences = check_if_msa(self.querynames, self.sequences)
        self.msa_seq = msa_seq[:]
        self.is_msa = is_msa
        if is_msa:
            self.sequences = sequences[:]
        
    ### SEG-HCA
    @property
    def domains(self):
        return self.__domains
    
    @property
    def clusters(self):
        return self.__clusters
    
    @property
    def seqbin(self):
        return self.__seqbin
    
    @property
    def clusters(self):
        return self.__clusters
    
    @property
    def scores(self):
        return self.__scores
    
    @domains.setter
    def domains(self, domains):
        self.__domains = domains
    
    @clusters.setter
    def clusters(self, clusters):
        self.__clusters = clusters
    
    @scores.setter
    def scores(self, scores):
        self.__scores = scores
        
    def segments(self, t=0.1, verbose=False):
        """ run the segmentation of protein sequences into HCA domains, store domain and cluster positions
        
        Parameters
        ----------
        t: float
            cutoff value used to segment the protein, trained on SCOP protein sequences
        verbose: bool
            print information
        
        Usage:
        ------
        
        >>> seq_str_lst = ["ATGYHVVLIVQEAGFHILLV", "AGVVLATGYHHILLVFHILLV"]
        >>> hca = HCA(seq=seq_str_lst, querynames=["my_query_1", "my_query_2"])
        >>> hca.segments()    
        >>>
        """
        t1 = time.time()
        if not self._segments_done or t != self._segments_done_with_t:
            if t != self._segments_done_with_t and verbose:
                print("Running segmentation with a different t value ({} -> {})".format(self._segments_done_with_t, t))
            for i in range(len(self.sequences)):
                seq = self.sequences[i]
                prot = self.querynames[i]
                annotations, seqbin = _annotation_aminoacids(seq, t=t, method="domain", return_seqbin=True, verbose=verbose)
                self.seqbin[prot] = seqbin
                self.domains[prot] = annotations["domain"]
                self.clusters[prot] = annotations["cluster"]
                self.scores[prot] = annotations["scores"]
                self._segments_done = True
                self._segments_done_with_t = t
        elif verbose:
            print("Warning, segmentations already performed return precomputed results, "
                    "used a different t value ({}) to get different results".format(self._segments_done_with_t), file=sys.stderr)
        t2 = time.time()
        if verbose:
            print("Segmentation done in {}".format(t2-t1))
    
    def get_seqbin(self, prot=None):
        """ function wrapper to return HCA binary sequence.
        If only one sequence was provided return a string
        If multiple sequences were provided return a dictionary with
        protein querynames as keys and the list of string as values.
        
        Parameters:
        -----------
        prot: None or string
            If None get all domains results of every proteins in a dictionary.
            If a string is provided and corresponds to an analysed proteins, 
            only return the list of domain of the protein.
            If only one sequence was analysed a list of domain is returned
        
        Usage:
        ------
        
        >>> seq_str_lst = ["ATGYHVVLIVQEAGFHILLV", "AGVVLATGYHHILLVFHILLV"]
        >>> hca = HCA(seq=seq_str_lst, querynames=["my_query_1", "my_query_2"])
        >>>        
        >>> domains = hca.get_seqbin("my_query_1")
        >>> print(domains)
        [<pyHCA.core.classHCA.DomHCA object at 0x7f604cab0b88>]
        >>> print(domains[0])
        domain  4       20      0.00025890542325190946  0.7058823529411765
        >>>
        
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.domains:
                raise KeyError("Error, unable to find proteins '{}' in domain results".format(prot))
            return self.seqbin[prot]
        elif self._number_of_sequences == 1:
            return self.seqbin[self.querynames[0]]
        else:
            return dict([(prot, self.seqbin[prot]) for prot in self.querynames])
    
    
    def get_domains(self, prot=None):
        """ function wrapper to return HCA domain annotation.
        If only one sequence was provided return a list of domains.
        If multiple sequences were provided return a dictionary with
        protein quernames as keys and the list of domains as values.
        
        Parameters:
        -----------
        prot: None or string
            If None get all domains results of every proteins in a dictionary.
            If a string is provided and corresponds to an analysed proteins, 
            only return the list of domain of the protein.
            If only one sequence was analysed a list of domain is returned
        
        Usage:
        ------
        
        >>> seq_str_lst = ["ATGYHVVLIVQEAGFHILLV", "AGVVLATGYHHILLVFHILLV"]
        >>> hca = HCA(seq=seq_str_lst, querynames=["my_query_1", "my_query_2"])
        >>> hca.segments()    
        >>> domains = hca.get_domains("my_query_1")
        >>> print(domains)
        [<pyHCA.core.classHCA.DomHCA object at 0x7f604cab0b88>]
        >>> print(domains[0])
        domain  4       20      0.00025890542325190946  0.7058823529411765
        >>>
        
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.domains:
                raise KeyError("Error, unable to find proteins '{}' in domain results".format(prot))
            return self.domains[prot]
        elif self._number_of_sequences == 1:
            return self.domains[self.querynames[0]]
        else:
            return dict([(prot, self.domains[prot]) for prot in self.querynames])
    
    def get_clusters(self, prot=None):
        """ function wrapper to return HCA cluster positions. 
        If only one sequence was provided return a list of clusters.
        If multiple sequences were provided return a dictionary with
        as protein querynames as keys and the list of clusters as values.
        
        Usage:
        ------
        
        >>> seq_str_lst = ["ATGYHVVLIVQEAGFHILLV", "AGVVLATGYHHILLVFHILLV"]
        >>> hca = HCA(seq=seq_str_lst, querynames=["my_query_1", "my_query_2"])
        >>> hca.segments()    
        >>> clusters = hca.get_clusters("my_query_1")
        >>> print(clusters)
        [<pyHCA.core.classHCA.HydroCluster at 0x7f604cab18d0>,
         <pyHCA.core.classHCA.HydroCluster at 0x7f604cab1940>]
        >>> print(clusters[0])
        cluster 4       10      1011111
        >>>
        
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.clusters:
                raise KeyError("Error, unable to find proteins '{}' in cluster results".format(prot))
            return self.clusters[prot]
        elif self._number_of_sequences == 1:
            return self.clusters[self.querynames[0]]
        else:
            return dict([(prot, self.clusters[prot]) for prot in self.querynames])
        
    @deprecated
    def get_scores(self, prot=None):
        """ function wrapper to return HCA scores of each domains. 
        If only one sequence was provided return a list of scores.
        If multiple sequences were provided return a dictionary with
        as protein querynames as keys and the list of scores as values.
        
        Usage:
        ------
        
        >>> seq_str_lst = ["ATGYHVVLIVQEAGFHILLV", "AGVVLATGYHHILLVFHILLV"]
        >>> hca = HCA(seq=seq_str_lst, querynames=["my_query_1", "my_query_2"])
        >>> hca.segments()    
        >>> scores = hca.get_scores("my_query_1")
        >>> print(scores[0])
        []
        >>>
        
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.scores:
                raise KeyError("Error, unable to find proteins '{}' in scores results".format(prot))
            return self.scores[prot]
        elif self._number_of_sequences == 1:
            return self.scores[self.querynames[0]]
        else:
            return dict([(prot, self.scores[prot]) for prot in self.querynames])
    
    def save_annotation(self, output):
        """ save the seg-HCA annotation to a file
        
        Parameter
        ---------
        output: string
            the path to the output file
            
        Usage
        -----
        
        >>> hca.save_tremolo("tremolo_results.txt")
        >>>
        
        """
        if not self._segments_done:
            raise RuntimeError("Error, no annotation to save, you must perform an HCA segmentation with the segments() method first")
        with open(output, "w") as outf:
            for i, prot in enumerate(self.querynames):
                outf.write(">{} {}\n".format(prot, len(self.sequences[i])))
                for domannot in self.domains[prot]:
                    outf.write("{}\n".format(str(domannot)))
                for clustannot in self.clusters[prot]:
                    outf.write("{}\n".format(str(clustannot)))
            
    ### DrAW-HCA
    def draw_hca2svg(self, external_annotation=dict(), show_hca_dom=False, outputfile=None):
        """ draw a HCA plot in svg of each sequence
        
        Parameters:
        -----------
        external_annotation: dict
            currently only support external domain annotation. 
            Format: 
            external_annotation[protein_A] = [(start_A1, stop_A1, ext_domain_A1), (start_A2, stop_A2, ext_domain_A2), ...]
            external_annotation[protein_B] = [(start_B1, stop_B1, ext_domain_B1), ...]
            ...
        show_hca_dom: bool
            if set to True display the computed hca domains on the hca plot
        outputfile: None or string
            if a path is provided save the plot in an svg format
        
        Usage:
        ------
        
        >>> svg = hca.draw_hca2svg()
        >>>
        
        """
        if self.is_msa:
            msa_conserved = compute_conserved_positions(dict(zip(self.querynames, self.sequences)), dict(zip(self.querynames, self.msa_seq)))
        
        self.all_svg = dict()
        max_aa = 0
        cnt = 0
        # create hca plot for each sequence
        for i in range(len(self.querynames)):
            prot = self.querynames[i]
            prot_seq = self.sequences[i]
            # read domain annotation if provided
            annotation = {"domains": list()}
            if prot in external_annotation and "domains" in external_annotation[prot]:
                for start, stop, dom in external_annotation[prot]:
                    annotation["domains"].append((start, stop, dom, "!", None))
            if show_hca_dom:
                for dom in self.domains[prot]:
                    start = dom.start
                    stop = dom.stop
                    annotation["domains"].append((start, stop, "domain", "!", None))
            # make svg
            
            if self.is_msa:
                annotation["positions"] = colorize_positions(self.msa_seq[i], prot_seq, msa_conserved[prot], method="rainbow")
            cur_svg, nbaa = make_svg(prot, prot_seq, annotation, cnt)
            self.all_svg[prot] = cur_svg
            if nbaa > max_aa:
                max_aa = nbaa
            cnt += 1
        # write in outputfile if provided
        if outputfile != None:
            svgheader = getSVGheader(max_aa, (cnt+1)*230)
            with open(outputfile, "w") as fout:
                fout.write(svgheader)
                for i in range(len(self.querynames)):
                    prot = self.querynames[i]
                    fout.write(self.all_svg[prot])
                fout.write("</svg>")
        # return the svg dictionary
        return self.all_svg

    def draw_hca2plot(self, external_annotation=dict(), ax=None, show_hca_dom=False, outputfile=None, show=False, **kwargs):
        """ draw a HCA plot in svg of each sequence
        
        Parameters:
        -----------
        external_annotation: dict
            currently only support external domain annotation. 
            Format: 
            external_annotation[protein_A] = [(start_A1, stop_A1, ext_domain_A1), (start_A2, stop_A2, ext_domain_A2), ...]
            external_annotation[protein_B] = [(start_B1, stop_B1, ext_domain_B1), ...]
            ...
        show_hca_dom: bool
            if set to True display the computed hca domains on the hca plot
        outputfile: None or string
            if a path is provided save the plot in an svg format
        
        Usage:
        ------
        
        >>> hca.draw_hca2plot(show=False, outputfile="hca.pdf", dpi=600)
        >>> 
        >>> fig, ax = plt.subplots()
        >>> hca.draw_hca2plot(ax=ax)
        >>> plt.show()
        
        """
        print("Not yet implemented")
    
    ## Tremolo-HCA 
    
    def tremolo(self, prot, start, stop, hhblitsdb="", hhblitsparams="", annotation_path="", annotation="Interpro", evalue=1e-3):
        """ run tremolo hca on a specific part of a protein sequence
        
        Parameters:
        -----------
        prot: string
            protein name, must be a part of input protein names
        start: int 
            domain start position to analysis, inclusive, positive
        stop: int 
            domain stop position to analysis, inclusive, lesser than sequence length
        hhblitsdb: string
            path to the HHblits database to use 
        hhblitsparams: string
            specific options to run HHblits with
        evalue: float
            evalue threshold to filter hhblits hits
        annotation: string
            hhblits target annotation to use: CDD (online annotation) or Interpro (local annotation)
            CDD used the NCBI webserver to annotate hhblits targets
            Interpro required a local file of protein annotations
            
        Usage:
        ------
        
        >>> tremolo_res = hca.tremolo("my_query_1", 1, 10, hhblitsdb="HHSuite/uniprot20_2016_02", annotation_path="protein2ipr.dat")
        >>> print(tremolo_res)
        ...
        
        """
        # check query sequence and domain boundaries
        if prot not in self.querynames:
            raise ValueError("Protein {} not found in list of possible protein queries [{}]".format(prot, ", ".join(self.querynames)))
        prot_idx =self.querynames.index(prot)
        sequence = self.sequences[prot_idx]
        start -= 1 # shift to 0 start
        if stop <= len(sequence) and start > -1:
            self.__tremolo_domains = [(start, stop)]
        else:
            raise ValueError("Invalid domain boundaries for protein {}, ({}, {}), sequence length is {}".format(prot, start+1, stop, len(sequence)))
        self.__tremolo_query = SeqHCA([prot], [""], [sequence], len(sequence))

        # prepare workdir
        self.__tremolo_workdir = tempfile.mkdtemp()
        
        # run hhblits
        self.__tremolo_targets, alltargetids = search_domains(self.__tremolo_query, self.__tremolo_domains, hhblitsdb, evalue, hhblitsparams, self.__tremolo_workdir )
        if alltargetids == []:
            msg = "Unable to find any targets with hhblits in database {}".format(hhblitsdb)
            msg+= "with parameters {}\n".format(hhblitsparams)
            msg+= "Please try less stringent parameters or a different database"
            raise RuntimeError(msg)
        
        
        # get domain from Interpro
        self.__tremolo_annotation = interpro_search(alltargetids, self.__tremolo_workdir , annotation_path)

        # store results in an easier format
        self.__tremolo_res = dict()
        for target in self.__tremolo_targets[0]:
            self.__tremolo_res[target] = {"hits": [], "domains": []}
            for hit in self.__tremolo_targets[0][target]:
                hit_res = dict()
                for key in ["Aligned_cols", "E-value", "Identities", "Probab", 
                            "Qali", "Qsize", "Qstart", "Qstop", 
                            "Score", "Similarity", "Sum_probs",
                            "Tali", "Tsize", "Tstart", "Tstop", 
                            "descr"]:
                    hit_res[key] = self.__tremolo_targets[0][target][hit].get(key, "")
                self.__tremolo_res[target]["hits"].append(hit_res)
            self.__tremolo_res[target]["domains"] = self.__tremolo_annotation[target]
        return self.__tremolo_res
    
    def save_tremolo(self, outputfile):
        """ save tremolo results to a text file
        
        Parameter
        ---------
        outputfile: string
            the path to the output file
            
        Usage
        -----
        
        >>> hca.save_tremolo("tremolo_results.txt")
        >>>
        """
        groups = group_resda(self.__tremolo_targets, self.__tremolo_annotation)
        write_tremolo_results(self.__tremolo_query, self.__tremolo_domains, 
                              self.__tremolo_targets, self.__tremolo_annotation, 
                              groups, outputfile)
    
def check_seq_type(seq):
    """ check that sequence is of correct type
    """
    queryname = "query"
    if isinstance(seq, str):
        seq = seq
    elif isinstance(seq, Seq.Seq):
        seq = str(seq)
    elif isinstance(seq, Bio.SeqRecord.SeqRecord):
        queryname = seq.id
        seq = str(seq.seq)
    else:
        raise TypeError("Error, the 'seq' argument must be either a string a Bio.Seq instance or a Bio.SeqRecord instance")
    return seq, queryname

def check_if_msa(querynames, sequences):
    """ check if provided sequences are from an MSA
    if yes, transform msa sequences to ungapped sequences for HCA analysis
    """
    msa_length, msa_seq, new_sequences = list(), list(), list()
    is_msa = False
    if hasattr(IUPAC, "protein"):
        prot_alphabet = set(IUPAC.protein.letters)
    else:
        prot_alphabet = set(IUPAC.protein_letters)
    gaps = set(["-", "."])
    for i, seq in enumerate(sequences):
        new_seq = ""
        for j, c in enumerate(seq):
            if c in prot_alphabet:
                new_seq += c
            else:
                if c in gaps:
                    is_msa = True
                else:
                    print("Invalid amino acids ({}, {}) in protein {}, replaced by X".format(j, c, querynames[i]), file=sys.stderr)
        new_sequences.append(new_seq)
        msa_seq.append(seq)
        msa_length.append(len(seq))
        
    if is_msa:
        # chec identical sequence lengths
        if len(set(msa_length)) != 1:
            raise ValueError("Error, MSA symbols detected but sequences have different lengths")
    return is_msa, msa_seq, new_sequences

