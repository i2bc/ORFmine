from __future__ import division, absolute_import, print_function

import os, sys, shutil
import numpy as np
import tempfile
import unittest

from pyHCA import HCA
from Bio.SeqRecord import SeqRecord

def random_seq(X):
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    return "".join(np.random.choice(aas, X))

def tmp_fasta(N):
    fileid, filepath = tempfile.mkstemp()
    sequences = list()
    for n in range(N):
        seq = random_seq(np.random.randint(50, 200))
        sequences.append(seq)
        txt = ">query_{}\n{}\n".format(n, seq)
        os.write(fileid, txt.encode())
    os.close(fileid)
    return filepath, sequences

class TestHCA(unittest.TestCase):
    """ test input arguments of the HCA class
    """
    def setup(self):
        pass
    
    def test_initseqstr(self):
        """ input is a sequences
        """
        for x in range(10):
            seq = random_seq(np.random.randint(50, 200))
            hca = HCA(seq=seq)
            self.assertEqual(hca._number_of_sequences, 1)
            self.assertEqual(hca.sequences[0], seq)
    
    def test_initseqlist(self):
        """ input is a list of sequences
        """
        for x in range(10):
            s = np.random.randint(2, 10)
            sequences = [random_seq(np.random.randint(50, 200)) for n in range(s)]
            hca = HCA(seq=sequences)
            self.assertEqual(hca._number_of_sequences, s)
        
    def test_initseqrec(self):
        """ input is a SeqRecord
        """
        for x in range(10):
            seq = random_seq(np.random.randint(50, 200))
            rec = SeqRecord(seq, "query_1")
            hca = HCA(seq=seq)
            self.assertEqual(hca._number_of_sequences, 1)
            self.assertEqual(hca.sequences[0], seq)
    
    def test_initseqlistrec(self):
        """ input is a list of SeqRecord
        """
        for x in range(10):
            sequences = list()
            s = np.random.randint(2, 10)
            for n in range(s):
                seq = random_seq(np.random.randint(50, 200))
                rec = SeqRecord(seq, "query_{}".format(n+1))
                sequences.append(rec)
            hca = HCA(seq=sequences)
            self.assertEqual(hca._number_of_sequences, s)
    
    def test_initseqfile(self):
        """ input is a fasta file
        """
        s = np.random.randint(2, 10)
        fasta_file, sequences = tmp_fasta(s)
        hca = HCA(seqfile=fasta_file)
        self.assertEqual(hca._number_of_sequences, s)
    
    #def test_initseqfile_err(self):
        #""" input is a fasta file
        #"""
        #s = np.random.randint(2, 10)
        #fasta_file, sequences = tmp_fasta(s)
        #try:
            #hca = HCA(seq=fasta_file)
            #self.assertEqual(hca._number_of_sequences, 0)
        #except:
            #pass
    
    
if __name__ == "__main__":
    unittest.main()
    

