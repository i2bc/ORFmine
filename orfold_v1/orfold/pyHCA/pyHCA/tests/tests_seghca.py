from __future__ import division, absolute_import, print_function

import os, sys, shutil
import numpy as np
import unittest

from pyHCA import HCA


class TestSegHCA(unittest.TestCase):
    
    def setUp(self):
        self.check_sequence = "NGRHTGFGRTCCDKGADHLKGEGHCCITLAKRGYFPCEPWCTLLFALNMFNMQNMMRQQFSDDHNNMGRLCQQTTHRFPFNSDNKEEYIWLYKVQRLGAW"
        self.check_domains = {0: [6, 100, 0.00222224]}        
        self.check_clusters = {0: [ 6, 12,  np.array([1, 0, 0, 0, 1, 1])],
                               1: [18, 19,  np.array([1])],
                               2: [24, 29,  np.array([1, 1, 1, 0, 1])],
                               3: [33, 35,  np.array([1, 1])],
                               4: [36, 37,  np.array([1])],
                               5: [39, 60,  np.array([1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1])],
                               6: [66, 71,  np.array([1, 0, 0, 1, 1])],
                               7: [77, 78,  np.array([1])],
                               8: [79, 80,  np.array([1])],
                               9: [87, 100, np.array([1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1])],
                               }
    def assert_domain(self, i, dom):
        self.assertEqual(self.check_domains[i][0], dom.start)
        self.assertEqual(self.check_domains[i][1], dom.stop)
        if np.isnan(dom.pvalue):
            assert(np.isnan(self.check_domains[i][2]))
        else:
            assert(np.allclose(self.check_domains[i][2], dom.pvalue))
        
    def test_domains(self):
        hca = HCA(seq=self.check_sequence)
        hca.segments()
        domains = hca.get_domains()
        self.assertEqual(len(domains), len(self.check_domains))
        for i, dom in enumerate(domains):
            self.assert_domain(i, dom)
    
    def assert_cluster(self, i, cluster):
        self.assertEqual(self.check_clusters[i][0], cluster.start)
        self.assertEqual(self.check_clusters[i][1], cluster.stop)
        assert(np.all(self.check_clusters[i][2] == cluster.hydro_cluster))
        
    def test_clusters(self):
        hca = HCA(seq=self.check_sequence)
        hca.segments()
        clusters = hca.get_clusters()
        self.assertEqual(len(clusters), len(self.check_clusters))
        for i, cluster in enumerate(clusters):
            self.assert_cluster(i, cluster)
    
    
if __name__ == "__main__":
    unittest.main()
