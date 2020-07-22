#!/usr/bin/env python
# Imports
import unittest
import os.path
import sys
import subprocess
import tempfile
import filecmp
from scipy import stats


path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','enrichm')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_annotate = os.path.join(path_to_data, 'enrichm_annotate')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from enrichm.enrichment import Test
from enrichm.databases import Databases

###############################################################################

class Tests(unittest.TestCase):
    
    genome_annotation_simple_example = {"genome_1": {"K00001":1,
                                                     "K00002":2,
                                                     "K00003":0},
                                        "genome_2": {"K00001":0,
                                                     "K00002":0,
                                                     "K00003":1},
                                        "genome_3": {"K00001":5,
                                                     "K00002":4,
                                                     "K00003":5}}
    genome_groups_simple_example = {"group_1": ["genome_1"],
                                    "group_2": ["genome_2", "genome_3"]}
    simple_test_object = Test(genome_annotation_simple_example,
                              genome_groups_simple_example,
                              "kegg",
                              0.1,
                              'fdr_bh',
                              1,
                              Databases())
    sample_to_annotation = {'sample_group_1': {'K00001': [16.0, 25.5, 30.1],
                                               'K00002': [14.0, 21.0, 24.2],
                                               'K00003': [15.5, 26.2, 31.1]},
                            'sample_group_2': {'K00001': [6.0, 6.5, 7.0],
                                               'K00002': [10.8, 12.4, 14.0],
                                               'K00003': [6.2, 5.4, 5.0]}}
    annotations = ["K00001", "K00002", "K00003"]

    def test_test_chooser(self):
        groups_1 = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]]
        groups_2 = [[1], [1, 2, 3, 4, 5]]
                
        test_instance_1 = self.simple_test_object.test_chooser(groups_1)
        test_instance_2 = self.simple_test_object.test_chooser(groups_2)
            
        self.assertEqual(test_instance_1[0], stats.fisher_exact)
        self.assertEqual(test_instance_1[1], stats.mannwhitneyu)
        self.assertEqual(test_instance_2[0], self.simple_test_object.PA)
        self.assertEqual(test_instance_2[1], stats.norm.cdf)

    def test_count(self):
        '''
        test both frequency and presence absence counting in Test.
        '''
        
        self.assertEqual(self.simple_test_object.count("K00001", "group_1", False), (1,0))
        self.assertEqual(self.simple_test_object.count("K00001", "group_1", True), ([1],0))
        self.assertEqual(self.simple_test_object.count("K00001", "group_2", False), (1,1))
        self.assertEqual(self.simple_test_object.count("K00001", "group_2", True), ([0, 5], 0))
        self.assertEqual(self.simple_test_object.count("K00003", "group_2", False), (2, 0))
        self.assertEqual(self.simple_test_object.count("K00003", "group_2", True), ([1, 5], 0))

    def test_gene_frequencies(self):
        expect_1 = [['K00002', 'group_1', 'group_2', [[2], 0], [[0, 4], 0]],
                    ['K00003', 'group_1', 'group_2', [[0], 0], [[1, 5], 0]],
                    ['K00001', 'group_1', 'group_2', [[1], 0], [[0, 5], 0]]]
        expect_2 = [['K00002', 'group_1', 'group_2', [1, 0], [1, 1]],
                    ['K00003', 'group_1', 'group_2', [0, 1], [2, 0]],
                    ['K00001', 'group_1', 'group_2', [1, 0], [1, 1]]]
        for result in self.simple_test_object.gene_frequencies("group_1", "group_2", True):
            if result in expect_1:
                expect_1.pop(expect_1.index(result))
        self.assertEqual(expect_1, list())
        
        for result in self.simple_test_object.gene_frequencies("group_1", "group_2", False):
            if result in expect_2:
                expect_2.pop(expect_2.index(result))
        self.assertEqual(expect_2, list())
    
    def test_test_weighted_abundances(self):
        expect = [[[['annotation', 'group_1', 'group_2', 'enriched_in', 'group_1_mean', 'group_2_mean', 'score', 'pvalue', 'corrected_pvalue', 'description'],
                    ['K00001', 'sample_group_1', 'sample_group_2', 'sample_group_1', '23.866666666666664', '6.5', 0.0, 0.04042779918502612, '0.060591636418731595', 'E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]'],
                    ['K00002', 'sample_group_1', 'sample_group_2', 'sample_group_1', '19.733333333333334', '12.4', 0.5, 0.060591636418731595, '0.060591636418731595', 'AKR1A1, adh; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]'],
                    ['K00003', 'sample_group_1', 'sample_group_2', 'sample_group_1', '24.26666666666667', '5.533333333333334', 0.0, 0.04042779918502612, '0.060591636418731595', 'hom; homoserine dehydrogenase [EC:1.1.1.3]']],
                     'sample_group_1_vs_sample_group_2_gvg_results.mannwhitneyu.tsv']]
        
        result = self.simple_test_object.test_weighted_abundances(self.sample_to_annotation,
                                                                  self.annotations)

        self.assertEqual(expect, result)

if __name__ == "__main__":
    unittest.main()
