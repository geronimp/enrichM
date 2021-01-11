#!/usr/bin/env python
# Imports
import unittest
import os.path
import sys
import subprocess
import tempfile
import filecmp

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','enrichm')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_annotate = os.path.join(path_to_data, 'enrichm_annotate')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from enrichm.enrichment import Enrichment

###############################################################################

class Tests(unittest.TestCase):

    genome_annotation_simple_example = {"genome_1": {"K00001":1,
                                                     "K00002":2},
                                        "genome_2": {"K00003":1},
                                        "genome_3": {"K00001":5,
                                                     "K00002":4,
                                                     "K00003":5}}

    genome_groups_simple_example = {"group_1": ["genome_1"],
                                "group_2": ["genome_2", "genome_3"]}

    sample_abundance = {"sample_1": {"genome_1":1.0, "genome_2":0.5, "genome_3":3.0},
                        "sample_2": {"genome_1":0.5, "genome_2":1.2, "genome_3":5.0},
                        "sample_3": {"genome_1":0.1, "genome_2":1.1, "genome_3":6.0},
                        "sample_4": {"genome_1":5.0, "genome_2":5.2, "genome_3":0.2},
                        "sample_5": {"genome_1":6.0, "genome_2":4.9, "genome_3":0.1},
                        "sample_6": {"genome_1":7.0, "genome_2":5.0, "genome_3":0.0}}

    sample_groups = {"sample_group_1": ["sample_1", "sample_2", "sample_3"],
                     "sample_group_2": ["sample_4", "sample_5", "sample_6"]}

    genomes = ["genome_1","genome_2","genome_3"]
    annotations = ["K00001", "K00002", "K00003"]
    enrichment_test_object = Enrichment()

    def test_check_annotation_type(self):
        pfam = ['PF10117']
        self.assertEqual(
            Enrichment().check_annotation_type(pfam), Enrichment.PFAM)
        cazy = ['GH42']
        self.assertEqual(
            Enrichment().check_annotation_type(cazy), Enrichment.CAZY)
        tigrfam = ['TIGR00008']
        self.assertEqual(
            Enrichment().check_annotation_type(tigrfam), Enrichment.TIGRFAM)
        ko = ['K00399']
        self.assertEqual(
            Enrichment().check_annotation_type(ko), Enrichment.KEGG)
        ec = ['1.2.3.4']
        self.assertEqual(
            Enrichment().check_annotation_type(ec), Enrichment.EC)

    def test_calculate_portions(self):
        expected = [['Annotation', 'group_1', 'group_2'],
                    ['K00001', '1.0', '0.5'],
                    ['K00002', '1.0', '0.5'],
                    ['K00003', '0.0', '1.0']]
        result = self.enrichment_test_object.calculate_portions(self.annotations,
                                                                self.genome_groups_simple_example,
                                                                self.genome_annotation_simple_example,
                                                                self.genomes,
                                                                1)
        self.assertEqual(result, expected)

    def test_enrichment_from_ko_matrix(self):

        tmp = tempfile.mkdtemp()
        expected_output = os.path.join(path_to_data, 'enrichm_enrichment_ko')
        metadata        = os.path.join(path_to_data, 'metadata.tsv')
        cmd             = '%s enrichment --annotate_output %s --metadata %s --output %s --force --ko --verbosity 1' \
                            % (path_to_script, path_to_annotate, metadata, tmp)
        subprocess.call(cmd, shell=True)

        self.assertTrue(filecmp.dircmp(tmp, expected_output))
        # The pvalues are never exact - cannot compare files directly
        #for file in os.listdir(tmp):
        #    if file         == 'enrichment.log': continue
        #    output_file     = os.path.join(tmp, file)
        #    expected_file   = os.path.join(expected_output, file)
        #    self.assertTrue(filecmp.cmp(output_file, expected_file))

    def test_enrichment_from_pfam_matrix(self):

        tmp = tempfile.mkdtemp()
        expected_output = os.path.join(path_to_data, 'enrichm_enrichment_pfam')
        metadata        = os.path.join(path_to_data, 'metadata.tsv')
        cmd             = '%s enrichment --annotate_output %s --metadata %s --output %s --force --pfam  --verbosity 1' \
                            % (path_to_script, path_to_annotate, metadata, tmp)

        subprocess.call(cmd, shell=True)

        self.assertTrue(filecmp.dircmp(tmp, expected_output))

        # The pvalues are never exact - cannot compare files directly
        #for file in os.listdir(tmp):
        #    if file         == 'enrichment.log': continue
        #    output_file     = os.path.join(tmp, file)
        #    expected_file   = os.path.join(expected_output, file)
        #    self.assertTrue(filecmp.cmp(output_file, expected_file))

    def test_weight_annotation_matrix(self):
        expected = {'sample_group_1': {'K00001': [16.0, 25.5, 30.1],
                                       'K00002': [14.0, 21.0, 24.2],
                                       'K00003': [15.5, 26.2, 31.1]},
                    'sample_group_2': {'K00001': [6.0, 6.5, 7.0],
                                       'K00002': [10.8, 12.4, 14.0],
                                       'K00003': [6.2, 5.4, 5.0]}}

        result = self.enrichment_test_object.weight_annotation_matrix(self.sample_abundance,
                                                             self.genome_annotation_simple_example,
                                                             self.sample_groups,
                                                             self.annotations)
        self.assertEqual(result, expected)

    def test_operon_enrichment(self):
        pass
        

if __name__ == "__main__":
    unittest.main()
