#!/usr/bin/env python
# Imports
import unittest
import os
import sys
import tempfile
import filecmp
import shutil

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'bin', 'enrichm')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')]+sys.path

from enrichm.classify_checks import ClassifyChecks
from enrichm.parser import RulesJson, RulesJsonVersion0, Parser

###############################################################################

class Tests(unittest.TestCase):

    def test_basic(self):
        parser = Parser
        rules = RulesJson()
        basic_rules = os.path.join(path_to_data, "rules/basic_rules_1.json")
        basic_gff = os.path.join(path_to_data, "gffs/basic_gff_1.gff")

        basic_rules_object = rules.load(basic_rules)
        classify_checks = ClassifyChecks(basic_rules_object)
        features, _ =  parser.parse_gff(basic_gff)

        self.assertFalse(classify_checks.check("sample_pathway", features['genome_1']))
        self.assertFalse(classify_checks.check("sample_pathway", features['genome_2']))
        self.assertTrue(classify_checks.check("sample_pathway", features['genome_3']))

    def test_split(self):
        parser = Parser
        rules = RulesJson()
        basic_rules = os.path.join(path_to_data, "rules/basic_rules_2.json")
        basic_gff = os.path.join(path_to_data, "gffs/basic_gff_2.gff")
        basic_rules_object = rules.load(basic_rules)
        classify_checks = ClassifyChecks(basic_rules_object)
        features, _ =  parser.parse_gff(basic_gff)

        self.assertTrue(classify_checks.check("sample_pathway", features['genome_1']))
        self.assertFalse(classify_checks.check("sample_pathway", features['genome_2']))

    def test_missing_genes(self):
        parser = Parser
        rules = RulesJson()
        basic_rules = os.path.join(path_to_data, "rules/basic_rules_3.json")
        basic_gff = os.path.join(path_to_data, "gffs/basic_gff_3.gff")
        basic_rules_object = rules.load(basic_rules)
        classify_checks = ClassifyChecks(basic_rules_object)
        features, _ =  parser.parse_gff(basic_gff)

        self.assertTrue(classify_checks.check("sample_pathway", features['genome_1']))
        self.assertFalse(classify_checks.check("sample_pathway", features['genome_2']))

    @unittest.skip("Known failure")
    def test_split_contig(self):
        parser = Parser
        rules = RulesJson()
        basic_rules = os.path.join(path_to_data, "rules/basic_rules_4.json")
        basic_gff = os.path.join(path_to_data, "gffs/basic_gff_4.gff")
        basic_rules_object = rules.load(basic_rules)
        classify_checks = ClassifyChecks(basic_rules_object)
        features, _ =  parser.parse_gff(basic_gff)

        self.assertTrue(classify_checks.check("sample_pathway", features['genome_1']))

    # TODO: Check overlapping contig names

if __name__ == "__main__":
    unittest.main()
