#!/usr/bin/env python
# Imports
import unittest
import filecmp
import tempfile
import os
import sys

PATH_TO_DATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from enrichm.network_analyzer import NetworkAnalyser

class Tests(unittest.TestCase):

    abundance_metadata = "abundances.metadata.tsv"
    abundance = "abundances.tsv"
    genome_annotations = "genome.annotations.tsv"
    genome_metadata = "genome.metadata.tsv"
    queries_file = "queries.txt"
    network_data = "network_data"
    expected_explore_output = "explore_test"
    expected_pathway_output = "pathway_test"

    abundance_metadata_path = os.path.join(PATH_TO_DATA, network_data, abundance_metadata)
    abundance_path = os.path.join(PATH_TO_DATA, network_data, abundance)
    expected_explore_output_path = os.path.join(PATH_TO_DATA, network_data, expected_explore_output)
    genome_annotations_path = os.path.join(PATH_TO_DATA, network_data, genome_annotations)
    genome_metadata_path = os.path.join(PATH_TO_DATA, network_data, genome_metadata)
    expected_pathway_output_path = os.path.join(PATH_TO_DATA, network_data, expected_pathway_output)
    queries_file_path = os.path.join(PATH_TO_DATA, network_data, queries_file)

    def test_hello_pathway(self):
        '''
        Testing normal pathway workflow, no abundance or transcriptome or metabolome data included
        '''

        tmp = tempfile.mkdtemp()

        network_analyser = NetworkAnalyser()
        network_analyser.network_pipeline(NetworkAnalyser.PATHWAY,
                                          self.genome_annotations_path, self.genome_metadata_path,
                                          None, None, # No transcriptome data
                                          None, None, # No metagenome data
                                          None, # No metabolome data
                                          None, # No fisher results
                                          None, # No depth (explore runs only)
                                          list(), # Dont filter anything out
                                          ["map00860"], # Filter to keep it small.
                                          None, # No queries
                                          tmp # Output directory
                                          )
        

        expected_files = sorted(os.listdir(tmp))
        observed_files = sorted(os.listdir(self.expected_pathway_output_path))

        self.assertEqual(len(expected_files), len(observed_files))

        for expected_file, observed_file in zip(expected_files, observed_files):
            expected_file_path = os.path.join(tmp, expected_file)
            observed_file_path = os.path.join(self.expected_pathway_output_path, observed_file)
            # Are all files present?
            self.assertEqual(expected_file, observed_file)

            # Do the all look the same on the inside?
            with open(expected_file_path) as expected_file_io:

                with open(observed_file_path) as observed_file_io:
                    expected_lines = sorted(expected_file_io.readlines())
                    observed_lines = sorted(observed_file_io.readlines())
                    self.assertEqual(len(expected_lines), len(observed_lines))

                    for expected_line, observed_line in zip(expected_lines, observed_lines):
                        self.assertEqual(expected_line, observed_line)

    def test_hello_query(self):
        '''
        Testing normal explore workflow, no abundance or transcriptome or metabolome data included
        '''

        tmp = tempfile.mkdtemp()

        network_analyser = NetworkAnalyser()
        network_analyser.network_pipeline(NetworkAnalyser.EXPLORE,
                                          self.genome_annotations_path, self.genome_metadata_path,
                                          None, None, # No transcriptome data
                                          None, None, # No metagenome data
                                          None, # No metabolome data
                                          None, # No fisher results
                                          1, # No depth (explore runs only)
                                          list(), # Dont filter anything out
                                          list(), # Filter to keep it small.
                                          self.queries_file_path, # No queries
                                          tmp # Output directory
                                          )

        expected_files = sorted(os.listdir(tmp))
        observed_files = sorted(os.listdir(self.expected_explore_output_path))

        self.assertEqual(len(expected_files), len(observed_files))

        for expected_file, observed_file in zip(expected_files, observed_files):
            expected_file_path = os.path.join(tmp, expected_file)
            observed_file_path = os.path.join(self.expected_explore_output_path, observed_file)
            # Are all files present?
            self.assertEqual(expected_file, observed_file)

            # Do the all look the same on the inside?
            with open(expected_file_path) as expected_file_io:

                with open(observed_file_path) as observed_file_io:
                    expected_lines = sorted(expected_file_io.readlines())
                    observed_lines = sorted(observed_file_io.readlines())
                    self.assertEqual(len(expected_lines), len(observed_lines))

                    for expected_line, observed_line in zip(expected_lines, observed_lines):
                        self.assertEqual(expected_line, observed_line)

if __name__ == "__main__":
    unittest.main()
