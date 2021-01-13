#!/usr/bin/env python
# Imports
import unittest
import os.path
import sys
import subprocess
import tempfile

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'bin', 'enrichm')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')]+sys.path

from enrichm.annotate import Annotate

###############################################################################

class Tests(unittest.TestCase):

    def test_hello_world_nucleic_dir(self):
        tmp = tempfile.mkdtemp()
        bin = os.path.join(path_to_data, 'test_nucleic_bin')
        cmd = '%s annotate \
                        --threads 10 \
                        --ko \
                        --pfam \
                        --tigrfam \
                        --genome_directory %s \
                        --output %s \
                        --force' % (path_to_script, bin, tmp)
        subprocess.call(cmd, shell=True)

    def test_hello_world_protein_dir(self):
        tmp = tempfile.mkdtemp()
        bin = os.path.join(path_to_data, 'test_protein_bin')

        cmd = '%s annotate \
                        --threads 10 \
                        --ko \
                        --pfam \
                        --tigrfam \
                        --protein_directory %s \
                        --output %s \
                        --force' % (path_to_script, bin, tmp)
        subprocess.call(cmd, shell=True)

    def test_hello_world_nucleic_file(self):
        tmp = tempfile.mkdtemp()
        genome_file = os.path.join(path_to_data, 'test_nucleic_bin', "GCF_001889405.1_ASM188940v1_subset.fna")
        cmd = '%s annotate \
                        --threads 10 \
                        --ko \
                        --pfam \
                        --tigrfam \
                        --genome_files %s \
                        --output %s \
                        --force' % (path_to_script, genome_file, tmp)
        subprocess.check_call(cmd, shell=True)

    def test_hello_world_protein_file(self):
        tmp = tempfile.mkdtemp()
        protein_file = os.path.join(path_to_data, 'test_protein_bin', "GCF_001889405.1_ASM188940v1_subset.faa")

        cmd = '%s annotate \
                        --threads 10 \
                        --ko \
                        --pfam \
                        --tigrfam \
                        --protein_files %s \
                        --output %s \
                        --force' % (path_to_script, protein_file, tmp)
        subprocess.check_call(cmd, shell=True)

    def test_hello_world_hypothetical_cluster(self):
        tmp = tempfile.mkdtemp()
        protein_file = os.path.join(path_to_data, 'test_protein_bin', "GCF_001889405.1_ASM188940v1_subset.faa")

        cmd = '%s annotate \
                        --verbosity 5 \
                        --threads 4 \
                        --clusters \
                        --protein_files %s \
                        --output %s \
                        --force' % (path_to_script, protein_file, tmp)
        subprocess.check_call(cmd, shell=True)

    def test_very_simple_orthology(self):
        tmp = tempfile.mkdtemp()
        protein_file = os.path.join(path_to_data, 'cluster_data', "very_simple_test.faa")

        cmd = '%s annotate \
                        --threads 4 \
                        --orthologs \
                        --protein_files %s \
                        --output %s \
                        --force' % (path_to_script, protein_file, tmp)
        subprocess.check_call(cmd, shell=True)

    def test_very_simple_homology(self):
        tmp = tempfile.mkdtemp()
        protein_file = os.path.join(path_to_data, 'cluster_data', "very_simple_test.faa")

        cmd = '%s annotate \
                        --threads 4 \
                        --cluster \
                        --protein_files %s \
                        --output %s \
                        --force' % (path_to_script, protein_file, tmp)
        subprocess.check_call(cmd, shell=True)

    def test(self):
        tmp = tempfile.mkdtemp()
        self.simple_annotate_instance \
            = Annotate(tmp,
                       True, True, True, True, True, True, True, True, True,# Annotate with all databases
                       1e-05, 0, 0.3, 0.7, 0.7, 0.7, # Runtime options
                       False, False, False, False, False, False, True, # Cutoffs
                       5, # Inflation
                       4, 2500, # chunks
                       False, 1, 1, '.fna', False)

if __name__ == "__main__":
    unittest.main()
