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
from enrichm.genome import Genome, Annotation, AnnotationParser

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
        subprocess.check_call(cmd, shell=True)

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
        subprocess.check_call(cmd, shell=True)

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
        
    @unittest.skip("No input data for this test")
    def test_very_simple_orthology(self):
        tmp = tempfile.mkdtemp()
        protein_file = os.path.join(path_to_data, 'test_protein_bin', "GCF_001889405.1_ASM188940v1_subset.faa")

        cmd = '%s annotate \
                        --threads 4 \
                        --orthologs \
                        --protein_files %s \
                        --output %s \
                        --force' % (path_to_script, protein_file, tmp)
        subprocess.check_call(cmd, shell=True)

    def test_very_simple_homology(self):
        tmp = tempfile.mkdtemp()
        protein_file = os.path.join(path_to_data, 'test_protein_bin', "GCF_001889405.1_ASM188940v1_subset.faa")

        cmd = '%s annotate \
                        --threads 4 \
                        --cluster \
                        --protein_files %s \
                        --output %s \
                        --force' % (path_to_script, protein_file, tmp)
        subprocess.check_call(cmd, shell=True)

    def test_genome_instances_have_isolated_sequences(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            genome_one = os.path.join(tmp_dir, "genome_one.faa")
            genome_two = os.path.join(tmp_dir, "genome_two.faa")

            with open(genome_one, 'w', encoding='utf-8') as handle:
                handle.write(">shared_seq\nMAAAA\n")
            with open(genome_two, 'w', encoding='utf-8') as handle:
                handle.write(">shared_seq\nMTTTT\n")

            genome_a = Genome(False, None, genome_one, None)
            genome_b = Genome(False, None, genome_two, None)

            self.assertIsNot(genome_a.sequences, genome_b.sequences)
            self.assertEqual(set(genome_a.sequences.keys()), {"shared_seq"})
            self.assertEqual(set(genome_b.sequences.keys()), {"shared_seq"})

            seq_a = next(iter(genome_a.sequences.values()))
            seq_b = next(iter(genome_b.sequences.values()))
            seq_a.annotations.append(Annotation("K00001", 0.0, {0}, AnnotationParser.KO))

            self.assertEqual(len(seq_a.annotations), 1)
            self.assertEqual(len(seq_b.annotations), 0)

    def test_genome_loads_sequences_in_order(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            protein_file = os.path.join(tmp_dir, "test_genome.faa")
            with open(protein_file, 'w') as fh:
                fh.write(">seq_a\nMAAAK\n>seq_b\nMTTTK\n>seq_c\nMCCCK\n")

            genome = Genome(False, None, protein_file, None)

            self.assertEqual(len(genome.sequences), 3)
            self.assertIn('seq_a', genome.sequences)
            self.assertIn('seq_b', genome.sequences)
            self.assertIn('seq_c', genome.sequences)
            self.assertEqual(genome.protein_ordered_dict, {0: 'seq_a', 1: 'seq_b', 2: 'seq_c'})

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
