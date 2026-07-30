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
from enrichm.genome import Genome, Annotation, AnnotationParser, Sequence

###############################################################################

class Tests(unittest.TestCase):

    def test_hello_world_nucleic_dir(self):
        tmp = tempfile.mkdtemp()
        bin = os.path.join(path_to_data, 'test_nucleic_bin')
        cmd = '%s annotate \
                        --threads 10 \
                        --ko_hmm \
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
                        --ko_hmm \
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
                        --ko_hmm \
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
                        --ko_hmm \
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
        
    @unittest.skip("Requires mcl which is not available in this environment")
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

    def test_overlapping_pfam_annotations_from_different_clans_are_kept_once(self):
        # Overlapping domains from different clans are all retained, but each hit
        # must only ever be added once. Adding it once per overlapping annotation
        # makes the annotation list grow exponentially with the number of hits.
        sequence = Sequence("seq_a")
        pfam2clan = {'PF%05i' % hit: 'CL%05i' % hit for hit in range(20)}

        for hit in range(20):
            sequence.add(['PF%05i.1' % hit], 1e-10, range(0, 100),
                         AnnotationParser.PFAM, pfam2clan=pfam2clan)

        self.assertEqual(len(sequence.annotations), 20)

    def test_pfam_annotation_beaten_within_its_clan_is_discarded(self):
        # An overlap with an unrelated clan must not rescue an annotation that
        # has already lost to a better description of the same region.
        sequence = Sequence("seq_a")
        pfam2clan = {'PF00001': 'CL0001', 'PF00002': 'CL0002', 'PF00003': 'CL0001'}

        sequence.add(['PF00001.1'], 1e-30, range(0, 100), AnnotationParser.PFAM, pfam2clan=pfam2clan)
        sequence.add(['PF00002.1'], 1e-10, range(50, 200), AnnotationParser.PFAM, pfam2clan=pfam2clan)
        # Overlaps both, but loses to PF00001.1 inside its own clan.
        sequence.add(['PF00003.1'], 1e-02, range(20, 120), AnnotationParser.PFAM, pfam2clan=pfam2clan)

        self.assertEqual(sorted(annotation.annotation for annotation in sequence.annotations),
                         ['PF00001.1', 'PF00002.1'])

    def test_annotation_must_beat_every_overlap_to_be_kept(self):
        # Annotation types without clans compete on e-value alone, so beating
        # only one of several overlapping annotations is not enough.
        sequence = Sequence("seq_a")

        sequence.add(['K00001'], 1e-30, range(0, 50), AnnotationParser.KO)
        sequence.add(['K00002'], 1e-02, range(60, 100), AnnotationParser.KO)
        # Better than K00002, worse than K00001, overlapping both.
        sequence.add(['K00003'], 1e-10, range(40, 70), AnnotationParser.KO)

        self.assertEqual(sorted(annotation.annotation for annotation in sequence.annotations),
                         ['K00001', 'K00002'])

    def test_overlapping_pfam_annotations_from_same_clan_keep_best_evalue(self):
        sequence = Sequence("seq_a")
        pfam2clan = {'PF00001': 'CL0001', 'PF00002': 'CL0001'}

        sequence.add(['PF00001.1'], 1e-10, range(0, 100), AnnotationParser.PFAM, pfam2clan=pfam2clan)
        sequence.add(['PF00002.1'], 1e-20, range(10, 90), AnnotationParser.PFAM, pfam2clan=pfam2clan)
        sequence.add(['PF00001.1'], 1e-05, range(20, 80), AnnotationParser.PFAM, pfam2clan=pfam2clan)

        self.assertEqual([annotation.annotation for annotation in sequence.annotations], ['PF00002.1'])

    def test_pfam_annotations_outside_of_clans_are_kept(self):
        # Pfams that belong to no clan are not alternative annotations of each
        # other, so overlapping hits are all retained.
        sequence = Sequence("seq_a")

        sequence.add(['PF00001.1'], 1e-10, range(0, 100), AnnotationParser.PFAM, pfam2clan={})
        sequence.add(['PF00002.1'], 1e-20, range(10, 90), AnnotationParser.PFAM, pfam2clan={})

        self.assertEqual(sorted(annotation.annotation for annotation in sequence.annotations),
                         ['PF00001.1', 'PF00002.1'])

    def test_called_proteins_are_uniquely_named(self):
        # Genome objects and HMM/DIAMOND searches key sequences on the first word
        # of the FASTA description, so gene ids must be unique within it.
        with tempfile.TemporaryDirectory() as tmp_dir:
            annotate = Annotate(tmp_dir,
                                False, False, True, False, False, False, False, False, False,
                                1e-05, 0, 0.3, 0.7, 0.7, 0.7,
                                False, False, False, False, False, False, True,
                                5, 4, 2500, False, 1, 1, '.fna', False)
            annotate.call_proteins(os.path.join(path_to_data, 'test_nucleic_bin'))

            protein_file = os.path.join(tmp_dir, Annotate.GENOME_PROTEINS,
                                        'GCF_001889405.1_ASM188940v1_subset.faa')
            descriptions = [line[1:].strip() for line in open(protein_file) if line.startswith('>')]
            names = [description.partition(' ')[0] for description in descriptions]

            self.assertGreater(len(names), 1)
            self.assertEqual(len(set(names)), len(names))

            genome = Genome(False, None, protein_file, None)
            self.assertEqual(len(genome.sequences), len(names))

            # Coordinates are parsed out of the description, and are needed to
            # write .gff files.
            sequence = genome.sequences[names[0]]
            self.assertEqual(int(sequence.finishpos) > int(sequence.startpos), True)
            self.assertIn(sequence.direction, ('1', '-1'))

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
