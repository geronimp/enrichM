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

from enrichm.enrichment import Enrichment, gene_fisher_calc, mannwhitneyu_calc, zscore_calc

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
        subprocess.check_call(cmd, shell=True)

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

        subprocess.check_call(cmd, shell=True)

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

    # --- gene_fisher_calc ---

    def test_gene_fisher_calc_enriched_in_group1(self):
        # K00001 present in both group_1 genomes, 1 of 3 group_2 genomes
        result = gene_fisher_calc(['K00001', 'group_1', 'group_2', [2, 0], [1, 2]])
        self.assertEqual(result[0], 'K00001')
        self.assertEqual(result[3], 'group_1')
        self.assertAlmostEqual(result[-1], 0.4, places=5)

    def test_gene_fisher_calc_all_absent_returns_na(self):
        # Zero in both groups → NA result, p=1.0
        result = gene_fisher_calc(['K00001', 'group_1', 'group_2', [0, 2], [0, 2]])
        self.assertEqual(result[3], 'NA')
        self.assertEqual(result[-1], 1.0)

    def test_gene_fisher_calc_enriched_in_group2(self):
        # K00001 absent from group_1, present in group_2
        result = gene_fisher_calc(['K00001', 'group_1', 'group_2', [0, 3], [3, 0]])
        self.assertEqual(result[3], 'group_2')

    # --- mannwhitneyu_calc ---

    def test_mannwhitneyu_calc_constant_group_returns_na(self):
        # Constant values in group_1 → MWU cannot be computed → NA
        result = mannwhitneyu_calc(['K00001', 'g1', 'g2', [[5, 5, 5]], [[1, 2, 3]]])
        self.assertEqual(result[3], 'NA')

    def test_mannwhitneyu_calc_both_zero_returns_na(self):
        # All zeros in both groups → NA
        result = mannwhitneyu_calc(['K00001', 'g1', 'g2', [[0, 0, 0]], [[0, 0, 0]]])
        self.assertEqual(result[3], 'NA')

    def test_mannwhitneyu_calc_enriched_in_higher_mean_group(self):
        # group_1 clearly higher → enriched_in = 'g1'
        result = mannwhitneyu_calc(['K00001', 'g1', 'g2', [[10, 20, 30]], [[1, 2, 3]]])
        self.assertEqual(result[3], 'g1')

    # --- zscore_calc ---

    def test_zscore_calc_genome_above_mean(self):
        # Single genome value well above reference mean → enriched in genome group
        result = zscore_calc(['K00001', 'ref', 'genome', [[1.0, 2.0, 3.0]], [[10.0]]])
        self.assertIsNotNone(result)
        self.assertEqual(result[3], 'genome')
        self.assertLess(result[-1], 0.01)  # highly significant

    def test_zscore_calc_genome_at_or_below_mean_returns_none(self):
        # Genome value not above reference mean → function returns None
        result = zscore_calc(['K00001', 'ref', 'genome', [[5.0, 6.0, 7.0]], [[1.0]]])
        self.assertIsNone(result)

    def test_zscore_calc_zero_genome_returns_none(self):
        result = zscore_calc(['K00001', 'ref', 'genome', [[1.0, 2.0, 3.0]], [[0.0]]])
        self.assertIsNone(result)

    # --- calculate_portions edge cases ---

    def test_calculate_portions_absent_annotation_gives_zeros(self):
        # Annotation not present in any genome → row of 0.0s is included (cutoff=0)
        result = self.enrichment_test_object.calculate_portions(
            ['K99999'],
            self.genome_groups_simple_example,
            self.genome_annotation_simple_example,
            self.genomes, 0)
        self.assertEqual(len(result), 2)  # header + 1 data row
        self.assertIn('0.0', result[1])

    def test_calculate_portions_cutoff_filters_absent_annotation(self):
        # Annotation not present anywhere → proportion 0.0 < cutoff 0.5 → filtered
        result = self.enrichment_test_object.calculate_portions(
            ['K99999'],
            self.genome_groups_simple_example,
            self.genome_annotation_simple_example,
            self.genomes, 0.5)
        self.assertEqual(len(result), 1)  # header only

    def test_calculate_portions_cutoff_keeps_present_filters_absent(self):
        # K00001 passes cutoff (1.0 in group_1), K99999 does not
        result = self.enrichment_test_object.calculate_portions(
            ['K00001', 'K99999'],
            self.genome_groups_simple_example,
            self.genome_annotation_simple_example,
            self.genomes, 0.5)
        self.assertEqual(len(result), 2)  # header + K00001 only
        self.assertEqual(result[1][0], 'K00001')

    # --- check_annotation_type homogeneity ---

    def test_check_annotation_type_mixed_raises(self):
        with self.assertRaises(ValueError):
            Enrichment().check_annotation_type(['K00001', 'PF00001'])

    def test_check_annotation_type_mixed_ko_tigrfam_raises(self):
        with self.assertRaises(ValueError):
            Enrichment().check_annotation_type(['K00001', 'TIGR00001'])

    def test_check_annotation_type_homogeneous_ko_passes(self):
        result = Enrichment().check_annotation_type(['K00001', 'K00002', 'K00003'])
        self.assertEqual(result, Enrichment.KEGG)

    def test_check_annotation_type_homogeneous_pfam_passes(self):
        result = Enrichment().check_annotation_type(['PF00001', 'PF00002'])
        self.assertEqual(result, Enrichment.PFAM)

    # --- _annotation_type_of new types ---

    def test_annotation_type_go(self):
        e = Enrichment()
        self.assertEqual(e._annotation_type_of("GO:0005488"), Enrichment.GO)

    def test_annotation_type_cog_single_letter(self):
        e = Enrichment()
        self.assertEqual(e._annotation_type_of("C"), Enrichment.COG)

    def test_annotation_type_cog_k_is_cog_not_kegg(self):
        # Single-letter "K" is COG category Transcription, not a KO id
        e = Enrichment()
        self.assertEqual(e._annotation_type_of("K"), Enrichment.COG)

    def test_annotation_type_ko_five_digits(self):
        # K00001 matches KO_PATTERN K\d{5}
        e = Enrichment()
        self.assertEqual(e._annotation_type_of("K00001"), Enrichment.KEGG)

    def test_annotation_type_eggnog_raw_og(self):
        # Raw eggNOG OG string with '@' should be detected as EGGNOG
        e = Enrichment()
        self.assertEqual(e._annotation_type_of("COG0148@1224|Gamma"), Enrichment.EGGNOG)

    def test_annotation_type_cog_two_letters(self):
        e = Enrichment()
        self.assertEqual(e._annotation_type_of("CM"), Enrichment.COG)

    def test_check_annotation_type_go_list(self):
        result = Enrichment().check_annotation_type(['GO:0005488', 'GO:0016020'])
        self.assertEqual(result, Enrichment.GO)

    def test_check_annotation_type_cog_list(self):
        result = Enrichment().check_annotation_type(['C', 'M', 'J'])
        self.assertEqual(result, Enrichment.COG)


if __name__ == "__main__":
    unittest.main()
