#!/usr/bin/env python
import unittest
import os
import sys
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from enrichm.parser import Parser

# Minimal emapper.annotations fixture content (21 tab-separated columns)
# Comment lines and header start with '#' and are skipped by the parser.
EMAPPER_HEADER = (
    "## emapper version 2.1.0\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
    "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t"
    "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
    "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
)

# genome1: seq1 has KOs K00001,K00002; seq2 has K00001 only
GENOME1_ROWS = (
    "genome1~seq_001\t-\t1e-50\t200\tCOG0148@1224|Gamma,OG6@1|root\tGamma\t"
    "CM\t-\t-\tGO:0005488,GO:0016020\t6.3.4.2\tko:K00001,ko:K00002\t"
    "-\t-\t-\t-\t-\t-\t-\t-\tPF00001,PF00002\n"
    "genome1~seq_002\t-\t1e-30\t150\tCOG0200@1|root\t1|root\t"
    "K\t-\t-\tGO:0005488\t-\tko:K00001\t"
    "-\t-\t-\t-\t-\t-\t-\t-\tPF00001\n"
)

# genome2: seq1 has K00003; no GOs; COG cat J
GENOME2_ROWS = (
    "genome2~seq_001\t-\t1e-10\t100\tCOG0001@1|root\t1|root\t"
    "J\t-\t-\t-\t1.1.1.1\tko:K00003\t"
    "-\t-\t-\t-\t-\t-\t-\t-\tPF00003\n"
)

###############################################################################

class Tests(unittest.TestCase):

    def test_parse_simple_matrix(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("ID\tgenome_1\tgenome_2\n")
            f.write("K00001\t3\t0\n")
            f.write("K00002\t1\t2\n")
            fname = f.name

        try:
            output_dict, colnames, rownames = Parser.parse_simple_matrix(fname)
            self.assertEqual(colnames, ['genome_1', 'genome_2'])
            self.assertEqual(rownames, ['K00001', 'K00002'])
            self.assertEqual(output_dict['genome_1']['K00001'], 3.0)
            self.assertEqual(output_dict['genome_2']['K00001'], 0.0)
            self.assertEqual(output_dict['genome_1']['K00002'], 1.0)
            self.assertEqual(output_dict['genome_2']['K00002'], 2.0)
        finally:
            os.unlink(fname)

    def test_parse_simple_matrix_single_genome(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("ID\tonly_genome\n")
            f.write("K00001\t5\n")
            fname = f.name

        try:
            output_dict, colnames, rownames = Parser.parse_simple_matrix(fname)
            self.assertEqual(colnames, ['only_genome'])
            self.assertEqual(rownames, ['K00001'])
            self.assertEqual(output_dict['only_genome']['K00001'], 5.0)
        finally:
            os.unlink(fname)

    def test_parse_metadata_matrix(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("genome_1\tgroup_a\n")
            f.write("genome_2\tgroup_b\n")
            f.write("genome_3\tgroup_a\n")
            fname = f.name

        try:
            cols_to_rows, nr_values, attribute_dict = Parser.parse_metadata_matrix(fname)
            self.assertEqual(nr_values, {'group_a', 'group_b'})
            self.assertEqual(attribute_dict['group_a'], {'genome_1', 'genome_3'})
            self.assertEqual(attribute_dict['group_b'], {'genome_2'})
            self.assertIn('group_a', cols_to_rows['genome_1'])
            self.assertIn('group_b', cols_to_rows['genome_2'])
        finally:
            os.unlink(fname)

    def test_parse_metadata_matrix_single_group(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("genome_1\tgroup_a\n")
            f.write("genome_2\tgroup_a\n")
            fname = f.name

        try:
            _, nr_values, attribute_dict = Parser.parse_metadata_matrix(fname)
            self.assertEqual(nr_values, {'group_a'})
            self.assertEqual(attribute_dict['group_a'], {'genome_1', 'genome_2'})
        finally:
            os.unlink(fname)

    def test_parse_taxonomy(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("genome_1\td__Bacteria;p__Firmicutes;c__Bacilli\n")
            f.write("genome_2\td__Archaea;p__Euryarchaeota\n")
            fname = f.name

        try:
            result = Parser.parse_taxonomy(fname)
            self.assertEqual(result['genome_1'], ['d__Bacteria', 'p__Firmicutes', 'c__Bacilli'])
            self.assertEqual(result['genome_2'], ['d__Archaea', 'p__Euryarchaeota'])
        finally:
            os.unlink(fname)

    def test_parse_taxonomy_single_rank(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("genome_1\td__Bacteria\n")
            fname = f.name

        try:
            result = Parser.parse_taxonomy(fname)
            self.assertEqual(result['genome_1'], ['d__Bacteria'])
        finally:
            os.unlink(fname)


class TestParseEmapperOutput(unittest.TestCase):

    def _make_emapper_dir(self):
        tmpdir = tempfile.mkdtemp()
        g1 = os.path.join(tmpdir, "genome1.emapper.annotations")
        g2 = os.path.join(tmpdir, "genome2.emapper.annotations")
        with open(g1, 'w') as f:
            f.write(EMAPPER_HEADER + GENOME1_ROWS)
        with open(g2, 'w') as f:
            f.write(EMAPPER_HEADER + GENOME2_ROWS)
        return tmpdir

    def tearDown(self):
        import shutil
        if hasattr(self, '_tmpdir') and os.path.isdir(self._tmpdir):
            shutil.rmtree(self._tmpdir)

    def _parse(self, annotation_type):
        self._tmpdir = self._make_emapper_dir()
        return Parser.parse_emapper_output(self._tmpdir, annotation_type)

    def test_ko_extraction(self):
        genome_ids, tables = self._parse("ko")
        self.assertEqual(sorted(genome_ids), ["genome1", "genome2"])
        idx1 = genome_ids.index("genome1")
        ko_dict = dict(zip(tables[idx1]["annotation"].to_list(),
                           tables[idx1]["count"].to_list()))
        self.assertEqual(ko_dict.get("K00001"), 2)
        self.assertEqual(ko_dict.get("K00002"), 1)
        self.assertNotIn("ko:K00001", ko_dict)

    def test_ko_genome2(self):
        genome_ids, tables = self._parse("ko")
        idx2 = genome_ids.index("genome2")
        ko_dict = dict(zip(tables[idx2]["annotation"].to_list(),
                           tables[idx2]["count"].to_list()))
        self.assertEqual(ko_dict.get("K00003"), 1)

    def test_cog_extraction(self):
        genome_ids, tables = self._parse("cog")
        idx1 = genome_ids.index("genome1")
        cog_dict = dict(zip(tables[idx1]["annotation"].to_list(),
                            tables[idx1]["count"].to_list()))
        # genome1 seq1 has "CM", seq2 has "K" -> C:1, M:1, K:1
        self.assertEqual(cog_dict.get("C"), 1)
        self.assertEqual(cog_dict.get("M"), 1)
        self.assertEqual(cog_dict.get("K"), 1)

    def test_cog_genome2(self):
        genome_ids, tables = self._parse("cog")
        idx2 = genome_ids.index("genome2")
        cog_dict = dict(zip(tables[idx2]["annotation"].to_list(),
                            tables[idx2]["count"].to_list()))
        self.assertEqual(cog_dict.get("J"), 1)

    def test_go_extraction(self):
        genome_ids, tables = self._parse("go")
        idx1 = genome_ids.index("genome1")
        go_dict = dict(zip(tables[idx1]["annotation"].to_list(),
                           tables[idx1]["count"].to_list()))
        # seq1 has GO:0005488,GO:0016020; seq2 has GO:0005488
        self.assertEqual(go_dict.get("GO:0005488"), 2)
        self.assertEqual(go_dict.get("GO:0016020"), 1)

    def test_go_missing_filtered(self):
        genome_ids, tables = self._parse("go")
        idx2 = genome_ids.index("genome2")
        # genome2 seq1 has '-' for GOs -> nothing returned
        self.assertEqual(len(tables[idx2]), 0)

    def test_eggnog_extraction(self):
        genome_ids, tables = self._parse("eggnog")
        idx1 = genome_ids.index("genome1")
        og_dict = dict(zip(tables[idx1]["annotation"].to_list(),
                           tables[idx1]["count"].to_list()))
        # seq1 last OG is "OG6@1|root" -> "OG6"; seq2 "COG0200@1|root" -> "COG0200"
        self.assertIn("OG6", og_dict)
        self.assertIn("COG0200", og_dict)
        self.assertNotIn("@", str(list(og_dict.keys())))

    def test_pfam_extraction(self):
        genome_ids, tables = self._parse("pfam")
        idx1 = genome_ids.index("genome1")
        pf_dict = dict(zip(tables[idx1]["annotation"].to_list(),
                           tables[idx1]["count"].to_list()))
        # seq1 has PF00001,PF00002; seq2 has PF00001
        self.assertEqual(pf_dict.get("PF00001"), 2)
        self.assertEqual(pf_dict.get("PF00002"), 1)

    def test_ec_extraction(self):
        genome_ids, tables = self._parse("ec")
        idx1 = genome_ids.index("genome1")
        ec_dict = dict(zip(tables[idx1]["annotation"].to_list(),
                           tables[idx1]["count"].to_list()))
        self.assertEqual(ec_dict.get("6.3.4.2"), 1)

    def test_ec_genome2(self):
        genome_ids, tables = self._parse("ec")
        idx2 = genome_ids.index("genome2")
        ec_dict = dict(zip(tables[idx2]["annotation"].to_list(),
                           tables[idx2]["count"].to_list()))
        self.assertEqual(ec_dict.get("1.1.1.1"), 1)

    def test_invalid_dir_raises(self):
        with self.assertRaises(Exception):
            Parser.parse_emapper_output("/nonexistent/path", "ko")

    def test_invalid_annotation_type_raises(self):
        tmpdir = tempfile.mkdtemp()
        try:
            with self.assertRaises(ValueError):
                Parser.parse_emapper_output(tmpdir, "invalid_type")
        finally:
            import shutil
            shutil.rmtree(tmpdir)

    def test_non_emapper_files_ignored(self):
        tmpdir = tempfile.mkdtemp()
        try:
            # Write a non-.emapper.annotations file
            with open(os.path.join(tmpdir, "genome1.txt"), 'w') as f:
                f.write("irrelevant\n")
            genome_ids, tables = Parser.parse_emapper_output(tmpdir, "ko")
            self.assertEqual(genome_ids, [])
            self.assertEqual(tables, [])
        finally:
            import shutil
            shutil.rmtree(tmpdir)


if __name__ == "__main__":
    unittest.main()
