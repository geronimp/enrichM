#!/usr/bin/env python
# Imports
import unittest
import os.path
import sqlite3
import sys
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')]+sys.path

from enrichm.data import Data

###############################################################################

HMM_FILE = '''HMMER3/f [3.3 | Nov 2019]
NAME  1-cysPrx_C
ACC   PF10417.15
DESC  C-terminal domain of 1-Cys peroxiredoxin
LENG  40
CL    CL0172
//
HMMER3/f [3.3 | Nov 2019]
NAME  Thioredoxin
ACC   PF00085.24
DESC  Thioredoxin
LENG  103
CL    CL0172
//
HMMER3/f [3.3 | Nov 2019]
NAME  Clanless
ACC   PF99999.1
DESC  A family that belongs to no clan
LENG  50
//
HMMER3/f [3.3 | Nov 2019]
NAME  NoDescription
ACC   PF88888.2
LENG  50
//
'''

class Tests(unittest.TestCase):

    def _hmm_file(self, tmp_dir):
        hmm_path = os.path.join(tmp_dir, 'test.hmm')

        with open(hmm_path, 'w', encoding='utf-8') as handle:
            handle.write(HMM_FILE)

        return hmm_path

    def test_hmm_metadata_is_keyed_on_accession(self):
        # Annotations are stored under the accession hmmsearch reports, so the
        # descriptions (which supply the frequency table rows) must match it.
        with tempfile.TemporaryDirectory() as tmp_dir:
            descriptions, _ = Data()._parse_hmm_file(self._hmm_file(tmp_dir))

            self.assertEqual(descriptions['PF10417.15'],
                             'C-terminal domain of 1-Cys peroxiredoxin')
            self.assertIn('PF00085.24', descriptions)
            self.assertNotIn('1-cysPrx_C', descriptions)

    def test_hmm_metadata_falls_back_to_model_name(self):
        # Every model must be represented, even without a DESC line, or it goes
        # missing from the frequency table entirely.
        with tempfile.TemporaryDirectory() as tmp_dir:
            descriptions, _ = Data()._parse_hmm_file(self._hmm_file(tmp_dir))

            self.assertEqual(descriptions['PF88888.2'], 'NoDescription')

    def test_clans_are_keyed_without_accession_version(self):
        # Sequence.same_clan strips the version before looking clans up.
        with tempfile.TemporaryDirectory() as tmp_dir:
            _, clans = Data()._parse_hmm_file(self._hmm_file(tmp_dir))

            self.assertEqual(clans, {'PF10417': 'CL0172', 'PF00085': 'CL0172'})

    def test_name_keyed_metadata_is_refreshed(self):
        # Databases built before this fix are keyed on model names, which match
        # nothing at annotation time and must be replaced rather than skipped.
        with tempfile.TemporaryDirectory() as tmp_dir:
            data = Data()
            conn = sqlite3.connect(os.path.join(tmp_dir, 'enrichm.db'))
            data._ensure_schema(conn)
            conn.execute('INSERT INTO pfam_descriptions VALUES (?, ?)',
                         ('1-cysPrx_C', 'C-terminal domain of 1-Cys peroxiredoxin'))
            conn.commit()

            data._populate_pfam_metadata(conn, self._hmm_file(tmp_dir))

            ids = {row[0] for row in conn.execute('SELECT pfam_id FROM pfam_descriptions')}
            self.assertNotIn('1-cysPrx_C', ids)
            self.assertIn('PF10417.15', ids)
            self.assertEqual(
                conn.execute('SELECT count(*) FROM pfam_clans').fetchone()[0], 2)
            conn.close()

    def test_pfam_clans_file_is_parsed(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            clans_path = os.path.join(tmp_dir, Data.PFAM_CLANS_FILENAME)

            with open(clans_path, 'w', encoding='utf-8') as handle:
                handle.write('PF10417\tCL0172\tRedoxin\t1-cysPrx_C\tC-terminal domain\n')
                handle.write('PF00085\tCL0172\tRedoxin\tThioredoxin\tThioredoxin\n')
                handle.write('PF99999\t\t\tClanless\tBelongs to no clan\n')
                handle.write('PF99998\t\\N\t\\N\tAlsoClanless\tBelongs to no clan\n')

            clans = Data()._parse_pfam_clans_file(clans_path)

            self.assertEqual(clans, {'PF10417': 'CL0172', 'PF00085': 'CL0172'})

    def test_clans_are_taken_from_tsv_when_hmm_has_none(self):
        # Current Pfam releases carry no CL lines in Pfam-A.hmm, so clan
        # membership has to come from Pfam-A.clans.tsv instead.
        with tempfile.TemporaryDirectory() as tmp_dir:
            hmm_path = os.path.join(tmp_dir, 'test.hmm')

            with open(hmm_path, 'w', encoding='utf-8') as handle:
                handle.write(HMM_FILE.replace('CL    CL0172\n', ''))

            with open(os.path.join(tmp_dir, Data.PFAM_CLANS_FILENAME), 'w',
                      encoding='utf-8') as handle:
                handle.write('PF10417\tCL0172\tRedoxin\t1-cysPrx_C\tC-terminal domain\n')

            data = Data()
            conn = sqlite3.connect(os.path.join(tmp_dir, 'enrichm.db'))
            data._ensure_schema(conn)
            data._populate_pfam_metadata(conn, hmm_path)

            self.assertEqual(conn.execute('SELECT pfam_id, clan_id FROM pfam_clans').fetchall(),
                             [('PF10417', 'CL0172')])
            conn.close()

    def test_accession_keyed_metadata_is_left_alone(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            data = Data()
            conn = sqlite3.connect(os.path.join(tmp_dir, 'enrichm.db'))
            data._ensure_schema(conn)
            conn.execute('INSERT INTO pfam_descriptions VALUES (?, ?)',
                         ('PF00001.1', 'only entry'))
            conn.commit()

            data._populate_pfam_metadata(conn, self._hmm_file(tmp_dir))

            ids = {row[0] for row in conn.execute('SELECT pfam_id FROM pfam_descriptions')}
            self.assertEqual(ids, {'PF00001.1'})
            conn.close()

if __name__ == "__main__":
    unittest.main()
