#!/usr/bin/env python3
import os
import sqlite3
import pickle
import logging
from enrichm.data import Data

###############################################################################

class Databases:
    """Access the enrichm reference database.

    Supports two on-disk formats:
      - New (SQLite): $ENRICHM_DB/enrichm.db  — created by 'enrichm data --create'
      - Legacy (pickle): $ENRICHM_DB/<version>/  — the old tarball-based layout

    The new format is preferred. The legacy format is used automatically as a
    fallback so that existing database installations continue to work until
    users migrate via 'enrichm data --create'.
    """

    def __init__(self):
        self.db_path = os.path.join(Data.DATABASE_DIR, Data.DB_FILENAME)
        self.hmm_dir = os.path.join(Data.DATABASE_DIR, 'hmm')
        self.KO_HMM_CUTOFFS = os.path.join(Data.DATABASE_DIR, 'ko_cutoffs.tsv')

        self._db_available = os.path.isfile(self.db_path)
        self._legacy = False

        # Fall back to old pickle-based database if SQLite DB doesn't exist
        if not self._db_available:
            version_file = os.path.join(Data.DATABASE_DIR, 'VERSION')
            if os.path.isfile(version_file):
                with open(version_file) as fh:
                    db_version = fh.readline().strip().replace('.tar.gz', '')
                cur_dir = os.path.join(Data.DATABASE_DIR, db_version)
                pickle_version_file = os.path.join(cur_dir, 'VERSION')
                if os.path.isfile(pickle_version_file):
                    with open(pickle_version_file) as fh:
                        self._pickle_version = fh.readline().strip()
                    self._legacy_dir = cur_dir
                    self._legacy = True
                    logging.debug(
                        'Using legacy pickle database at %s '
                        '(run "enrichm data --create" to migrate)',
                        cur_dir,
                    )
                    # Legacy HMM paths (old layout)
                    legacy_ref = os.path.join(cur_dir, 'databases')
                    self.hmm_dir = legacy_ref
                    self.KO_HMM_CUTOFFS = os.path.join(cur_dir, 'ko_cutoffs.tsv')

        # HMM database paths
        self.KO_HMM_DB  = os.path.join(self.hmm_dir, 'ko.hmm')
        self.PFAM_DB    = os.path.join(self.hmm_dir, 'pfam.hmm')
        self.TIGRFAM_DB = os.path.join(self.hmm_dir, 'tigrfam.hmm')
        self.CAZY_DB    = os.path.join(self.hmm_dir, 'cazy.hmm')
        # Diamond databases (legacy layout; not created by new 'enrichm data')
        self.KO_DB = os.path.join(self.hmm_dir, 'uniref100.KO.dmnd')
        self.EC_DB = os.path.join(self.hmm_dir, 'uniref100.EC.dmnd')

        self.signature_modules = set([
            'M00611', 'M00612', 'M00613', 'M00614',
            'M00617', 'M00618', 'M00615', 'M00616',
            'M00363', 'M00542', 'M00574', 'M00575',
            'M00564', 'M00660', 'M00664', 'M00625',
            'M00627', 'M00745', 'M00651', 'M00652',
            'M00704', 'M00725', 'M00726', 'M00730',
            'M00744', 'M00718', 'M00639', 'M00641',
            'M00642', 'M00643', 'M00769', 'M00649',
            'M00696', 'M00697', 'M00698', 'M00700',
            'M00702', 'M00714', 'M00705', 'M00746',
        ])

    # -------------------------------------------------------------------------
    # Internal helpers
    # -------------------------------------------------------------------------

    def _require_db(self):
        if not self._db_available and not self._legacy:
            raise Exception(
                f"\nNo enrichm database found at {Data.DATABASE_DIR}.\n"
                f"Have you:\n"
                f"- Installed the EnrichM database using 'enrichm data --create'?\n"
                f"- Set the ENRICHM_DB environment variable? "
                f"(Currently looking here: {Data.DATABASE_DIR})"
            )

    def _query_dict(self, sql, params=()):
        """Return a dict from a two-column SELECT (key, value)."""
        self._require_db()
        with sqlite3.connect(self.db_path) as conn:
            return dict(conn.execute(sql, params).fetchall())

    def _query_list(self, sql, params=()):
        """Return a flat list from a single-column SELECT."""
        self._require_db()
        with sqlite3.connect(self.db_path) as conn:
            return [row[0] for row in conn.execute(sql, params).fetchall()]

    def _load_pickle(self, name):
        """Load a pickle file from the legacy database directory."""
        path = '.'.join([
            os.path.join(self._legacy_dir, name),
            self._pickle_version,
            'pickle',
        ])
        with open(path, 'rb') as fh:
            return pickle.load(fh)

    # -------------------------------------------------------------------------
    # Description lookups
    # -------------------------------------------------------------------------

    def m2def(self):
        logging.debug('Loading module definitions')
        if self._legacy:
            return self._load_pickle('module_to_definition')
        return self._query_dict('SELECT module_id, definition FROM module_definitions')

    def m(self):
        logging.debug('Loading module descriptions')
        if self._legacy:
            return self._load_pickle('module_descriptions')
        return self._query_dict('SELECT module_id, description FROM module_descriptions')

    def k(self):
        logging.debug('Loading KO descriptions')
        if self._legacy:
            return self._load_pickle('ko_descriptions')
        return self._query_dict('SELECT ko_id, description FROM ko_descriptions')

    def pfam2clan(self):
        logging.debug('Loading Pfam clan membership')
        if self._legacy:
            return self._load_pickle('pfam_to_clan')
        return self._query_dict('SELECT pfam_id, clan_id FROM pfam_clans')

    def pfam2description(self):
        logging.debug('Loading Pfam descriptions')
        if self._legacy:
            return self._load_pickle('pfam_to_description')
        return self._query_dict('SELECT pfam_id, description FROM pfam_descriptions')

    def ec2description(self):
        logging.debug('Loading EC descriptions')
        if self._legacy:
            return self._load_pickle('ec_to_description')
        return self._query_dict('SELECT ec_id, description FROM ec_descriptions')

    def tigrfamdescription(self):
        logging.debug('Loading TIGRFAM descriptions')
        if self._legacy:
            return self._load_pickle('tigrfam_descriptions')
        return self._query_dict(
            'SELECT tigrfam_id, description FROM tigrfam_descriptions')

    # -------------------------------------------------------------------------
    # ID lists for MatrixGenerator (replaces IDS_DIR flat files)
    # -------------------------------------------------------------------------

    def get_all_ko_ids(self):
        if self._legacy:
            ids_dir = os.path.join(self._legacy_dir, 'ids')
            return [x.strip() for x in open(os.path.join(ids_dir, 'KO_IDS.txt'))]
        return self._query_list('SELECT ko_id FROM ko_descriptions ORDER BY ko_id')

    def get_all_ec_ids(self):
        if self._legacy:
            ids_dir = os.path.join(self._legacy_dir, 'ids')
            return [x.strip() for x in open(os.path.join(ids_dir, 'EC_IDS.txt'))]
        return self._query_list('SELECT ec_id FROM ec_descriptions ORDER BY ec_id')

    def get_all_pfam_ids(self):
        if self._legacy:
            ids_dir = os.path.join(self._legacy_dir, 'ids')
            return [x.strip() for x in open(os.path.join(ids_dir, 'PFAM_IDS.txt'))]
        return self._query_list(
            'SELECT pfam_id FROM pfam_descriptions ORDER BY pfam_id')

    def get_all_tigrfam_ids(self):
        if self._legacy:
            ids_dir = os.path.join(self._legacy_dir, 'ids')
            return [x.strip() for x in open(os.path.join(ids_dir, 'TIGRFAM_IDS.txt'))]
        return self._query_list(
            'SELECT tigrfam_id FROM tigrfam_descriptions ORDER BY tigrfam_id')

    def get_all_cazy_ids(self):
        if self._legacy:
            ids_dir = os.path.join(self._legacy_dir, 'ids')
            cazy_path = os.path.join(ids_dir, 'CAZY_IDS.txt')
            if os.path.isfile(cazy_path):
                return [x.strip() for x in open(cazy_path)]
        names = []
        if os.path.isfile(self.CAZY_DB):
            with open(self.CAZY_DB) as fh:
                for line in fh:
                    if line.startswith('NAME'):
                        names.append(line.split(None, 1)[1].strip())
        return sorted(names)

    # -------------------------------------------------------------------------
    # KOfam cutoff file
    # -------------------------------------------------------------------------

    def parse_ko_cutoffs(self):
        cut_ko = {}
        with open(self.KO_HMM_CUTOFFS) as fh:
            fh.readline()  # header
            for line in fh:
                sline = line.strip().split('\t')
                if len(sline) < 3:
                    continue
                if sline[1] == '-':
                    cut_ko[sline[0]] = [0.0, 'NA']
                else:
                    cut_ko[sline[0]] = [float(sline[1]), sline[2]]
        return cut_ko

    # -------------------------------------------------------------------------
    # Database metadata
    # -------------------------------------------------------------------------

    def download_date(self):
        """Return the ISO-8601 download date, or None for legacy installs."""
        if self._legacy or not self._db_available:
            return None
        with sqlite3.connect(self.db_path) as conn:
            row = conn.execute(
                "SELECT value FROM metadata WHERE key='download_date'"
            ).fetchone()
        return row[0] if row else None
