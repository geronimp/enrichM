#!/usr/bin/env python3
import os
import re
import ssl
import sqlite3
import urllib.request
import gzip
import tarfile
import shutil
import logging
import time
from datetime import datetime, timezone
from pathlib import Path

from enrichm.toolbox import run_command

###############################################################################

class Data:
    db_var = "ENRICHM_DB"
    DB_FILENAME = 'enrichm.db'
    KEGG_BASE = 'https://rest.kegg.jp'

    # Canonical upstream sources for HMM annotation databases.
    # kofam entry: [profiles_url, ko_list_url]
    HMM_DATABASES = {
        'pfam':     ['https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'],
        'kofam':    ['https://www.genome.jp/ftp/db/kofam/profiles.tar.gz',
                     'https://www.genome.jp/ftp/db/kofam/ko_list.gz'],
        'hmm_pgap': ['https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz'],
        'dbcan3':   ['https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2_9-13-2025/dbCAN.hmm'],
    }

    # Clan membership is not carried in the CL lines of current Pfam-A.hmm
    # releases, so it is taken from this file instead.
    PFAM_CLANS_URL = ('https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/'
                      'Pfam-A.clans.tsv.gz')
    PFAM_CLANS_FILENAME = 'Pfam-A.clans.tsv'

    if db_var in os.environ:
        DATABASE_DIR = os.environ[db_var]
    else:
        DATABASE_DIR = os.path.join(str(Path.home()), 'databases')

    def __init__(self, db_path=None):
        if db_path is not None:
            self.DATABASE_DIR = db_path
        self._ssl_ctx = self._make_ssl_context()

    @staticmethod
    def _make_ssl_context():
        """Return an SSL context with a valid CA bundle.

        conda/pixi environments often lack access to the system CA store.
        Prefer certifi's bundle if available, otherwise fall back to the
        default context (which works fine outside isolated environments).
        """
        try:
            import certifi
            return ssl.create_default_context(cafile=certifi.where())
        except ImportError:
            return ssl.create_default_context()

    def _download(self, url, dest):
        """Download url to dest using the shared SSL context."""
        logging.debug('Downloading %s', url)
        with urllib.request.urlopen(url, context=self._ssl_ctx) as resp, \
                open(dest, 'wb') as fh:
            shutil.copyfileobj(resp, fh)

    # -------------------------------------------------------------------------
    # SQLite schema helpers
    # -------------------------------------------------------------------------

    @staticmethod
    def _ensure_schema(conn):
        """Create all tables if they do not already exist."""
        c = conn.cursor()
        c.execute('CREATE TABLE IF NOT EXISTS metadata '
                  '(key TEXT PRIMARY KEY, value TEXT)')
        c.execute('CREATE TABLE IF NOT EXISTS ko_descriptions '
                  '(ko_id TEXT PRIMARY KEY, description TEXT)')
        c.execute('CREATE TABLE IF NOT EXISTS module_descriptions '
                  '(module_id TEXT PRIMARY KEY, description TEXT)')
        c.execute('CREATE TABLE IF NOT EXISTS module_definitions '
                  '(module_id TEXT PRIMARY KEY, definition TEXT)')
        c.execute('CREATE TABLE IF NOT EXISTS ec_descriptions '
                  '(ec_id TEXT PRIMARY KEY, description TEXT)')
        c.execute('CREATE TABLE IF NOT EXISTS pfam_descriptions '
                  '(pfam_id TEXT PRIMARY KEY, description TEXT)')
        c.execute('CREATE TABLE IF NOT EXISTS pfam_clans '
                  '(pfam_id TEXT PRIMARY KEY, clan_id TEXT)')
        c.execute('CREATE TABLE IF NOT EXISTS tigrfam_descriptions '
                  '(tigrfam_id TEXT PRIMARY KEY, description TEXT)')
        conn.commit()

    @staticmethod
    def _table_has_data(conn, table):
        """Return True if table contains at least one row."""
        row = conn.execute(f'SELECT 1 FROM {table} LIMIT 1').fetchone()
        return row is not None

    # -------------------------------------------------------------------------
    # HMM readiness check
    # -------------------------------------------------------------------------

    @staticmethod
    def _hmm_ready(hmm_path):
        """Return True if hmmpress has already been run (*.h3m marker exists)."""
        return os.path.isfile(hmm_path + '.h3m')

    # -------------------------------------------------------------------------
    # KEGG REST API helpers
    # -------------------------------------------------------------------------

    def _kegg_fetch(self, endpoint):
        url = f'{self.KEGG_BASE}/{endpoint}'
        logging.debug('KEGG API: %s', url)
        with urllib.request.urlopen(url, context=self._ssl_ctx) as resp:
            return resp.read().decode('utf-8')

    def _populate_ko_descriptions(self, conn):
        """Fetch KO descriptions and write to DB. Skips if table already has data."""
        if self._table_has_data(conn, 'ko_descriptions'):
            logging.info('  KO descriptions already present — skipping')
            return
        text = self._kegg_fetch('list/ko')
        rows = []
        for line in text.strip().split('\n'):
            if '\t' in line:
                ko_id, desc = line.split('\t', 1)
                rows.append((ko_id.replace('ko:', ''), desc.strip()))
        conn.executemany('INSERT OR IGNORE INTO ko_descriptions VALUES (?, ?)', rows)
        conn.commit()
        logging.info('  Fetched %d KO descriptions', len(rows))

    def _populate_ec_descriptions(self, conn):
        """Fetch EC descriptions and write to DB. Skips if table already has data."""
        if self._table_has_data(conn, 'ec_descriptions'):
            logging.info('  EC descriptions already present — skipping')
            return
        text = self._kegg_fetch('list/ec')
        rows = []
        for line in text.strip().split('\n'):
            if '\t' in line:
                ec_id, desc = line.split('\t', 1)
                rows.append((ec_id.replace('ec:', ''), desc.strip()))
        conn.executemany('INSERT OR IGNORE INTO ec_descriptions VALUES (?, ?)', rows)
        conn.commit()
        logging.info('  Fetched %d EC descriptions', len(rows))

    def _populate_module_info(self, conn):
        """Fetch module descriptions and definitions, writing to DB incrementally.

        Module descriptions are fetched in one call and written atomically.
        Module definitions are fetched in batches of 10 and committed after
        each successful batch so a restart can skip already-fetched modules.
        """
        # --- module descriptions ---
        if self._table_has_data(conn, 'module_descriptions'):
            logging.info('  Module descriptions already present — skipping list fetch')
            # Still need the full ID list to drive definitions below
            module_ids = [row[0] for row in
                          conn.execute('SELECT module_id FROM module_descriptions').fetchall()]
        else:
            text = self._kegg_fetch('list/module')
            module_ids = []
            desc_rows = []
            for line in text.strip().split('\n'):
                if '\t' in line:
                    mid, desc = line.split('\t', 1)
                    mid = mid.replace('md:', '')
                    module_ids.append(mid)
                    desc_rows.append((mid, desc.strip()))
            conn.executemany('INSERT OR IGNORE INTO module_descriptions VALUES (?, ?)',
                             desc_rows)
            conn.commit()
            logging.info('  Fetched %d module descriptions', len(desc_rows))

        # --- module definitions (incremental by batch) ---
        already_done = set(
            row[0] for row in
            conn.execute('SELECT module_id FROM module_definitions').fetchall()
        )
        pending = [mid for mid in module_ids if mid not in already_done]

        if not pending:
            logging.info('  Module definitions already complete — skipping')
            return

        batch_size = 10
        n_total = len(module_ids)
        n_done = len(already_done)
        n_pending = len(pending)
        n_batches = (n_pending + batch_size - 1) // batch_size
        logging.info('  Fetching module definitions: %d/%d remaining (%d batches)',
                     n_pending, n_total, n_batches)

        for i in range(0, n_pending, batch_size):
            batch = pending[i:i + batch_size]
            try:
                text = self._kegg_fetch('get/' + '+'.join(batch))
                definitions = {}
                current_id = None
                current_def_lines = []
                in_def = False
                for line in text.split('\n'):
                    if line.startswith('ENTRY'):
                        if current_id and current_def_lines:
                            definitions[current_id] = ' '.join(current_def_lines)
                        current_id = line.split()[1]
                        current_def_lines = []
                        in_def = False
                    elif line.startswith('DEFINITION'):
                        in_def = True
                        current_def_lines = [line[12:].strip()]
                    elif in_def and line.startswith(' '):
                        current_def_lines.append(line.strip())
                    elif in_def:
                        in_def = False
                if current_id and current_def_lines:
                    definitions[current_id] = ' '.join(current_def_lines)
                conn.executemany('INSERT OR IGNORE INTO module_definitions VALUES (?, ?)',
                                 definitions.items())
                conn.commit()
                n_done += len(definitions)
            except Exception as exc:
                logging.warning('Failed to fetch module batch %s: %s', batch, exc)
            time.sleep(0.34)  # KEGG rate limit: ~3 requests/second

        total_defs = conn.execute(
            'SELECT COUNT(*) FROM module_definitions').fetchone()[0]
        logging.info('  Module definitions in database: %d', total_defs)

    # -------------------------------------------------------------------------
    # HMM file parsing
    # -------------------------------------------------------------------------

    def _parse_hmm_file(self, hmm_path):
        """Parse ACC, NAME, DESC, and CL lines from an HMM file.

        Entries are keyed on the accession (ACC) declared by each model, because
        that is the identifier hmmsearch reports in its domtblout output and so
        the identifier annotations are stored under. Models that declare no
        accession fall back to their name. Descriptions default to the model
        name so that every model in the library is represented.

        Returns (descriptions, clans). Clan keys have any accession version
        suffix stripped, matching how Sequence.same_clan looks them up. Clan
        entries are only present for HMMs that belong to a Pfam clan.
        """
        descriptions = {}
        clans = {}
        name = accession = description = clan = None

        def store():
            if name is None:
                return
            identifier = accession or name
            descriptions[identifier] = description or name
            if clan:
                clans[identifier.split('.')[0]] = clan

        with open(hmm_path) as fh:
            for line in fh:
                line = line.rstrip()

                if line.startswith('NAME'):
                    # Start of a new model, so record the previous one.
                    store()
                    name = line.split(None, 1)[1].strip()
                    accession = description = clan = None
                elif line.startswith('ACC') and name:
                    accession = line.split(None, 1)[1].strip()
                elif line.startswith('DESC') and name:
                    description = line.split(None, 1)[1].strip()
                elif line.startswith('CL') and name:
                    clan = line.split(None, 1)[1].strip()
        store()

        logging.info('  Parsed %d entries (%d with clan) from %s',
                     len(descriptions), len(clans), os.path.basename(hmm_path))
        return descriptions, clans

    @staticmethod
    def _table_keyed_by_accession(conn, table, column):
        """Return True if a table's ids look like HMM accessions.

        Databases built before accession keying was fixed are keyed on model
        names instead, which silently matches nothing at annotation time.
        """
        row = conn.execute(f'SELECT {column} FROM {table} LIMIT 1').fetchone()

        if row is None:
            return False

        return bool(re.match(r'^(PF|NF|TIGR)\d', row[0]))

    def _parse_pfam_clans_file(self, clans_path):
        """Parse clan membership from a Pfam-A.clans.tsv file.

        The file is tab separated as (accession, clan, clan name, family name,
        description), with an empty clan field for families that belong to no
        clan. Accessions are unversioned, matching how Sequence.same_clan looks
        them up.
        """
        clans = {}

        with open(clans_path) as fh:
            for line in fh:
                fields = line.rstrip('\n').split('\t')

                if len(fields) < 2:
                    continue

                pfam_id, clan_id = fields[0].split('.')[0].strip(), fields[1].strip()

                if clan_id and clan_id != '\\N':
                    clans[pfam_id] = clan_id

        logging.info('  Parsed %d clan memberships from %s',
                     len(clans), os.path.basename(clans_path))
        return clans

    def _pfam_clans(self, hmm_dir):
        """Return {pfam_id: clan_id}, downloading Pfam-A.clans.tsv if needed.

        Returns an empty dict if the file cannot be fetched, so that a missing
        network does not fail the whole database build.
        """
        clans_path = os.path.join(hmm_dir, self.PFAM_CLANS_FILENAME)

        if not os.path.isfile(clans_path):
            clans_gz = clans_path + '.gz'
            try:
                logging.info('  Downloading Pfam clan membership')
                self._download(self.PFAM_CLANS_URL, clans_gz)
                with gzip.open(clans_gz, 'rb') as gz_in, open(clans_path, 'wb') as out:
                    shutil.copyfileobj(gz_in, out)
                os.remove(clans_gz)
            except Exception as exc:
                logging.warning('  Could not fetch Pfam clan membership from %s: %s',
                                self.PFAM_CLANS_URL, exc)
                return {}

        return self._parse_pfam_clans_file(clans_path)

    def _populate_pfam_metadata(self, conn, hmm_path):
        """Parse Pfam HMM metadata and write to DB. Skips if already present."""
        if self._table_has_data(conn, 'pfam_descriptions'):
            if self._table_keyed_by_accession(conn, 'pfam_descriptions', 'pfam_id'):
                logging.info('  Pfam metadata already present — skipping')
                return
            logging.info('  Pfam metadata is keyed on model names — refreshing')
            conn.execute('DELETE FROM pfam_descriptions')
            conn.execute('DELETE FROM pfam_clans')

        pfam_desc, pfam_clans = self._parse_hmm_file(hmm_path)

        if not pfam_clans:
            pfam_clans = self._pfam_clans(os.path.dirname(hmm_path))

        if not pfam_clans:
            logging.warning('  No Pfam clan membership available. Overlapping Pfam '
                            'domains from the same clan cannot be resolved without '
                            'it.')

        conn.executemany('INSERT OR REPLACE INTO pfam_descriptions VALUES (?, ?)',
                         pfam_desc.items())
        conn.executemany('INSERT OR REPLACE INTO pfam_clans VALUES (?, ?)',
                         pfam_clans.items())
        conn.commit()

    def _populate_tigrfam_metadata(self, conn, hmm_path):
        """Parse TIGRFAM HMM metadata and write to DB. Skips if already present."""
        if self._table_has_data(conn, 'tigrfam_descriptions'):
            if self._table_keyed_by_accession(conn, 'tigrfam_descriptions', 'tigrfam_id'):
                logging.info('  TIGRFAM metadata already present — skipping')
                return
            logging.info('  TIGRFAM metadata is keyed on model names — refreshing')
            conn.execute('DELETE FROM tigrfam_descriptions')

        tigrfam_desc, _ = self._parse_hmm_file(hmm_path)
        conn.executemany('INSERT OR REPLACE INTO tigrfam_descriptions VALUES (?, ?)',
                         tigrfam_desc.items())
        conn.commit()

    # -------------------------------------------------------------------------
    # HMM database downloads
    # -------------------------------------------------------------------------

    def _download_hmm_databases(self, dbs, hmm_dir):
        """Download, decompress, and hmmpress each HMM database.

        Each database is skipped when its *.h3m marker file already exists,
        making re-runs safe after a partial failure.

        Returns {name: path} for all HMM files (downloaded or pre-existing).
        """
        hmm_paths = {}

        for name, urls in dbs.items():

            if name == 'pfam':
                dest = os.path.join(hmm_dir, 'pfam.hmm')
                if self._hmm_ready(dest):
                    logging.info('Skipping pfam — already downloaded and pressed')
                else:
                    logging.info('Downloading pfam')
                    dest_gz = os.path.join(hmm_dir, 'pfam.hmm.gz')
                    self._download(urls[0], dest_gz)
                    with gzip.open(dest_gz, 'rb') as gz_in, open(dest, 'wb') as out:
                        shutil.copyfileobj(gz_in, out)
                    os.remove(dest_gz)
                    run_command(f'hmmpress {dest}')
                hmm_paths['pfam'] = dest

            elif name == 'dbcan3':
                dest = os.path.join(hmm_dir, 'cazy.hmm')
                if self._hmm_ready(dest):
                    logging.info('Skipping dbcan3 — already downloaded and pressed')
                else:
                    logging.info('Downloading dbcan3')
                    self._download(urls[0], dest)
                    run_command(f'hmmpress {dest}')
                hmm_paths['cazy'] = dest

            elif name == 'kofam':
                ko_hmm = os.path.join(hmm_dir, 'ko.hmm')
                ko_cutoffs = os.path.join(self.DATABASE_DIR, 'ko_cutoffs.tsv')
                if self._hmm_ready(ko_hmm) and os.path.isfile(ko_cutoffs):
                    logging.info('Skipping kofam — already downloaded and pressed')
                else:
                    logging.info('Downloading kofam')
                    if not self._hmm_ready(ko_hmm):
                        profiles_gz = os.path.join(hmm_dir, 'kofam_profiles.tar.gz')
                        self._download(urls[0], profiles_gz)
                        with tarfile.open(profiles_gz) as tar:
                            tar.extractall(hmm_dir)
                        os.remove(profiles_gz)

                        profiles_dir = os.path.join(hmm_dir, 'profiles')
                        with open(ko_hmm, 'wb') as out:
                            for profile in sorted(os.listdir(profiles_dir)):
                                if profile.endswith('.hmm'):
                                    with open(os.path.join(profiles_dir, profile), 'rb') as pf:
                                        shutil.copyfileobj(pf, out)
                        shutil.rmtree(profiles_dir)
                        run_command(f'hmmpress {ko_hmm}')

                    if not os.path.isfile(ko_cutoffs):
                        kolist_gz = os.path.join(self.DATABASE_DIR, 'ko_list.gz')
                        self._download(urls[1], kolist_gz)
                        with gzip.open(kolist_gz, 'rb') as gz_in, open(ko_cutoffs, 'wb') as out:
                            shutil.copyfileobj(gz_in, out)
                        os.remove(kolist_gz)

                hmm_paths['ko'] = ko_hmm

            elif name == 'hmm_pgap':
                tigrfam_hmm = os.path.join(hmm_dir, 'tigrfam.hmm')
                if self._hmm_ready(tigrfam_hmm):
                    logging.info('Skipping hmm_pgap — already downloaded and pressed')
                else:
                    logging.info('Downloading hmm_pgap')
                    dest_tgz = os.path.join(hmm_dir, 'hmm_PGAP.tgz')
                    self._download(urls[0], dest_tgz)
                    with tarfile.open(dest_tgz) as tar:
                        tar.extractall(hmm_dir)
                    os.remove(dest_tgz)

                    with open(tigrfam_hmm, 'wb') as out:
                        for root, _, files in os.walk(hmm_dir):
                            if root == hmm_dir:
                                continue  # only recurse into extracted subdirs
                            for fname in sorted(files):
                                if fname.upper().endswith('.HMM'):
                                    with open(os.path.join(root, fname), 'rb') as pf:
                                        shutil.copyfileobj(pf, out)

                    # Clean up extracted subdirectories
                    for entry in os.listdir(hmm_dir):
                        entry_path = os.path.join(hmm_dir, entry)
                        if os.path.isdir(entry_path):
                            shutil.rmtree(entry_path)

                    run_command(f'hmmpress {tigrfam_hmm}')
                hmm_paths['tigrfam'] = tigrfam_hmm

        return hmm_paths

    # -------------------------------------------------------------------------
    # Public entry point
    # -------------------------------------------------------------------------

    def do(self, uninstall, create):
        logging.info('Database location: %s', self.DATABASE_DIR)

        if uninstall:
            for entry in os.listdir(self.DATABASE_DIR):
                entry_path = os.path.join(self.DATABASE_DIR, entry)
                if os.path.isdir(entry_path):
                    shutil.rmtree(entry_path)
                else:
                    os.remove(entry_path)
            os.rmdir(self.DATABASE_DIR)
            logging.info('Database removed')

        elif create:
            os.makedirs(self.DATABASE_DIR, exist_ok=True)
            hmm_dir = os.path.join(self.DATABASE_DIR, 'hmm')
            os.makedirs(hmm_dir, exist_ok=True)

            db_path = os.path.join(self.DATABASE_DIR, self.DB_FILENAME)
            conn = sqlite3.connect(db_path)
            self._ensure_schema(conn)

            # Clear any previously recorded completion date so a partial
            # re-run does not appear finished until all steps succeed.
            conn.execute("DELETE FROM metadata WHERE key='download_date'")
            conn.commit()

            logging.info('Step 1/4: Downloading HMM databases')
            hmm_paths = self._download_hmm_databases(self.HMM_DATABASES, hmm_dir)

            logging.info('Step 2/4: Fetching KEGG data')
            self._populate_ko_descriptions(conn)
            self._populate_ec_descriptions(conn)
            self._populate_module_info(conn)

            logging.info('Step 3/4: Parsing HMM metadata')
            if 'pfam' in hmm_paths:
                self._populate_pfam_metadata(conn, hmm_paths['pfam'])
            if 'tigrfam' in hmm_paths:
                self._populate_tigrfam_metadata(conn, hmm_paths['tigrfam'])

            logging.info('Step 4/4: Finalising database')
            conn.execute('INSERT OR REPLACE INTO metadata VALUES (?, ?)',
                         ('download_date', datetime.now(timezone.utc).isoformat()))
            conn.commit()
            conn.close()
            logging.info('Database created successfully: %s', db_path)
