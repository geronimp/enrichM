#!/usr/bin/env python3
# Imports
import os
import urllib.request
import shutil
import logging
from pathlib import Path
from enrichm.toolbox import run_command
# Local
###############################################################################

class Data:
    '''
    Utilities for downloading and updating databases.
    '''
    db_var = "ENRICHM_DB"
    VERSION = 'VERSION'
    ARCHIVE_SUFFIX = '.tar.gz'

    if db_var in os.environ:
        DATABASE_DIR = os.environ[db_var]
    else:
        DATA_PATH = str(Path.home())
        DATABASE_DIR = os.path.join(DATA_PATH, 'databases')

    def __init__(self):
        self.ftp = 'https://data.ace.uq.edu.au/public/enrichm/'

    def _download_db(self, new_db_file):
        '''
        Download and decompress a new database file.

        Parameters
        ----------
        new_db_file - String. Filename of the new database tarball to download.
        '''
        new_db_path_archive = os.path.join(self.DATABASE_DIR, new_db_file)

        logging.info('Downloading new database: %s', new_db_file)
        urllib.request.urlretrieve(self.ftp + new_db_file, new_db_path_archive)
        urllib.request.urlretrieve(
            self.ftp + self.VERSION,
            os.path.join(self.DATABASE_DIR, self.VERSION)
        )

        logging.info('Decompressing new database')
        run_command(f'tar -xvzf {new_db_path_archive} -C {self.DATABASE_DIR} > /dev/null')

        logging.info('Cleaning up')
        os.remove(new_db_path_archive)

    def do(self, uninstall, create):
        '''
        Install, update, or remove the EnrichM database.
        '''
        logging.info(f"Database location: {self.DATABASE_DIR}")

        if uninstall:
            for entry in os.listdir(self.DATABASE_DIR):
                entry_path = os.path.join(self.DATABASE_DIR, entry)
                if os.path.isdir(entry_path):
                    shutil.rmtree(entry_path)
                else:
                    os.remove(entry_path)
            os.rmdir(self.DATABASE_DIR)

        elif create:
            try:
                version_remote = urllib.request.urlopen(
                    self.ftp + self.VERSION
                ).readline().strip().decode('utf-8')
            except Exception:
                raise Exception(
                    "Unable to fetch the current EnrichM database VERSION. "
                    "Please report the issue at https://github.com/geronimp/enrichM"
                )

            if os.path.isdir(self.DATABASE_DIR):
                version_local_path = os.path.join(self.DATABASE_DIR, self.VERSION)

                if os.path.isfile(version_local_path):
                    with open(version_local_path) as fh:
                        version_local = fh.readline().strip()

                    if version_local != version_remote:
                        logging.info('New database version available. Removing old database.')
                        shutil.rmtree(os.path.join(
                            self.DATABASE_DIR,
                            version_local.replace(self.ARCHIVE_SUFFIX, '')
                        ))
                        self._download_db(version_remote)
                    else:
                        logging.info('Database is up to date!')
                else:
                    logging.info(
                        'EnrichM database not detected in %s. Downloading database.',
                        self.DATABASE_DIR
                    )
                    self._download_db(version_remote)
            else:
                logging.info('Creating folder to store databases.')
                os.makedirs(self.DATABASE_DIR)
                self._download_db(version_remote)
