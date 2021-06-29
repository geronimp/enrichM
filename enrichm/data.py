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
    Utilities for archiving, downloading and updating databases.
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

    def _remove_db(self, old_db_file):
        '''
        Archive an old database file

        Parameters
        ----------
        old_db_file	- String. File name of old database file to archive
        '''
        logging.info('Removing old database.')

        old_db_path = os.path.join(self.DATABASE_DIR, old_db_file)
        shutil.rmtree(old_db_path)

    def _download_db(self, new_db_file):
        '''
        Download and decompress a new database file

        Parameters
        ----------
        new_db_file	- String. File name of new database to download and decompress.
        '''
        
        new_db_path_archive = os.path.join(self.DATABASE_DIR, new_db_file)
        
        logging.info('Downloading new database: %s', new_db_file)
        cmd = f'wget \
                    -q {self.ftp + new_db_file} \
                    -O {new_db_path_archive}'
        run_command(cmd)
        
        cmd = f'wget \
                    -q {self.ftp + self.VERSION} \
                    -O {os.path.join(self.DATABASE_DIR, self.VERSION)}'
        run_command(cmd)

        logging.info('Decompressing new database')
        cmd = 'tar -xvzf %s -C %s > /dev/null' % (new_db_path_archive, self.DATABASE_DIR)
        run_command(cmd)

        logging.info('Cleaning up')
        os.remove(new_db_path_archive)

    def do(self, uninstall, create):
        '''
        Check database versions, if they're out of date, archive the old and download the new.
        '''
        logging.info(f"Database location: {self.DATABASE_DIR}")
        
        if uninstall:

            for file in os.listdir(self.DATABASE_DIR):
                file_path = os.path.join(self.DATABASE_DIR, file)

                if os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                else:
                    os.remove(file_path)

            os.rmdir(self.DATABASE_DIR)

        elif create:
            try:
                version_remote = urllib.request.urlopen(self.ftp + self.VERSION).readline().strip().decode("utf-8")
            except:
                raise Exception(
                    "Unable to find most current EnrichM database VERSION in ftp. Please complain at https://github.com/geronimp/enrichM")

            if os.path.isdir(self.DATABASE_DIR):
                version_local_path = os.path.join(self.DATABASE_DIR, self.VERSION)

                if os.path.isfile(version_local_path):
                    version_local = open(version_local_path).readline().strip()
                else:
                    logging.info(f'EnrichM database not detected in database directory ({self.DATABASE_DIR}). Downloading database.')
                    self._download_db(version_remote)

                if version_local!=version_remote:
                    logging.info('New database found. Archiving old database.')
                    self._remove_db(version_local.replace(self.ARCHIVE_SUFFIX,''))
                    self._download_db(version_remote)
                else:
                    logging.info('Database is up to date!')

            else:
                logging.info('Creating folder to store databases.')
                os.makedirs(self.DATABASE_DIR)
                self._download_db(version_remote)

