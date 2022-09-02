#!/usr/bin/env python3
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################
# Imports
import os
import urllib
import shutil
import subprocess
import logging
from pathlib import Path
# Local
###############################################################################

class Data:
    '''
    Utilities for archiving, downloading and updating databases.
    '''
    db_var = "ENRICHM_DB"
    if db_var in os.environ:
        DATABASE_DIR 	= os.environ[db_var]
    else:
        DATA_PATH 		= str(Path.home())
        DATABASE_DIR	= os.path.join(DATA_PATH, 'databases')
    VERSION 		= 'VERSION'
    ARCHIVE_SUFFIX 	= '.tar.gz'

    def __init__(self):
        self.ftp = 'https://data.ace.uq.edu.au/public/enrichm/'

    def _archive_db(self, old_db_file):
        '''
        Archive an old database file

        Parameters
        ----------
        old_db_file	- String. File name of old database file to archive
        '''
        from enrichm.databases import Databases

        if not os.path.isdir(Databases.OLD_DATABASE_PATH):
            logging.info('Creating directory to store databases: %s' % (Databases.OLD_DATABASE_PATH))
            os.makedirs(Databases.OLD_DATABASE_PATH)

        old_db_path_archive \
            = os.path.join(Databases.OLD_DATABASE_PATH, old_db_file + self.ARCHIVE_SUFFIX)
        old_db_path \
            = os.path.join(self.DATABASE_DIR, old_db_file)

        logging.info('Compressing old database')
        cmd = "tar -cvzf %s %s > /dev/null" % (old_db_path_archive, old_db_path)
        subprocess.call(cmd, shell=True)

        logging.info('Cleaning up')
        shutil.rmtree(old_db_path)

    def _download_db(self, new_db_file):
        '''
        Download and decompress a new database file

        Parameters
        ----------
        new_db_file	- String. File name of new database to download and decompress.
        '''

        new_db_path_archive \
            = os.path.join(self.DATABASE_DIR, new_db_file)
        logging.info('Downloading new database: %s' % new_db_file)
        cmd = 'wget -c -q %s -O %s' % (self.ftp + new_db_file, new_db_path_archive)
        subprocess.call(cmd, shell = True)
        cmd = 'wget -c -q %s -O %s' % (self.ftp + self.VERSION, os.path.join(self.DATABASE_DIR, self.VERSION))
        subprocess.call(cmd, shell = True)

        logging.info('Decompressing new database')
        cmd = 'tar -xvzf %s -C %s > /dev/null' % (new_db_path_archive, self.DATABASE_DIR)
        subprocess.call(cmd, shell = True)

        logging.info('Cleaning up')
        os.remove(new_db_path_archive)

    def do(self, uninstall, dry):
        '''
        Check database versions, if they're out of date, archive the old and download the new.
        '''

        if uninstall:

            for file in os.listdir(self.DATABASE_DIR):
                file_path = os.path.join(self.DATABASE_DIR, file)

                if os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                else:
                    os.remove(file_path)

            os.rmdir(self.DATABASE_DIR)

        else:

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
                    raise Exception(
                        "Unable to locate enrichM database! Please specify its location by creating a local BASH variable called ENRICHM_DIR: export ENRICHM_DB=/path/to/database/")

                if version_local!=version_remote:
                    logging.info('New database found. Archiving old database.')
                    self._archive_db(version_local.replace(self.ARCHIVE_SUFFIX,''))
                    self._download_db(version_remote)
                else:
                    logging.info('Database is up to date!')

            else:
                logging.info('Creating file to store databases.')
                os.makedirs(self.DATABASE_DIR)
                self._download_db(version_remote)

