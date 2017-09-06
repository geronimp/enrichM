#!/usr/bin/env python
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

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

###############################################################################
# Imports
import os
import urllib
import shutil
import subprocess

from databases import Databases

###############################################################################

class Data:
	'''
	Utilities for archiving, downloading and updating databases.
	'''

	VERSION ='VERSION'
	ARCHIVE_SUFFIX = '.tar.gz'
	
	def __init__(self):
		self.ftp = 'https://data.ace.uq.edu.au/public/enrichm/'

	def _archive_db(self, old_db_file):
		'''
		Archive an old database file
		
		Parameters
		----------
		old_db_file	- String. File name of old database file to archive
		'''
		
		if not os.path.isdir(Databases.OLD_DATABASE_PATH):
			logging.info('Creating directory to store databases: %s' % (Databases.OLD_DATABASE_PATH))
			os.makedirs(Databases.OLD_DATABASE_PATH)
	
		old_db_path_archive \
			= os.path.join(Databases.OLD_DATABASE_PATH, old_db_file + self.ARCHIVE_SUFFIX) 
		old_db_path \
			= os.path.join(Databases.DATABASE_PATH, old_db_file) 

		logging.info('Compressing old database')
		cmd = "tar -cvzf %s %s" % (old_db_path_archive, old_db_path)
		subprocess.call(cmd, shell=True)

		logging.info('Cleaning up')
		shutil.rmtree(old_db_path)
		subprocess.call(cmd, shell=True)

	def _download_db(self, new_db_file):
		'''
		Download and decompress a new database file
		
		Parameters
		----------
		new_db_file	- String. File name of new database to download and decompress.
		'''
		new_db_path_archive \
			= os.path.join(Databases.DATABASE_PATH, new_db_file)
		new_db_path \
			= os.path.join(Databases.DATABASE_PATH, new_db_file.replace(self.ARCHIVE_SUFFIX, ''))
		
		logging.info('Downloading new database')
		urllib.urlretrieve(self.ftp + new_db_file, new_db_path_archive)
		
		logging.info('Decompressing new database')
		cmd = 'tar xvzf file.tar.gz -C %s' % (new_db_path)
		subprocess.call(cmd, shell = True)

		logging.info('Cleaning up')
		shutil.remove(new_db_path_archive)
	
	def do(self):
		'''
		Check database versions, if they're out of date, archive the old and download the new.
		'''
		version_remote = urllib.request.urlopen(self.ftp + self.VERSION).readline().strip()
		if os.path.isdir(Databases.DATABASE_DIR):
			version_local  = open(os.path.join(Databases.DATABASE_DIR, Databases.VERSION)).readline().strip()
			if version_local!=version_remote:
				logging.info('New database found. Archiving old database.')
				self._archive_db(version_local)
			else:
				logging.info('Database is up to date!')
		else:
			os.makedirs(Databases.DATABASE_DIR)
		self._download_db(version_remote)
		