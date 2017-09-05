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

from databases import Databases

###############################################################################

class Data:
	VERSION = 'VERSION'
	def __init__(self):
		
		self.ftp = 'https://data.ace.uq.edu.au/public/enrichm/'
		self.d = Databases()
		self.version_file = os.path.join(self.DATABASE_DIR, Databases.VERSION)

	def _rversion(self):
		version=float(open(self.version_file).readline())
		return version

	def create(self):
		logging.debug('Creating Database directory')
		os.makedirs(Databases.DATABASE_DIR)

		logging.debug('Downloading release version')	
		urllib.urlretrieve(self.ftp + Databases.VERSION, self.version_file)

		logging.debug('Downloading UniRef100 database')		
		urllib.urlretrieve(self.ftp + Databases.KO_DB_NAME, self.d.KO_DB)
		logging.debug('Downloading PFAM database')		
		urllib.urlretrieve(self.ftp + Databases.PFAM_DB_NAME, self.d.PFAM_DB)
		logging.debug('Downloading TIGRFAM database')		
		urllib.urlretrieve(self.ftp + Databases.TIGRFAM_DB_NAME, self.d.TIGRFAM_DB)

	def update(self):
		# Check current database isnt being used

		version_remote = float(urllib.request.urlopen(self.ftp + Databases.VERSION).readline().strip())
		version_local = self._rversion()

		if version_local<version_remote:
			if not os.path.isdir(Databases.OLD_DATABASE_PATH):
				os.makedirs(Databases.OLD_DATABASE_PATH)
		
			old_db_directory = os.path.join(Databases.OLD_DATABASE_PATH, str(version_local)) 
			os.makedirs(old_db_directory)

			for old_db_file in os.path.listdir(Databases.DATA_PATH):
				old_db_path 			= os.path.join(Databases.DATA_PATH, old_db_file)
				old_db_path_archived 	= os.path.join(old_db_directory, old_db_file)
				shutil.move(old_db_path, old_db_path_archived)

			logging.debug('Downloading UniRef100 database')		
			urllib.urlretrieve(self.ftp + Databases.KO_DB_NAME, self.d.KO_DB)
			logging.debug('Downloading PFAM database')		
			urllib.urlretrieve(self.ftp + Databases.PFAM_DB_NAME, self.d.PFAM_DB)
			logging.debug('Downloading TIGRFAM database')		
			urllib.urlretrieve(self.ftp + Databases.TIGRFAM_DB_NAME, self.d.TIGRFAM_DB)
		else:
			logging.info('Databases are up to date!')
	
	def do(self, create, update):
		
		if create:
			logging.info('Generating backend databases')
			self.create()
		
		elif update:
			logging.info('Updating backend databases')
			self.update()
