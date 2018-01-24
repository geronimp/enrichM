#!/usr/bin/env python
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
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
import logging
import pickle
import inspect

###############################################################################


class Databases:
	"""Databases yo"""

	DATA_PATH   			= os.path.join(os.path.dirname(inspect.stack()[-1][1]), '..', 'share', 'enrichm')
	DATABASE_DIR			= os.path.join(DATA_PATH, 'databases')
	if os.path.isfile(os.path.join(DATABASE_DIR, 'VERSION')):
		DB_VERSION				= open(os.path.join(DATABASE_DIR, 'VERSION')).readline().strip().replace('.tar.gz','')
		CUR_DATABASE_DIR		= os.path.join(DATABASE_DIR, DB_VERSION)
		PICKLE_VERSION			= open(os.path.join(CUR_DATABASE_DIR, 'VERSION')).readline().strip()
		OLD_DATABASE_PATH		= os.path.join(DATA_PATH, 'databases', 'old')
		IDS_DIR					= os.path.join(CUR_DATABASE_DIR, 'ids')
		REF_DIR					= os.path.join(CUR_DATABASE_DIR, 'databases')
		PICKLE					= 'pickle'	
		HMM_SUFFIX 				= '.hmm'
		DMND_SUFFIX				= '.dmnd'
		KO_DB_NAME				= 'uniref100'
		PFAM_DB_NAME			= 'pfam'
		TIGRFAM_DB_NAME			= 'tigrfam'

		M2DEF       			= os.path.join(CUR_DATABASE_DIR, 'module_to_definition')
		M           			= os.path.join(CUR_DATABASE_DIR, 'module_descriptions')
		COMPOUND_DESC_PICKLE	= os.path.join(CUR_DATABASE_DIR, 'br08001')    
		R2RPAIR 				= os.path.join(CUR_DATABASE_DIR, 'reaction_to_rpair')
		R2K     				= os.path.join(CUR_DATABASE_DIR, 'reaction_to_orthology')
		R2C 					= os.path.join(CUR_DATABASE_DIR, 'reaction_to_compound')
		R2M 					= os.path.join(CUR_DATABASE_DIR, 'reaction_to_module')
		M2R 					= os.path.join(CUR_DATABASE_DIR, 'module_to_reaction')
		R2P 					= os.path.join(CUR_DATABASE_DIR, 'reaction_to_pathway')
		P2R 					= os.path.join(CUR_DATABASE_DIR, 'pathway_to_reaction')
		C2R 					= os.path.join(CUR_DATABASE_DIR, 'compound_to_reaction')
		C   					= os.path.join(CUR_DATABASE_DIR, 'compound_descriptions')    
		R   					= os.path.join(CUR_DATABASE_DIR, 'reaction_descriptions')
		P   					= os.path.join(CUR_DATABASE_DIR, 'pathway_descriptions')

		PFAM2CLAN				= os.path.join(CUR_DATABASE_DIR, 'pfam_to_clan')
		CLAN2NAME				= os.path.join(CUR_DATABASE_DIR, 'clan_to_name')
		PFAM2NAME				= os.path.join(CUR_DATABASE_DIR, 'pfam_to_name')
		PFAM2DESCRIPTION		= os.path.join(CUR_DATABASE_DIR, 'pfam_to_description')
		CLAN2PFAM				= os.path.join(CUR_DATABASE_DIR, 'clan_to_pfam')

	def __init__(self):

		self.signature_modules = set(['M00611', 'M00612', 'M00613', 'M00614',
		'M00617', 'M00618', 'M00615', 'M00616', 'M00363', 'M00542', 'M00574',
		'M00575', 'M00564', 'M00660', 'M00664', 'M00625', 'M00627', 'M00745',
		'M00651', 'M00652', 'M00704', 'M00725', 'M00726', 'M00730', 'M00744',
		'M00718', 'M00639', 'M00641', 'M00642', 'M00643', 'M00769', 'M00649',
		'M00696', 'M00697', 'M00698', 'M00700', 'M00702', 'M00714', 'M00705',
		'M00746'])

		logging.info("Loading databases")
		logging.debug("Loading module definitions")
		self.m2def = pickle.load(open('.'.join([self.M2DEF,
		                                         self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading module descriptions")
		self.m = pickle.load(open('.'.join([self.M,
		                                    self.PICKLE_VERSION, self.PICKLE])))	
		logging.debug("Loading reaction to pathway information")
		self.r2p = pickle.load(open('.'.join([self.R2P, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading pathway to reaction information")
		self.p2r = pickle.load(open('.'.join([self.P2R, 
		                                         self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading reaction to orthology information")
		self.r2k = pickle.load(open('.'.join([self.R2K, self.PICKLE_VERSION, 
		                                      self.PICKLE])))
		logging.debug("Loading reaction to module information")
		self.r2m = pickle.load(open('.'.join([self.R2M, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading reaction to module information")
		self.m2r = pickle.load(open('.'.join([self.M2R, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading reaction to compound information")
		self.r2c = pickle.load(open('.'.join([self.R2C, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading compound to reaction information")
		self.c2r = pickle.load(open('.'.join([self.C2R,
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading compound descriptions")
		self.c   = pickle.load(open('.'.join([self.C, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading pathway descriptions")
		self.p   = pickle.load(open('.'.join([self.P, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading reaction descriptions")
		self.r   = pickle.load(open('.'.join([self.R, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading compound classifications")
		self.compound_desc_dict \
		         = pickle.load(open('.'.join([self.COMPOUND_DESC_PICKLE, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading pfam to clan information")
		self.pfam2clan \
		         = pickle.load(open('.'.join([self.PFAM2CLAN, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading clan descriptions")
		self.clan2name \
		         = pickle.load(open('.'.join([self.CLAN2NAME, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading pfam descriptions")
		self.pfam2name \
		         = pickle.load(open('.'.join([self.PFAM2NAME, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading pfam")
		self.pfam2description \
		         = pickle.load(open('.'.join([self.PFAM2DESCRIPTION, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.debug("Loading pfam")
		self.clan2pfam \
		         = pickle.load(open('.'.join([self.CLAN2PFAM, 
		                                      self.PICKLE_VERSION, self.PICKLE])))
		logging.info("Loading reference db paths")		
		self.KO_DB 			= os.path.join(self.REF_DIR, self.KO_DB_NAME + self.DMND_SUFFIX)
		self.PFAM_DB 		= os.path.join(self.REF_DIR, self.PFAM_DB_NAME + self.HMM_SUFFIX)
		self.TIGRFAM_DB 	= os.path.join(self.REF_DIR, self.TIGRFAM_DB_NAME + self.HMM_SUFFIX)
		self.PFAM_CLAN_DB 	= os.path.join(self.IDS_DIR, 'PFAM_CLANS.txt')
