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

__author__      = "Joel Boyd"
__copyright__   = "Copyright 2017"
__credits__     = ["Joel Boyd"]
__license__     = "GPL3"
__version__     = "0.0.7"
__maintainer__  = "Joel Boyd"
__email__       = "joel.boyd near uq.net.au"
__status__      = "Development"
 
###############################################################################
# Imports
import os
import logging
import pickle
# Local
from enrichm.data import Data
###############################################################################


class Databases:
	
	if os.path.isfile(os.path.join(Data.DATABASE_DIR, 'VERSION')):
		DB_VERSION			= open(os.path.join(Data.DATABASE_DIR, 'VERSION')).readline().strip().replace('.tar.gz','')
		CUR_DATABASE_DIR	= os.path.join(Data.DATABASE_DIR, DB_VERSION)
		PICKLE_VERSION		= open(os.path.join(CUR_DATABASE_DIR, 'VERSION')).readline().strip()
		OLD_DATABASE_PATH	= os.path.join(Data.DATABASE_DIR, 'old')
		IDS_DIR				= os.path.join(CUR_DATABASE_DIR, 'ids')
		REF_DIR				= os.path.join(CUR_DATABASE_DIR, 'databases')
		GTDB_DIR			= os.path.join(CUR_DATABASE_DIR, 'gtdb')
		KO_HMM_CUTOFFS 		= os.path.join(CUR_DATABASE_DIR, 'ko_cutoffs.tsv')

		PICKLE				= 'pickle'	
		HMM_SUFFIX 			= '.hmm'
		DMND_SUFFIX			= '.dmnd'
		KO_DB_NAME			= 'uniref100.KO'
		EC_DB_NAME			= 'uniref100.EC'
		PFAM_DB_NAME		= 'pfam'
		KO_HMM_DB_NAME		= 'ko'
		TIGRFAM_DB_NAME		= 'tigrfam'
		CAZY_DB_NAME		= 'cazy'
		GTDB_DB_NAME		= 'GTDB_R80_DB'
		
		GTDB_CAZY 			= os.path.join(GTDB_DIR, "gtdb_cazy.tsv")
		GTDB_KO 			= os.path.join(GTDB_DIR, "gtdb_ko.tsv")
		GTDB_PFAM 			= os.path.join(GTDB_DIR, "gtdb_pfam.tsv")
		GTDB_TIGRFAM 		= os.path.join(GTDB_DIR, "gtdb_tigrfam.tsv")
		GTDB_EC 			= os.path.join(GTDB_DIR, "gtdb_ec.tsv")

		TAXONOMY			= os.path.join(CUR_DATABASE_DIR, 'taxonomy_gtdb.tsv')
		M2DEF				= os.path.join(CUR_DATABASE_DIR, 'module_to_definition')
		M					= os.path.join(CUR_DATABASE_DIR, 'module_descriptions')
		COMPOUND_DESC		= os.path.join(CUR_DATABASE_DIR, 'br08001')    
		R2K					= os.path.join(CUR_DATABASE_DIR, 'reaction_to_orthology')
		R2C 				= os.path.join(CUR_DATABASE_DIR, 'reaction_to_compound')
		R2M					= os.path.join(CUR_DATABASE_DIR, 'reaction_to_module')
		M2R					= os.path.join(CUR_DATABASE_DIR, 'module_to_reaction')
		M2C					= os.path.join(CUR_DATABASE_DIR, 'module_to_cpd')
		R2P					= os.path.join(CUR_DATABASE_DIR, 'reaction_to_pathway')
		P2R					= os.path.join(CUR_DATABASE_DIR, 'pathway_to_reaction')
		C2R					= os.path.join(CUR_DATABASE_DIR, 'compound_to_reaction')
		C					= os.path.join(CUR_DATABASE_DIR, 'compound_descriptions')    
		R					= os.path.join(CUR_DATABASE_DIR, 'reaction_descriptions')
		P					= os.path.join(CUR_DATABASE_DIR, 'pathway_descriptions')
		K					= os.path.join(CUR_DATABASE_DIR, 'ko_descriptions')

		PFAM2CLAN			= os.path.join(CUR_DATABASE_DIR, 'pfam_to_clan')
		CLAN2NAME			= os.path.join(CUR_DATABASE_DIR, 'clan_to_name')
		PFAM2NAME			= os.path.join(CUR_DATABASE_DIR, 'pfam_to_name')
		PFAM2DESCRIPTION	= os.path.join(CUR_DATABASE_DIR, 'pfam_to_description')
		EC2DESCRIPTION	= os.path.join(CUR_DATABASE_DIR, 'ec_to_description')
		TIGRFAM2DESCRIPTION= os.path.join(CUR_DATABASE_DIR, 'tigrfam_descriptions')
		CLAN2PFAM			= os.path.join(CUR_DATABASE_DIR, 'clan_to_pfam')

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
		self.m2def = self.load_pickle(self.M2DEF)
		logging.debug("Loading module descriptions")
		self.m = self.load_pickle(self.M)
		logging.debug("Loading reaction to pathway information")
		self.r2p = self.load_pickle(self.R2P)
		logging.debug("Loading pathway to reaction information")
		self.p2r = self.load_pickle(self.P2R)
		logging.debug("Loading reaction to orthology information")
		self.r2k = self.load_pickle(self.R2K)
		logging.debug("Loading reaction to module information")
		self.r2m = self.load_pickle(self.R2M)
		logging.debug("Loading module to reaction information")
		self.m2r = self.load_pickle(self.M2R)
		logging.debug("Loading module to compound information")
		self.m2c = self.load_pickle(self.M2C)
		logging.debug("Loading reaction to compound information")
		self.r2c = self.load_pickle(self.R2C)
		logging.debug("Loading compound to reaction information")
		self.c2r = self.load_pickle(self.C2R)
		logging.debug("Loading compound descriptions")
		self.c = self.load_pickle(self.C)
		logging.debug("Loading pathway descriptions")
		self.p = self.load_pickle(self.P)
		logging.debug("Loading reaction descriptions")
		self.r = self.load_pickle(self.R)
		logging.debug("Loading ko descriptions")
		self.k = self.load_pickle(self.K)
		logging.debug("Loading compound classifications")
		self.compound_desc_dict = self.load_pickle(self.COMPOUND_DESC)
		logging.debug("Loading pfam to clan information")
		self.pfam2clan = self.load_pickle(self.PFAM2CLAN)
		logging.debug("Loading clan descriptions")
		self.clan2name = self.load_pickle(self.CLAN2NAME)
		logging.debug("Loading pfam names")
		self.pfam2name = self.load_pickle(self.PFAM2NAME)
		logging.debug("Loading pfam descriptions")
		self.pfam2description = self.load_pickle(self.PFAM2DESCRIPTION)
		logging.debug("Loading ec descriptions")
		self.ec2description = self.load_pickle(self.EC2DESCRIPTION)
		logging.debug("Loading pfam hierarchy")
		self.clan2pfam = self.load_pickle(self.CLAN2PFAM)
		logging.debug("Loading tigrfam descriptions")
		self.tigrfamdescription = self.load_pickle(self.TIGRFAM2DESCRIPTION)
		logging.info("Loading reference db paths")		


		self.k2r = dict()
		for reaction, kos in self.r2k.items():
			for ko in kos:
				if ko not in self.k2r:
					self.k2r[ko] = list()
				self.k2r[ko].append(reaction)

		self.taxonomy 		= self.parse_taxonomy(self.TAXONOMY)

		self.KO_DB 			= os.path.join(self.REF_DIR, self.KO_DB_NAME + self.DMND_SUFFIX)
		self.EC_DB 			= os.path.join(self.REF_DIR, self.EC_DB_NAME + self.DMND_SUFFIX)
		self.GTDB_DB		= os.path.join(self.REF_DIR, self.GTDB_DB_NAME)

		self.PFAM_DB 		= os.path.join(self.REF_DIR, self.PFAM_DB_NAME + self.HMM_SUFFIX)
		self.KO_HMM_DB 		= os.path.join(self.REF_DIR, self.KO_HMM_DB_NAME + self.HMM_SUFFIX)
		self.TIGRFAM_DB 	= os.path.join(self.REF_DIR, self.TIGRFAM_DB_NAME + self.HMM_SUFFIX)
		self.CAZY_DB 		= os.path.join(self.REF_DIR, self.CAZY_DB_NAME + self.HMM_SUFFIX)
		self.PFAM_CLAN_DB 	= os.path.join(self.IDS_DIR, 'PFAM_CLANS.txt')


	def load_pickle(self, file):

		with open('.'.join([file, self.PICKLE_VERSION, self.PICKLE]), 'rb') as file_io:
			loaded_pickle = pickle.load(file_io)

		return loaded_pickle

	def parse_taxonomy(self, taxonomy_path):
		
		output_taxonomy_dictionary = dict()

		for line in open(taxonomy_path):
			genome, taxonomy_string = line.strip().split('\t')
			output_taxonomy_dictionary[genome] = taxonomy_string.split(';')
			
		return output_taxonomy_dictionary

	def parse_ko_cutoffs(self):

		cut_ko = dict()
		out_io = open(self.KO_HMM_CUTOFFS)
		header = out_io.readline()
		for line in out_io:
			sline = line.strip().split('\t')
			if sline[1]=='-':
				cut_ko[sline[0]] = [0.0, "NA"]
			else:
				cut_ko[sline[0]] = [float(sline[1]), sline[2]]
		return cut_ko
