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
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
 
###############################################################################
# Imports

import os
import logging
import pickle

###############################################################################


class Databases():
	"""docstring for Databases"""

	DATA_PATH   			= os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data')
	DATABASE_DIR			= os.path.join(DATA_PATH, 'databases')
	VERSION     			= open(os.path.join(DATA_PATH, 'VERSION')).readline().strip()
	PICKLE 					= 'pickle'	

	M2DEF       			= os.path.join(DATA_PATH, 'module_to_definition')
	M           			= os.path.join(DATA_PATH, 'module_descriptions')
	COMPOUND_DESC_PICKLE 	= os.path.join(DATA_PATH, 'br08001')    
	R2RPAIR 				= os.path.join(DATA_PATH, 'reaction_to_rpair')
	R2K     				= os.path.join(DATA_PATH, 'reaction_to_orthology')
	R2C 					= os.path.join(DATA_PATH, 'reaction_to_compound')
	R2M 					= os.path.join(DATA_PATH, 'reaction_to_module')
	M2R 					= os.path.join(DATA_PATH, 'module_to_reaction')
	R2P 					= os.path.join(DATA_PATH, 'reaction_to_pathway')
	P2R 					= os.path.join(DATA_PATH, 'pathway_to_reaction')
	C2R 					= os.path.join(DATA_PATH, 'compound_to_reaction')
	C   					= os.path.join(DATA_PATH, 'compound_descriptions')    
	R   					= os.path.join(DATA_PATH, 'reaction_descriptions')
	P   					= os.path.join(DATA_PATH, 'pathway_descriptions')

	def __init__(self):
		logging.info("Loading databases")
		logging.debug("Loading module definitions")
		self.m2def = pickle.load(open('.'.join([self.M2DEF,
		                                         self.VERSION, self.PICKLE])))
		logging.debug("Done!")
		logging.debug("Loading module descriptions")
		self.m = pickle.load(open('.'.join([self.M,
		                                    self.VERSION, self.PICKLE])))
		logging.debug("Done!")	
		logging.debug("Loading reaction to pathway information")
		self.r2p = pickle.load(open('.'.join([self.R2P, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading pathway to reaction information")
		self.p2r = pickle.load(open('.'.join([self.P2R, 
		                                         self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading reaction to orthology information")
		self.r2k = pickle.load(open('.'.join([self.R2K, self.VERSION, 
		                                      self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading reaction to module information")
		self.r2m = pickle.load(open('.'.join([self.R2M, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading reaction to module information")
		self.m2r = pickle.load(open('.'.join([self.M2R, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading reaction to compound information")
		self.r2c = pickle.load(open('.'.join([self.R2C, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading reaction to rpair information")
		self.r2rpair = pickle.load(open('.'.join([self.R2RPAIR, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading compound to reaction information")
		self.c2r = pickle.load(open('.'.join([self.C2R,
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading compound descriptions")
		self.c   = pickle.load(open('.'.join([self.C, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading pathway descriptions")
		self.p   = pickle.load(open('.'.join([self.P, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading reaction descriptions")
		self.r   = pickle.load(open('.'.join([self.R, 
		                                      self.VERSION, self.PICKLE])))
		logging.debug("Done")
		logging.debug("Loading compound classifications")
		self.compound_desc_dict \
		         = pickle.load(open('.'.join([self.COMPOUND_DESC_PICKLE, 
		                                      self.VERSION, self.PICKLE])))
		logging.info("Done")
		
		logging.info("Loading reference db paths")		
		self.KO_DB 			= os.path.join(self.DATABASE_DIR, 'uniref100.dmnd')
		self.PFAM_DB 		= os.path.join(self.DATABASE_DIR, 'pfam.hmm')
		self.PFAM_CLAN_DB 	= os.path.join(self.DATABASE_DIR, 'pfam_clans.txt')
		self.TIGRFAM_DB 	= os.path.join(self.DATABASE_DIR, 'tigrfam.hmm')
		logging.info('Done')