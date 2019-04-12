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
import logging
import os
import pickle
import multiprocessing as mp
# Local
from enrichm.annotate import Annotate
################################################################################

def parse_genomes(path):
	logging.info('Loading: %s' % (os.path.basename(path)))
	genome = pickle.load(open(path, 'rb'))
	return genome

################################################################################

class ParseAnnotate:
	
	"""docstring for ParseAnnotate"""

	def __init__(self, enrichm_annotate_output, processes):
		self.path = enrichm_annotate_output
		# Parse genome objects
		self.genome_pickle_file_path \
			= os.path.join(enrichm_annotate_output, Annotate.GENOME_OBJ)
		self.processes \
			= processes

		ko = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_KO)
		if os.path.isfile(ko):
			self.ko = ko
		else:
			self.ko = None
		ko_hmm = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_KO_HMM)
		if os.path.isfile(ko_hmm):
			self.ko_hmm = ko_hmm
		else:
			self.ko_hmm = None
		pfam = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_PFAM)
		if os.path.isfile(pfam):
			self.pfam = pfam
		else:
			self.pfam = None
		tigrfam = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_TIGRFAM)
		if os.path.isfile(tigrfam):
			self.tigrfam = tigrfam
		else:
			self.tigrfam = None
		cazy = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_CAZY)
		if os.path.isfile(cazy):
			self.cazy = cazy
		else:
			self.cazy = None
		ec = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_EC)
		if os.path.isfile(ec):
			self.ec = ec
		else:
			self.ec = None
		hypothetical_cluster = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_HYPOTHETICAL_CLUSTER)
		if os.path.isfile(hypothetical_cluster):
			self.hypothetical_cluster = hypothetical_cluster
		else:
			self.hypothetical_cluster = None
		hypothetical_ortholog = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_HYPOTHETICAL_ORTHOLOG)
		if os.path.isfile(hypothetical_ortholog):
			self.hypothetical_ortholog = hypothetical_ortholog
		else:
			self.hypothetical_ortholog = None

	def parse_pickles(self, path, genome_list):
		'''
		Opens the pickled genome objects from a previous run of enrichm 
		annotate

		Parameters
		----------
		enrichm_annotate_output - String. Output directory of a previous run 
								  of enrichm annotate (At lease version 0.0.7)
								  or above
		Outputs
		-------
		List of Genome objects
		'''	
		logging.info('Loading genome pickles')
		
		output_genome_list 	= list()
		paths 				= list()

		for pickled_genome in genome_list:
			pickled_genome_path = os.path.join(path, pickled_genome + '.pickle')
			if os.path.isfile(pickled_genome_path):
				paths.append(pickled_genome_path)
		
		self.pool = mp.Pool(processes = self.processes)
		output_genome_list = self.pool.map_async(parse_genomes, paths)
		output_genome_list.wait()
		genome_objects = output_genome_list.get()
		self.pool.close()
		return genome_objects
