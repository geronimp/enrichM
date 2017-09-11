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

import logging
import os
import pickle

from genome import Genome, AnnotationParser
from annotate import Annotate

###############################################################################

class Compare:


	def _parse_pickles(self, enrichm_annotate_output):
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
		output_genome_list = []

		genome_pickle_file_path \
			= os.path.join(enrichm_annotate_output, Annotate.GENOME_OBJ)

		for pickled_genome in os.listdir(genome_pickle_file_path):
			pickled_genome_path = os.path.join(genome_pickle_file_path, pickled_genome)
			logging.debug('Parsing genome: %s' % (pickled_genome_path))
			output_genome_list.append(pickle.load(open(pickled_genome_path)))

		return output_genome_list

	def do(self, enrichm_annotate_output):
		'''
		### ~ TODO: Not sure what this does yet.		

		Parameters
		----------
		enrichm_annotate_output - String. Output directory of a previous run 
								  of enrichm annotate (At lease version 0.0.7)
								  or above
		'''	
		logging.info('Parsing pickled genomes from previous enrichm run: %s' \
						% (enrichm_annotate_output))
		genome_list = self._parse_pickles(enrichm_annotate_output)
		
