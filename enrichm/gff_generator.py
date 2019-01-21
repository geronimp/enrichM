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
# System imports

class GffGenerator():

	def write(self, genome, output_file):
		'''
		Writes a gff file for a Genome object
		
		Parameters
		----------
		genome 			- Genome object
		output_file 	- string. file name to output results to.
		'''

		with open(output_file, 'w') as out_io:
			for sequence in genome.ordered_sequences():
				contig_id = '_'.join(sequence.seqname.split('_')[:-1])
				
				features  = ['seq_id=%s'		% sequence.seqname,
						 	 'prodigal_id=%s' 	% sequence.prod_id,
							 'partial=%s' 	 	% sequence.partial,
							 'start_type=%s'  	% sequence.starttype,
							 'rbs_motif=%s'   	% sequence.rbs_motif,
							 'rbs_spacer=%s'  	% sequence.rbs_spacer,
							 'gc=%s'    		% sequence.gc]

				if len(sequence.all_annotations())>0:
					features.append('annotations=%s' % ','.join(sequence.all_annotations()))
				else:
					features.append('annotations=hypothetical_protein')

				line = [contig_id, 
						'prodigal',
						'CDS',
						sequence.startpos,
						sequence.finishpos,
						'.',
						('-' if sequence.direction == '-1'
						 else '+'),
						'0',
						';'.join(features)]
						
				out_io.write('\t'.join(line) + '\n')