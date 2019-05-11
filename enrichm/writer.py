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

import logging

###############################################################################

class Writer:
    '''
	A collection of functions to write files in various formats.
	'''

    @staticmethod
    def write(output_lines_list, output_path):
        '''
        
        Parameters
        ----------
        output_lines_list   - List. A list of lists. Each sublist is a line, 
                              where each entry is a column entry
        output_path         - String. Path to write output lines to.
        Output
        ------
        '''
        logging.info("Writing results to file: %s" % output_path)

        with open(output_path, 'w') as out_io:

            for output_line_list in output_lines_list:
                output_line_string = '\t'.join([str(column_entry) for column_entry in output_line_list]) + '\n'
                out_io.write(output_line_string)
    
    @staticmethod
	def write_gff(genome, output_file):
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

				features = ['seq_id=%s' % sequence.seqname,
                                    'prodigal_id=%s' % sequence.prod_id,
                                    'partial=%s' % sequence.partial,
                                    'start_type=%s' % sequence.starttype,
                                    'rbs_motif=%s' % sequence.rbs_motif,
                                    'rbs_spacer=%s' % sequence.rbs_spacer,
                                    'gc=%s' % sequence.gc]

				if len(sequence.all_annotations()) > 0:
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
