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

################################################################################

class Parser:
	
	@staticmethod
	def parse_genome_and_annotation_file_lf(self, genome_and_annotation_file):
		genome_to_annotation_sets = dict()
		
		for line in open(genome_and_annotation_file):
			try:
				genome, annotation = line.strip().split("\t")
			except:
				raise Exception("Input genomes/annotation file error on %s" % line)
			
			if genome not in genome_to_annotation_sets:
				genome_to_annotation_sets[genome] = set()

			genome_to_annotation_sets[genome].add(annotation)
		
		return genome_to_annotation_sets
	
	@staticmethod
	def parse_genome_and_annotation_file_matrix(self, genome_and_annotation_matrix):
		genome_and_annotation_matrix_io = open(genome_and_annotation_matrix)
		headers=genome_and_annotation_matrix_io.readline().strip().split('\t')[1:]
		genome_to_annotation_sets = {genome_name:set() for genome_name in headers}

		for line in genome_and_annotation_matrix_io:
			sline = line.strip().split('\t')
			annotation, entries = sline[0], sline[1:]
			for genome_name, entry in zip(headers, entries):
				if float(entry) > 0:
					genome_to_annotation_sets[genome_name].add(annotation)

		return genome_to_annotation_sets
	
	@staticmethod
	def parse_taxonomy(self, taxonomy_path):
		
		output_taxonomy_dictionary = dict()

		for line in open(taxonomy_path):
			genome, taxonomy_string = line.strip().split('\t')
			output_taxonomy_dictionary[genome] = taxonomy_string.split(';')
			
		return output_taxonomy_dictionary

	
	@staticmethod
    def _parse_matrix(self, matrix):
        
        output_dict = {}
        
        for idx, line in enumerate(open(matrix)):
            sline = line.strip().split('\t')
            if idx==0:
                self.sample_names = sline[1:]
                    
                for sample in self.sample_names: 
                    output_dict[sample] = {}
            else:   
                ko_id      = sline[0]
                abundances = sline[1:]
                
                for abundance, sample in zip(abundances, self.sample_names):
                    output_dict[sample][ko_id] = float(abundance)
        return output_dict
    
    @staticmethod
    def _parse_matrix(self, matrix_file_io, colnames):
        for line in matrix_file_io:
            sline = line.strip().split('\t')
            rowname, entries = sline[0], sline[1:]
            for colname, entry in zip(colnames, entries):
                yield colname, entry, rowname