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
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

###############################################################################

import logging
import os
from itertools import product

###############################################################################

class Matrix:
    suffix_to_delim = {'.tsv':'\t',
                       '.csv':','}
    
    def __init__(self, matrix):
        self.rownames = []
        self.colnames = []
        self.matrix = {}
        self._parse_matrix(open(matrix))
        
    def _parse_matrix(self, matrix_io):
        delim = self.suffix_to_delim[os.path.splitext(matrix_io.name)[-1]]
        headers = matrix_io.readline().strip().split(delim)[1:]
        for header in headers:
            self.matrix[header]={}
            if header in self.colnames:
                raise Exception("Duplicate column names found in %s" % \
                                                         (matrix_io.name))
            self.colnames.append(header)
        for row in matrix_io:
            split_row         = row.strip().split(delim)
            row_name, entries = split_row[0], split_row[1:]
            for column_entry in split_row:
                for column_name, entry in zip(self.colnames, entries):
                    self.matrix[column_name][row_name]=entry
            if row_name in self.rownames:
                raise Exception("Duplicate row names found in %s" % \
                                                         (matrix_io.name))
            self.rownames.append(row_name)

    def get_row(self, rowname):
        if rowname in self.rownames:
            row = [rowname]
            for colname in self.colnames: 
                row.append(self.matrix[colname][rowname])
        return row
    
    def get_col(self, colname):      
        if colname in self.colnames:
            
            return self.matrix[colname].values()
        else:
            return 0 # raise Exception("Colname %s not in matrix" % colname)

    def get_entry(self, colname, rowname):
        if colname in self.colnames:
            if rowname in self.rownames:
                return self.matrix[colname][rowname]
            else:
                return 0 # raise Exception("Rowname %s not in matrix" % rowname)
        else:
            return 0 # raise Exception("Colname %s not in matrix" % colname)
    
    def filter_by_cols(self, set):
        current_list = []
        for row in self.rownames:
            if tuple(self.get_row(row)[1:]) == set:
                current_list.append(row)
        return current_list
        
class BuildEncrichmentMatrix:
    
    BACKGROUND      = 'background'
    MATRIX_SUFFIX   = '_enrichment_matrix.tsv'
    SAMPLE_MATRIX_SUFFIX   = '_sample_enrichment_matrix.tsv'
    COMPARE_SUFFIX  = '_compare_matrix.tsv'    
    
    def _parse_annotations(self, annotations_path):
        modules              = set()
        genomes              = set()
        genome_to_annotation = {}
        
        for line in open(annotations_path):
            sline = line.strip().split('\t')
            genome_id, module_id = sline[:2]
            modules.add(module_id)
            genomes.add(genome_id)
            if genome_id in genome_to_annotation:
                genome_to_annotation[genome_id].add(module_id)
            else:
                genome_to_annotation[genome_id] = set([module_id])
        return genome_to_annotation, modules, genomes
        
    def main(self, annotations, abundances, metadata, subset_modules, 
             output_path):
        annotations_dict, modules, genomes = self._parse_annotations(annotations)
        
        if subset_modules:
            modules = subset_modules
        
        logging.info("Comparing sets of genomes")
        output_enrichment_matrix = output_path + self.COMPARE_SUFFIX
        metadata = Matrix(metadata)
        metadata_value_lists = \
            [set(metadata.get_col(col)) for col in metadata.colnames]
        
        combination_dict = {self.BACKGROUND: genomes}
        for combination in product(*metadata_value_lists):
            genome_list = metadata.filter_by_cols(combination)
            combination_dict['_'.join(combination)]=genome_list

        output_lines = ['\t'.join(['Module'] + combination_dict.keys())]
        for module in modules:
            output_line = [module]
            for group_name, genome_list in combination_dict.items():
                if len(genome_list)>0:
                    coverage = len([genome for genome in genome_list 
                                    if module in annotations_dict[genome]])
                    all      = float(len(genome_list))
                    output_line.append(str(coverage/all))
                else:
                    output_line.append('0.0')
            output_lines.append('\t'.join(output_line))
            
        logging.info("Writing results to file: %s" % output_enrichment_matrix)
        with open(output_enrichment_matrix, 'w') as output_enrichment_matrix_io:
            output_enrichment_matrix_io.write('\n'.join(output_lines))
        logging.info('Done!')

        if abundances:
            logging.info("Generating enrichment matrix")
            abundances = Matrix(abundances)
            
            output_sample_enrichment_matrix = \
                                        output_path + self.SAMPLE_MATRIX_SUFFIX
            output_lines = ['\t'.join(['Module'] + abundances.colnames)]
            for module in modules:
                output_line = [module]
                for sample in abundances.colnames:
                    
                    module_prevalence = 0.0
                    for genome in genomes:
                        if module in annotations_dict[genome]:
                            genome_abundance = abundances.get_entry(sample, genome)
                            module_prevalence+=float(genome_abundance)
                    output_line.append(str(module_prevalence))
                output_lines.append('\t'.join(output_line))
            
            logging.info("Writing results to file: %s" \
                            % output_sample_enrichment_matrix)
            with open(output_sample_enrichment_matrix, 'w') as output_enrichment_matrix_io:
                output_enrichment_matrix_io.write('\n'.join(output_lines))
            logging.info('Done!')
            
            output_enrichment_matrix = output_path + self.MATRIX_SUFFIX
            output_lines = ['\t'.join(['Module', 'Sample', 
                                       'Group', 'Abundance'])]
            for module in modules:
                for sample in abundances.colnames:
                    for group in product(*metadata_value_lists):
                        group_name='_'.join(group)
                        genomes_in_group_list = metadata.filter_by_cols(group)                        
                        module_prevalence = 0.0
                        for genome in genomes_in_group_list:
                            if module in annotations_dict[genome]:
                                genome_abundance = abundances.get_entry(sample, genome)
                                module_prevalence+=float(genome_abundance)
                        output_line = [module, sample, 
                                       group_name, str(module_prevalence)]
                        
                        output_lines.append('\t'.join(output_line))
                
            logging.info("Writing results to file: %s" % output_enrichment_matrix)
            with open(output_enrichment_matrix, 'w') as output_enrichment_matrix_io:
                output_enrichment_matrix_io.write('\n'.join(output_lines))
            logging.info('Done!')
            

            
        
