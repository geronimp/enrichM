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
    
    def get_entry(self, colname, rowname):
        import IPython ; IPython.embed()
        return self.matrix[colname][rowname]


class BuildEncrichmentMatrix:
    
    MATRIX_SUFFIX   = '_enrichment_matrix.tsv'    
    COMPARE_SUFFIX  = '_compare_matrix.tsv'    
    
    def _parse_annotations(self, annotations_path):
        pass
        
    def main(self, annotations, abundances, metadata, output_path):
        annotations_dict = self._parse_annotations(annotations)

        if metadata:
            logging.info("Comparing sets of genomes")
            output_enrichment_matrix = output_path + self.COMPARE_SUFFIX
            metadata = Matrix(metadata)
            logging.info('Done!')

        if abundances:
            logging.info("Generating enrichment matrix")
            output_enrichment_matrix = output_path + self.MATRIX_SUFFIX
            abundances_dict = Matrix(abundances)        
            logging.info('Done!')
        
