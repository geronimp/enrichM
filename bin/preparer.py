#!/usr/bin/env python2
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
from network_analyzer import NetworkAnalyser

###############################################################################

class Preparer:
    
    REFERENCE_PATH = '/srv/db/uniprot/uniref100/19-09-2016/idmapping.KO.tsv'
    MATRIX_SUFFIX  = '_matrix.tsv'
    
    def __init__(self):        
        logging.info("Loading uniref to orthology information")
        self.reference_dictionary = {}
        for line in open(self.REFERENCE_PATH):
            protein_id, ko_id  = line.strip().split()
            self.reference_dictionary[protein_id] = ko_id
        logging.info('Done!')

    def blast_output(self, 
                     blast_output_paths, 
                     evalue_cutoff, 
                     bitscore_cutoff, 
                     percent_id_cutoff, 
                     output_matrix_path):
        logging.info("Parsing %i blast output files" % \
                                    (len(blast_output_paths)))
        possible_kos = list(set(self.reference_dictionary.values()))
        output_dictionary = {}
        for file in blast_output_paths:
            # Setting up sample column
            output_dictionary[file] = {}
            for ko in possible_kos:
                output_dictionary[file][ko]=0
            # Filling in column
            for line in open(file):
                sline = line.strip().split()
                
                if(float(sline[10]) <= evalue_cutoff and
                   float(sline[11]) >= bitscore_cutoff and
                   float(sline[2])  >= percent_id_cutoff and                   
                   float(sline[2])  >= percent_id_cutoff):
                        hit_ko = self.reference_dictionary[sline[1]]
                        output_dictionary[file][hit_ko] += 1

        logging.info("Writing results to file: %s" % output_matrix_path)
        with open(output_matrix_path, 'w') as output_matrix_path_io:
            header = 'KO\t%s\n' % '\t'.join([os.path.basename(path) for 
                                             path in blast_output_paths])
            output_matrix_path_io.write(header)
            for ko in possible_kos:
                output_line = [ko]
                
                for file in blast_output_paths:
                    output_line.append( str(output_dictionary[file][ko]) )
                
                output_matrix_path_io.write("%s\n" % '\t'.join(output_line))
                    
    def main(self, args):
        output_path = args.output_prefix + self.MATRIX_SUFFIX 
        if args.subparser_name == NetworkAnalyser.MATRIX:
            if args.blast_outputs:
                if(args.percent_aln_query_cutoff or 
                   args.percent_aln_reference_cutoff):
                    logging.warning("Cannot calculate the percent of the query \
aligned to the reference and vice versa with a blast output.")

                self.blast_output(args.blast_outputs,
                                  args.evalue_cutoff,
                                  args.bitscore_cutoff,
                                  args.percent_id_cutoff,
                                  output_path)
