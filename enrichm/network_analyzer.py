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
# Local
from enrichm.kegg_matrix import KeggMatrix
from enrichm.network_builder import NetworkBuilder
###############################################################################

class NetworkAnalyser:
    
    MATRIX          = 'matrix'
    NETWORK         = 'network'
    EXPLORE         = 'explore'
    DEGRADE         = 'degrade'
    PATHWAY         = 'pathway'
    ANNOTATE        = 'annotate'
    ENRICHMENT      = 'enrichment'
    MODULE_AB       = 'module_ab'
    TRAVERSE        = 'traverse'


    NETWORK_OUTPUT_FILE  = 'network.tsv'
    METADATA_OUTPUT_FILE = 'metadata.tsv'    
    TRAVERSE_OUTPUT_FILE = 'traverse.tsv'    

    def __init__(self, metadata):
        self.metadata = {}
        for line in open(metadata):
            if line.startswith('#'):continue
            
            split_line = line.strip().split('\t')
            
            if len(split_line) == 1:
                raise Exception("Only one column detected in metadata file, please check that your file is tab separated")
            else:
                sample_id, group = split_line 

            if group in self.metadata:
                self.metadata[group].append(sample_id)
            else:
                self.metadata[group] = [sample_id]

    def _write_results(self, output_path, output_lines):
        '''
        Parameters
        ----------
        output_path: string
            Path to non-existent file to write output lines to
        output_lines: list
            list containing lines to write to output path
        '''
        logging.info('Writing results to file: %s' % output_path)
        with open(output_path, 'w') as output_path_io: 
            output_path_io.write('\n'.join(output_lines))
            output_path_io.flush()
        
    def do(self, matrix, transcriptome, metabolome, depth, filter, limit, queries, 
           subparser_name, starting_compounds, steps, number_of_queries, output_directory):
        '''
        Parameters
        ----------
        depth
        filter
        limit
        metabolome
        queries

        subparser_name
        transcriptome
        output_directory

        '''

        nb = NetworkBuilder(self.metadata.keys())
        km = KeggMatrix(matrix, transcriptome)

        abundances_metagenome = \
                {key:km.group_abundances(self.metadata[key],
                                         km.reaction_matrix) 
                 for key in self.metadata.keys()}

        if transcriptome:
            abundances_transcriptome = \
                    {key:km.group_abundances(self.metadata[key],
                                             km.reaction_matrix_transcriptome) 
                     for key in self.metadata.keys()}            
            abundances_expression = \
                    {key:km.group_abundances(self.metadata[key],
                                             km.reaction_matrix_expression) 
                     for key in self.metadata.keys()}
        else:
            abundances_transcriptome = None
            abundances_expression    = None

        if metabolome:

            abundances_metabolome = km._parse_matrix(metabolome)
            ### ~ TODO: This is a TEMPORARY WORKAROUND
            ### ~ TODO: I've added a note in the help for network analyzer 
            ### ~ TODO: that warns the user about this.
        else:
            abundances_metabolome = None

        if subparser_name==self.TRAVERSE:
            logging.info('Traversing network')
            output_lines = \
                            nb.traverse(abundances_metagenome,
                                        abundances_transcriptome,
                                        limit,
                                        filter,
                                        starting_compounds,
                                        steps,
                                        number_of_queries)
            self._write_results(os.path.join(output_directory, self.TRAVERSE_OUTPUT_FILE), output_lines)

        elif subparser_name==self.EXPLORE:
            logging.info("Using supplied queries (%s) to explore network" \
                                                        % queries)
            network_lines, node_metadata = \
                            nb.query_matrix(abundances_metagenome, 
                                            abundances_transcriptome,
                                            abundances_expression,
                                            queries,
                                            depth)

            self._write_results(os.path.join(output_directory, self.NETWORK_OUTPUT_FILE), network_lines)
            self._write_results(os.path.join(output_directory, self.METADATA_OUTPUT_FILE), node_metadata)

        elif subparser_name==self.PATHWAY:
            logging.info('Generating pathway network')

            network_lines, node_metadata = \
                            nb.pathway_matrix(abundances_metagenome, 
                                              abundances_transcriptome,
                                              abundances_expression,
                                              abundances_metabolome,
                                              limit,
                                              filter)

            self._write_results(os.path.join(output_directory, self.NETWORK_OUTPUT_FILE), network_lines)
            self._write_results(os.path.join(output_directory, self.METADATA_OUTPUT_FILE), node_metadata)