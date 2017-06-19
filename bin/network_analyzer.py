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
from kegg_matrix import KeggMatrix
from network_builder import NetworkBuilder
from build_enrichment_matrix import Matrix

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


    NETWORK_SUFFIX  = '_network.tsv'
    METADATA_SUFFIX = '_metadata.tsv'    
    TRAVERSE_SUFFIX = '_traverse.tsv'    

    def __init__(self, metadata):
        self.metadata = {}
        for line in open(metadata):
            sample_id, group = line.strip().split('\t')
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
        
    def do(self, depth, filter, limit, metabolome, number_of_queries, queries, 
           starting_compounds, steps, subparser_name, transcriptome, output_prefix):
        '''
        Parameters
        ----------
        depth
        filter
        limit
        metabolome
        number_of_queries
        queries
        starting_compounds
        steps
        subparser_name
        transcriptome
        output_prefix

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
            abundances_metabolome = Matrix(metabolome)
        else:
            abundances_compounds     = None

        if subparser_name==self.TRAVERSE:
            output_lines = \
                            nb.traverse(abundances_metagenome,
                                        abundances_transcriptome,
                                        limit,
                                        filter,
                                        starting_compounds,
                                        steps,
                                        number_of_queries)
            self._write_results(output_prefix + self.TRAVERSE_SUFFIX, output_lines)

        elif subparser_name==self.EXPLORE:
            logging.info("Using supplied queries (%s) to explore network" \
                                                        % queries)
            network_lines, node_metadata = \
                            nb.query_matrix(abundances_metagenome, 
                                            abundances_transcriptome,
                                            abundances_expression,
                                            queries,
                                            depth)

            self._write_results(output_prefix + self.NETWORK_SUFFIX, network_lines)
            self._write_results(output_prefix + self.METADATA_SUFFIX, node_metadata)

