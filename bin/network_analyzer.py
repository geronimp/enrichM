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

from kegg_matrix import KeggMatrix
from network_builder import NetworkBuilder
import argparse
import logging
import os

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################
   
    
class NetworkAnalyser:
    
    NETWORK_SUFFIX  = '_network.tsv'
    METADATA_SUFFIX = '_metadata.tsv'
    
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
        
    def main(self, matrix_path, queries, depth, transcriptome, output):
        '''
        Parameters
        ----------
        matrix_path: string
            Path to file containing a KO matrix build from either metagenomic
            or metatranscriptomic data.
        queries: string
            Path to file containing query compounds
        depth: integer
            Number of steps into metabolism to take if 'queries' is provided
        transcriptome: string
            Path to file containing a KO matrix of transcriptomic data
        '''
        
        nb = NetworkBuilder()
        
        if transcriptome:
            km = KeggMatrix(matrix_path, transcriptome)
        else:
            km = KeggMatrix(matrix_path)

        if len(self.metadata.keys())==2:
            group1_abundances = \
                    km.group_abundances(self.metadata[self.metadata.keys()[0]],
                                        km.reaction_matrix)
            group2_abundances = \
                    km.group_abundances(self.metadata[self.metadata.keys()[1]],
                                        km.reaction_matrix)
            if transcriptome:
                
                group1_transcriptome_abundances = \
                    km.group_abundances(self.metadata[self.metadata.keys()[0]],
                                        km.reaction_matrix_transcriptome)
                group2_transcriptome_abundances = \
                    km.group_abundances(self.metadata[self.metadata.keys()[1]],
                                        km.reaction_matrix_transcriptome)
                
                group1_expression_abundances = \
                    km.group_abundances(self.metadata[self.metadata.keys()[0]],
                                        km.reaction_matrix_expression)
                group2_expression_abundances = \
                    km.group_abundances(self.metadata[self.metadata.keys()[1]],
                                        km.reaction_matrix_expression)
            else:
                group1_transcriptome_abundances=None
                group2_transcriptome_abundances=None
                group1_expression_abundances=None
                group2_expression_abundances=None
        else:
            raise Exception("network_analyzer does not yet handle comparisons \
between >2 groups")
        logging.info("Constructing network for input matrix")
        if queries:
            logging.info("Using supplied queries (%s) to explore network" \
                                                    % queries)
            output_lines, node_metadata = nb.query_matrix(group1_abundances, 
                                           group2_abundances,
                                           queries,
                                           depth,
                                           group1_expression_abundances,
                                           group2_expression_abundances,
                                           group1_transcriptome_abundances,
                                           group2_transcriptome_abundances,
                                           self.metadata.keys()[0],
                                           self.metadata.keys()[1])
        else:
            output_lines, node_metadata = nb.all_matrix(group1_abundances,
                                         group2_abundances)
        self._write_results(output + self.NETWORK_SUFFIX, output_lines)
        self._write_results(output + self.METADATA_SUFFIX, node_metadata)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Build a metabolic matrix''')
    parser.add_argument('--matrix', required=True,
                        help='KO matrix')
    parser.add_argument('--transcriptome', 
                        help='metagenome to ko matrix')
    parser.add_argument('--queries', 
                        help='query compounds')
    parser.add_argument('--depth', type=int, default=2,
                        help='depth')
    parser.add_argument('--metadata', required=True,
                        help='description of samples')
    parser.add_argument('--log',
                        help='output logging information to this file.')
    parser.add_argument('--verbosity', type = int, default = 4,
                        help='Level of verbosity (1 - 5 - default = 4) \
5 = Very verbose, 1 = Silent')
    parser.add_argument('--output', default = 'network_analyser_output',
                        help='Output directory or file')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite previous run')
    args = parser.parse_args()
        
    if args.log:
        if os.path.isfile(args.log): 
            raise Exception("File %s exists" % args.log)
        logging.basicConfig(filename=args.log, level=debug[args.verbosity], 
                            format='%(asctime)s %(levelname)s: %(message)s', 
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=debug[args.verbosity], 
                            format='%(asctime)s %(levelname)s: %(message)s', 
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    
    if os.path.isfile(args.output):
        if args.force:
            logging.warning("Removing existing file with name: %s" \
                                                                % args.output)
            os.remove(args.output)
        else:
            raise Exception("File %s exists" % args.output)
    
    na=NetworkAnalyser(args.metadata)
    na.main(args.matrix, args.queries, args.depth, args.transcriptome, 
            args.output)