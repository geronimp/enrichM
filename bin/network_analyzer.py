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

def check_args(args):
    if not(args.queries):
        if args.depth:
            logging.warning("--depth argument ignored without --queries flag")
    if os.path.isfile(args.output):
        if args.force:
            logging.warning("Removing existing file with name: %s" \
                                                                % args.output)
            os.remove(args.output)
        else:
            raise Exception("File %s exists" % args.output)

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
        
        nb = NetworkBuilder(self.metadata.keys())
        km = KeggMatrix(matrix_path, transcriptome)

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

        
        logging.info("Constructing network for input matrix")
        if queries:
            logging.info("Using supplied queries (%s) to explore network" \
                                                    % queries)
            network_lines, node_metadata = \
                            nb.query_matrix(abundances_metagenome, 
                                            abundances_transcriptome,
                                            abundances_expression,

                                            queries,
                                            depth)
        else:
            logging.info("Constructing entire metabolic network. Be patient.")
            network_lines, node_metadata = \
                            nb.all_matrix(abundances_metagenome,
                                          abundances_transcriptome,
                                          abundances_expression)
                    
        self._write_results(output + self.NETWORK_SUFFIX, network_lines)
        self._write_results(output + self.METADATA_SUFFIX, node_metadata)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Build a metabolic matrix''')
    parser.add_argument('--matrix', required=True,
                        help='KO matrix')
    parser.add_argument('--metadata', required=True,
                        help='description of samples')
    parser.add_argument('--transcriptome', 
                        help='metagenome to ko matrix')
    parser.add_argument('--queries', 
                        help='query compounds')
    parser.add_argument('--depth', type=int, default=2,
                        help='depth')
    parser.add_argument('--log',
                        help='output logging information to this file.')
    parser.add_argument('--verbosity', type = int, default = 4,
                        help='Level of verbosity (1 - 5 - default = 4) \
5 = Very verbose, 1 = Silent')
    parser.add_argument('--output', default = 'network_analyser_output',
                        help='Output file')
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
    
    check_args(args)
    
    na=NetworkAnalyser(args.metadata)
    na.main(args.matrix, args.queries, args.depth, args.transcriptome, 
            args.output)