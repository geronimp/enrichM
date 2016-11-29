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
import pickle
import shutil

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################

def check_args(args):
    '''
    '''    
    # Check output file doesn't already exist. Deleted if present and --force
    # is specified
    if (os.path.isfile(args.output) or os.path.isdir(args.output)):
        if args.force:
            logging.warning("Overwriting output file that already exists with \
name %s" % args.output)
            shutil.rmtree(args.output)
        else:
            raise Exception("Output file already exists: %s" % args.output)
    
    
class NetworkAnalyser:
    
    def __init__(self, metadata, output):
        
        self.metagenome_matrix    = 'metagenome.tsv'
        self.transcriptome_matrix = 'transcriptome.tsv'
        self.combined_matrix      = 'combined.tsv'
        
        self.output_directory     = output
        os.mkdir(self.output_directory)
        
        self.metadata = {}
        
        for line in open(metadata):
            sample_id, group = line.strip().split('\t')
            if group in self.metadata:
                self.metadata[group].append(sample_id)
            else:
                self.metadata[group] = [sample_id]

    def _write_results(self, output_path, output_lines):
        '''
        '''
        with open(os.path.join(self.output_directory,
                                output_path), 'w') as output_path_io: 
            output_path_io.write('\n'.join(output_lines))
            output_path_io.flush()
        


    def main(self, m_matrix_path, queries, depth):
        '''
        '''
        
        nb = NetworkBuilder()
        
        logging.info("Parsing input matrix: %s" % m_matrix_path)
        m_km = KeggMatrix(m_matrix_path)
        group1_abundances_mg = \
                        m_km.group_abundances(KeggMatrix.REACTION, 
                                              self.metadata['Eriophorum'])
        group2_abundances_mg = \
                        m_km.group_abundances(KeggMatrix.REACTION, 
                                              self.metadata['Sphagnum'])
        logging.info("Constructing matrix for %s metagenome abundances" \
                                                    % KeggMatrix.REACTION)
        if queries:
            logging.info("Using supplied queries: %s" \
                                                    % queries)
            output_lines = nb.query_matrix(group1_abundances_mg, 
                                           group2_abundances_mg,
                                           queries,
                                           depth)
        else:
            output_lines = nb.all_matrix(group1_abundances_mg,
                                         group2_abundances_mg)
        self._write_results('_'.join(['all', self.metagenome_matrix]), 
                            output_lines)

        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''Build a metabolic matrix''')
    parser.add_argument('--matrix', required=True,
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
        raise Exception("File %s exists" % args.output)
    
    check_args(args)

    na=NetworkAnalyser(args.metadata, args.output)
    na.main(args.matrix, args.queries, args.depth)