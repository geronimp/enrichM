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
import math
import numpy as np
# Local
from enrichm.kegg_module_grabber import KeggModuleGrabber
from enrichm.module_description_parser import ModuleDescription
###############################################################################

class MetagenomeAnalyzer:
    
    MODULE_AB = 'module_ab'
    KO_PREFIX = 'K'
    
    def __init__(self):
        self.module_ab_output_suffix = '_module_abundances.tsv'
    
    def _extract_kos(self, string):
        output_list = []
        for idx, char in enumerate(string):
            if char == self.KO_PREFIX:
                output_list.append(string[idx:idx+6])
                
        return output_list
    
    def calculate_module_abundance(self, matrix_path, module_to_def, 
                                   needed_kos=None):
        logging.info("Loading KO matrix: %s" % (matrix_path))
        ko_matrix_object = Matrix(matrix_path)
        output_lines = ['\t'.join(["Module", "Sample", 
                                   "Mean", "SD", "SEM"])]
        needed_kos = (needed_kos if needed_kos else ko_matrix_object.rownames)

        logging.info("Generating list of possible KOs for each sample")
        sample_to_possible_kos = {sample:set() for sample in 
                                  ko_matrix_object.colnames}

        for ko in needed_kos:
            for sample in ko_matrix_object.colnames:                
                if float(ko_matrix_object.get_entry(sample, ko))>=0.0:       
                    sample_to_possible_kos[sample].add(ko)
        
        for module, definition in module_to_def.items():
            module_is_possible = False
            output_line_base = [module]
            sample_batches = []
            definition_object = ModuleDescription(definition)
            
            for sample in ko_matrix_object.colnames:
                ko_values = []
                possible_kos_in_sample = set()
                for ko in needed_kos:
                    if float(ko_matrix_object.get_entry(sample, ko))>=0.0:       
                        possible_kos_in_sample.add(ko)
                if(definition_object.num_steps()==definition_object.num_covered_steps(list(sample_to_possible_kos[sample]))):
                    module_is_possible=True
                    for ko in definition_object.kos():  
                        if ko in needed_kos:
                            ko_values.append(float(ko_matrix_object.get_entry(sample, ko)))
                    ko_values_array = np.array(ko_values)
                    mu=ko_values_array.mean()
                    sd=ko_values_array.std()
                    se=sd / math.sqrt(len(ko_values))
                    sample_batches.append([sample, str(mu), 
                                           str(sd), str(se)])
                    module_is_possible = True
                else:
                    sample_batches.append([sample, '0.0', '0.0', '0.0'])
            
            if module_is_possible:
                for batch in sample_batches:
                    line = output_line_base + batch
                    output_lines.append('\t'.join(line))
        return output_lines

    def main(self, args):
        if args.subparser_name == self.MODULE_AB:
            output_path = args.output_prefix + self.module_ab_output_suffix
            kmg = KeggModuleGrabber()
            if(args.modules or args.module_list_file):
                module_to_def = {}
                needed_kos=set()
                if args.modules:
                    for module in args.modules:
                        module_to_def[module] = kmg.m2def[module]
                        for ko in self._extract_kos(kmg.m2def[module]):
                            needed_kos.add(ko)
                elif args.module_list_file:
                    for module in [x.strip() 
                                   for x in open(args.module_list_file)]:
                        module_to_def[module] = kmg.m2def[module]
                        for ko in self._extract_kos(kmg.m2def[module]):
                            needed_kos.add(ko)
            else:
                module_to_def = kmg.m2def
                needed_kos=False
            
            output_lines = self.calculate_module_abundance(args.matrix,
                                                           module_to_def,
                                                           needed_kos)
            
            logging.info('Writing output to file: %s' % (output_path))
            with open(output_path, 'w') as out_io:
                out_io.write('\n'.join(output_lines) + '\n')
            
            
        