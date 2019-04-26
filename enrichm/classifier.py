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

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2017"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near uq.edu.au"
__status__ = "Development"

###############################################################################
# Imports
import os
import re
import logging
from itertools import chain
# Local
from enrichm.databases import Databases
from enrichm.module_description_parser import ModuleDescription
from enrichm.parser import Parser
###############################################################################

class Classify:
    
    KO_OUTPUT       = "module_completeness.tsv"
    MODULE_PATHS    = "module_paths.tsv"
    AGGREGATE_OUTPUT = "aggregate_output.tsv"

    def __init__(self):

        d=Databases()
        self.ko_re              = re.compile('^K\d+$')
        self.signature_modules  = d.signature_modules
        self.m2def              = d.m2def
        self.m                  = d.m

    def _update_with_custom_modules(self, custom_modules):
        custom_modules_dict = dict()
        for line in open(custom_modules):
                custom_modules_dict[line.split('\t')[0]] = line.strip().split('\t')[1]
        self.m2def.update(custom_modules_dict)
        
        for key in custom_modules_dict.keys():
            self.m[key] = 'Custom'
    
    def write(self, lines, output_path):
        
        logging.info('Writing results to file: %s' % output_path)

        with open(output_path, 'w') as output_path_io:
            for line in lines:
                output_path_io.write(line)

    def do(self, custom_modules, cutoff, aggregate, genome_and_annotation_file, 
           genome_and_annotation_matrix, output_directory):
        '''
        
        Parameters
        ----------
        custom_modules                  - string. Path to file containing custom module definitions, 
                                          consistent with KEGG module nomenclature 
                                          (http://www.genome.jp/kegg/module.html)
        cutoff                          - float. Fraction of a module needed in order to be included
                                          in the output.
        genome_and_annotation_file      - string. Path to file containing genome - annotation file. This file 
                                          contains two columns, the first with the genome name, the 
                                          second with a annotation annotation within that genome 
        genome_and_annotation_matrix    - string. Path to file containing genome - annotation matrix
        output                          - string. Path to file to output results to.

        '''
        
        pathway             = dict()
        genome_output_lines = list()

        if custom_modules:
            logging.info('Reading in custom modules: %s' % custom_modules)
            self._update_with_custom_modules(custom_modules)
        
        if genome_and_annotation_file:
            genome_to_annotation_sets = Parser.parse_genome_and_annotation_file_lf(genome_and_annotation_file)
        
        elif genome_and_annotation_matrix:
            genome_to_annotation_sets = Parser.parse_genome_and_annotation_file_matrix(genome_and_annotation_matrix)

        if aggregate:
            logging.info('Reading in abundances: %s' %
                         (genome_and_annotation_matrix))
            abundances = Parser.parse_simple_matrix(genome_and_annotation_matrix, True)
            abundance_result = dict()

            
        logging.info("Read in annotations for %i genomes" % len(genome_to_annotation_sets))
        
        output_lines = ['\t'.join(["Genome_name", "Module_id", "Module_name", "Steps_found", 
                             "Steps_needed", "Percent_steps_found"]) + '\n'] # "KO_found", "KO_needed", "Percent_KO_found"
        
        genome_output_lines = ['\t'.join(["Genome_name", "Module_id", "Module_name"]) + '\n']

        for name, pathway_string in self.m2def.items():
        
            if name not in self.signature_modules:   
                path = ModuleDescription(pathway_string)
                pathway[name] = path

                for genome, annotations in genome_to_annotation_sets.items():    

                    num_covered, _, _, ko_path = path.num_covered_steps(annotations)
                    num_all         = path.num_steps()
                    perc_covered    = num_covered / float(num_all)
                    ko_path_list    = list(chain(*ko_path.values()))

                    if perc_covered >= cutoff:

                        if path.is_single_step:
                        
                            if perc_covered!=1:

                                if cutoff < 1:
                                    num_all = 1
                                    num_covered = 0
                                    perc_covered = 0.0

                                else:
                                    continue
                            
                            else:
                                num_all = 1
                                num_covered = 1

                        if aggregate:
                            
                            if genome not in abundance_result:
                                abundance_result[genome] = dict()
                            pathway_abundance = [abundances[genome][ko] for ko in ko_path_list]
                            pathway_average_abundance = sum(pathway_abundance) / len(pathway_abundance)
                            abundance_result[genome][name] = pathway_average_abundance

                        genome_output_lines.append('\t'.join([genome, name, self.m[name], 
                                                              ','.join(ko_path_list) + '\n']))
                        output_line = "\t".join([genome, name, self.m[name],
                                                  str(num_covered), 
                                                  str(num_all),
                                                  str(round(perc_covered * 100, 2))
                                                  ]) 
                        output_lines.append(output_line + '\n') 
                        
        self.write(output_lines, os.path.join(output_directory, self.KO_OUTPUT))
        self.write(genome_output_lines, os.path.join(output_directory, self.MODULE_PATHS)) 
        
        if aggregate:
            samples = list(abundance_result.keys() )
            output_lines = ['\t'.join(["ID"] + samples) + '\n']
            for module in self.m2def.keys():
            
                if module not in self.signature_modules:
                    ol = [module]

                    for sample in samples:
            
                        if module in abundance_result[sample]:
                            ol.append(str(abundance_result[sample][module]))
            
                        else:
                            ol.append('0.0')
                    output_lines.append('\t'.join(ol) + '\n')
            self.write(output_lines, os.path.join(output_directory, self.AGGREGATE_OUTPUT))



