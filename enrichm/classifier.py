#!/usr/bin/env python3
# pylint: disable=line-too-long
"""
Various functions that take genome annotations, and determine the metabolic pathways they encode.
"""
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

# Imports
import os
import logging
from itertools import chain
# Local
from enrichm.databases import Databases
from enrichm.module_description_parser import ModuleDescription
from enrichm.parser import Parser, RulesJson
from enrichm.writer import Writer
from enrichm.toolbox import get_present_annotations
from enrichm.classify_checks import ClassifyChecks
###############################################################################

class Classify:
    '''
    Determines which metabolic pathways are encoded by different MAGs
    '''

    def __init__(self):
        databases = Databases()
        self.signature_modules = databases.signature_modules
        self.m2def = databases.m2def()
        self.modules = databases.m()
        self.ko_output = "module_completeness.tsv"
        self.module_paths = "module_paths.tsv"
        self.aggregate_output = "aggregate_output.tsv"

    def update_with_custom_modules(self, custom_modules):
        '''
        open a file
        :param custom_modules: Path to a file containing custom modules to add to the list of 
                               default modules to use.
        :type custom_modules: str
        '''
        custom_modules_dict = dict()

        with open(custom_modules) as custom_modules_io:
            for line in custom_modules_io:
                custom_modules_dict[line.split('\t')[0]] = line.strip().split('\t')[1]

        for key in custom_modules_dict:
            self.modules[key] = 'Custom'

        return custom_modules_dict

    def classify_pipeline(self, custom_modules, cutoff, aggregate, genome_and_annotation_matrix,
                          module_rules_json, gff_files, output_directory):
        '''

        Parameters
        ----------
        custom_modules                  - string. Path to file containing custom module definitions,
                                          consistent with KEGG module nomenclature
                                          (http://www.genome.jp/kegg/module.html)
        cutoff                          - float. Fraction of a module needed in order to be included
                                          in the output.
        genome_and_annotation_matrix    - string. Path to file containing genome - annotation matrix
        output                          - string. Path to file to output results to.

        '''
        pathway = dict()
        genome_output_lines = list()

        if module_rules_json:
            cc = ClassifyChecks(RulesJson().load(module_rules_json))

        if custom_modules:
            logging.info('Reading in custom modules: %s' % custom_modules)
            modules_to_classify = self.update_with_custom_modules(custom_modules)
        else:
            modules_to_classify = self.m2def

        if gff_files:
            genome_to_annotation_sets = dict()
            annotation_results = dict()
            for gff in gff_files:
                feature_dict, genome_to_annotations_set = Parser.parse_gff(gff)
                genome_to_annotation_sets[gff] = genome_to_annotations_set
                annotation_results[gff] = feature_dict
        else:
            genome_to_annotation_sets, _, _ = Parser.parse_simple_matrix(genome_and_annotation_matrix)

        if aggregate:
            logging.info('Reading in abundances: %s' % (genome_and_annotation_matrix))
            abundances, _, _ = Parser.parse_simple_matrix(genome_and_annotation_matrix)
            abundance_result = dict()

        logging.info("Read in annotations for %i genomes" % len(genome_to_annotation_sets))

        output_lines = [["Genome_name", "Module_id", "Module_name", "Steps_found",
                                   "Steps_needed", "Percent_steps_found"]]

        genome_output_lines = [["Genome_name", "Module_id", "Module_name"]]

        for name, pathway_string in modules_to_classify.items():

            if name not in self.signature_modules:
                path = ModuleDescription(pathway_string)
                pathway[name] = path

                for genome, annotation_frequency in genome_to_annotation_sets.items():
                    annotations = get_present_annotations(annotation_frequency)
                    num_covered, _, _, ko_path = path.num_covered_steps(annotations)
                    num_all = path.num_steps()
                    perc_covered = num_covered / float(num_all)
                    ko_path_list = list(chain(*ko_path.values()))

                    if perc_covered >= cutoff:

                        rule_check_result = cc.check(name, annotation_results[genome])

                        if rule_check_result:

                            if path.is_single_step:

                                if perc_covered != 1:

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

                            if len(pathway_abundance)>0:
                                pathway_average_abundance = sum(pathway_abundance) / len(pathway_abundance)
                            else:
                                pathway_average_abundance = 0
                            abundance_result[genome][name] = pathway_average_abundance

                        genome_output_lines.append([genome, name, self.modules[name], ','.join(ko_path_list)])
                        output_line = [genome, name, self.modules[name], str(num_covered), str(num_all), str(round(perc_covered * 100, 2))]
                        output_lines.append(output_line)

        Writer.write(output_lines, os.path.join(output_directory, self.ko_output))
        Writer.write(genome_output_lines, os.path.join(output_directory, self.module_paths))

        if aggregate:
            samples = list(abundance_result.keys())
            output_lines = ['\t'.join(["ID"] + samples) + '\n']

            for module in modules_to_classify.keys():

                if module not in self.signature_modules:
                    output_line = [module]

                    for sample in samples:

                        if module in abundance_result[sample]:
                            output_line.append(str(abundance_result[sample][module]))

                        else:
                            output_line.append('0.0')
                    output_lines.append(output_line)

            Writer.write(output_lines, os.path.join(output_directory, self.aggregate_output))
