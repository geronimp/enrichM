#!/usr/bin/env python3
# pylint: disable=line-too-long
"""
Various functions that apply extra annotation filters to the annotation process
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
import logging
from itertools import product
from enrichm.toolbox import window, cluster

class ClassifyChecks:
    '''
    Applies extra, manually defined filters to the 'Classify' annotation process.
    '''
    def __init__(self, classify_checks):
        self.classify_checks = classify_checks

    def check(self, module_name, genome_gff):
        """ Description
        :type self:
        :param self:

        :type module_name:
        :param module_name:

        :type genome_gff:
        :param genome_gff:

        :raises:

        :rtype:
        """
        if module_name in self.classify_checks.members:
            if not self.check_synteny(module_name, genome_gff):
                logging.debug(f"Module {module_name} in genome failed synteny check")
                return False
            elif not self.check_required(module_name, genome_gff):
                logging.debug(f"Module {module_name} in genome failed required check")
                return False
            else:
                logging.debug(f"Module {module_name} passed all checks in genome")
                return True # in the clear

        else:
            return True

    def check_synteny(self, module_name, genome_gff): # TODO: How to account for spread across multiple contigs?
    
        """ Description
        :type self:
        :param self:
    
        :type module_name:
        :param module_name:
    
        :type genome_gff:
        :param genome_gff:
    
        :raises:
    
        :rtype:
        """
        synteny_rules = self.classify_checks.get_synteny_rules(module_name)
        result_list = list()
        
        for parameters in synteny_rules:
            
            synteny_range = parameters['range']
            if parameters['strict']: # TODO: this and others like it shouldn't be a string
                candidates = list()
                for idx, (gene_1, gene_2) in enumerate(window(parameters['genes'], 2)):
                    
                    gene_1_positions = genome_gff[gene_1]
                    gene_2_positions = genome_gff[gene_2]

                    if idx!=0:

                        if len(candidates)==0:
                            break

                    for left_gene_position in gene_1_positions:

                        if idx==0:
                            ends = [left_gene_position[1]] # Extract end position of first gene
                            logging.debug(f"Number of candidates for synteny: {len(ends)} in step {idx+1}")

                        else:
                            logging.debug(f"Number of candidates for synteny: {len(candidates)} in step {idx+1}")
                            ends = candidates
                            candidates = list()

                        for right_gene_position in gene_2_positions:
                            start = right_gene_position[0] # Extract start position of second gene

                            for end in ends:

                                if ((end-start) < synteny_range and (end-start) > -synteny_range):

                                    if right_gene_position[1] not in candidates:
                                        candidates.append(right_gene_position[1]) # remember the possible candidates for linkages

                    if len(candidates)==0:
                        result = False
                    else:
                        result = True

            else:
                gene_positions = [genome_gff[gene] for gene in parameters['genes']]
                for combination in product(*gene_positions):
                    endings = [gene[1] for gene in combination]
                    clustered_group = cluster(endings, synteny_range)
                    if len(clustered_group) == 1:
                        result = True
                        break
                    else:
                        result = False

            result_list.append(result)

        if all(result_list):
            return True
        else:
            return False

    def check_required(self, module_name, genome_gff):
    
        """ Description
        :type self:
        :param self:
    
        :type module_name:
        :param module_name:
    
        :type genome_gff:
        :param genome_gff:
    
        :raises:
    
        :rtype:
        """
        required_rules = self.classify_checks.get_required_rules(module_name)
        rule_satisfied = False
        for required_data in required_rules:

            for required_gene_set in required_data:

                number_found = len([required_gene for required_gene in required_gene_set if required_gene in genome_gff])

                if number_found == len(required_gene_set):
                    logging.debug(f"Required gene for {module_name} ({', '.join(required_gene_set)}) found in genome")
                    rule_satisfied = True
                    break

        return rule_satisfied
