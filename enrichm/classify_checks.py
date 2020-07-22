#!/usr/bin/env python3
# pylint: disable=line-too-long
"""
Various functions that apply extra annotation filters to the annotation process
"""
import logging
from itertools import product, chain
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
        logging.debug(f"Running rule checks for {module_name}")
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

    def check_synteny(self, module_name, genome_gff):
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

            genes_on_contigs = dict()

            for gene in parameters['genes']:

                for contig, positions in genome_gff.items():

                    if gene in positions:

                        if contig in genes_on_contigs:
                            genes_on_contigs[contig].append(genome_gff[contig][gene])

                        else:
                            genes_on_contigs[contig] = [genome_gff[contig][gene]]

            result = False
            broke = False

            for contig, gene_positions in genes_on_contigs.items():
                
                if broke:
                    break
                for combination in product(*gene_positions):

                    endings = [gene[1] for gene in combination]
                    
                    clustered_group = list()
                    for subblock in cluster(endings, synteny_range):
                        if len(subblock) >= parameters['minsubblocksize']:
                            clustered_group.append(subblock)

                    if len(clustered_group) <= (parameters['breaks']+1):
                        if len(list(chain(*clustered_group))) >= parameters['min']:
                            result = True
                            broke = True
                            logging.debug(f"Valid syntenous block found for {module_name} ({', '.join(parameters['genes'])}) on contig {contig}")
                            break

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
        genes_in_genome = set(list(chain(*[list(sub.keys()) for sub in genome_gff.values()])))
        required_rules = self.classify_checks.get_required_rules(module_name)
        rule_satisfied = False

        for required_data in required_rules:

            for required_gene_set in required_data:
                number_found = 0
                for required_gene in required_gene_set:
                    if required_gene in genes_in_genome:
                        number_found += 1

                if number_found == len(required_gene_set):
                    logging.debug(f"Required gene for {module_name} ({', '.join(required_gene_set)}) found in genome")
                    rule_satisfied = True
                    break

        return rule_satisfied
