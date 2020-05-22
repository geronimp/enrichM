#!/usr/bin/env python3
import logging
import os
import statistics
import itertools
from enrichm.network_builder import NetworkBuilder
from enrichm.databases import Databases
from enrichm.parser import Parser
from enrichm.writer import Writer
###############################################################################

class NetworkAnalyser:
    """
    Prepare metagenome, metatranscriptome, metabolomic data for constructing 
    SIF network files.
    """
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

    def __init__(self):
        self.databases = Databases()
        self.reactions = self.databases.r()
        self.reaction_to_ko = self.databases.r2k()

    def average(self, input_dictionary):
        '''
        Take the average of the values of a dictionary of dictionaries
        '''
        
        for sample_group, group_dict in input_dictionary.items():

            for group, reaction_dict in group_dict.items():

                for reaction, value in reaction_dict.items():
                    input_dictionary[sample_group][group][reaction] = sum(value) / len(value)

        return input_dictionary

    def median_genome_abundance(self, sample_abundance_dict, sample_metadata):
        """
        Create a dictionary with the median abundance in sample_abundance_dict
        using sample_metadata as a reference.
        """
        
        median_sample_abundance = dict()

        for group, samples in sample_metadata.items():
            median_sample_abundance[group] = dict()
            sample_dictionaries = [sample_abundance_dict[sample] for sample in samples]
            genomes = set(list(itertools.chain(*[list(sample_dictionary.keys()) for sample_dictionary in sample_dictionaries])))

            for genome in genomes:
                abundances = [sample_dictionary[genome] for sample_dictionary in sample_dictionaries]
                median_sample_abundance[group][genome] = statistics.median(abundances)

        return median_sample_abundance

    def normalise_by_abundance(self, median_sample_abundances, reaction_abundance_dict, group_to_genome, genome_to_group, genome_groups):

        normalised_abundance_dict = dict()
        for sample_group in list(median_sample_abundances.keys()):
            normalised_abundance_dict[sample_group] = dict()

            for genome_group in genome_groups:
                normalised_abundance_dict[sample_group][genome_group] = dict()

        for sample_group, genome_abundances in median_sample_abundances.items():

            for genome, genome_abundance in genome_abundances.items():

                if(genome in genome_to_group and
                   genome in reaction_abundance_dict):

                    for reaction in list(reaction_abundance_dict[genome].keys()):

                        normalised_value = reaction_abundance_dict[genome][reaction]*genome_abundance

                        genome_group = next(iter(genome_to_group[genome]))

                        if reaction in normalised_abundance_dict[sample_group][genome_group]:
                            normalised_abundance_dict[sample_group][genome_group][reaction].append( normalised_value )
                        else:
                            normalised_abundance_dict[sample_group][genome_group][reaction] = [normalised_value]

        return normalised_abundance_dict

    def average_tpm_by_sample(self, tpm_results, sample_metadata):
        output_dict = dict()
        tpm_dict, annotations, genomes = tpm_results

        for group, samples in sample_metadata.items():
            output_dict[group] = dict()

            for sample in samples:

                for annotation in annotations:

                    if str.encode(sample) in tpm_dict:

                        for genome in genomes:

                            if genome not in output_dict[group]:
                                output_dict[group][genome] = dict()

                            if annotation not in output_dict[group][genome]:
                                output_dict[group][genome][annotation] = list()

                            if genome in tpm_dict[str.encode(sample)]:

                                if annotation in tpm_dict[str.encode(sample)][genome]:
                                    output_dict[group][genome][annotation].append(tpm_dict[str.encode(sample)][genome][annotation])
                                else:
                                    output_dict[group][genome][annotation].append(0.0)

                            else:
                                output_dict[group][genome][annotation].append(0.0)

            for genome, values in output_dict[group].items():

                for annotation in values:
                    output_dict[group][genome][annotation] = sum(output_dict[group][genome][annotation])/len(output_dict[group][genome][annotation])

        return output_dict

    def average_tpm_values(self, transriptome_abundance_dict, group_metadata):
        output_dict = dict()
        reactions = list(self.reactions.keys())
        
        for genome_group_name, group_reaction_abundance_dict in transriptome_abundance_dict.items():
            output_dict[genome_group_name] = dict()

            for group, members in group_metadata.items():
                output_dict[genome_group_name][group] = dict()

                for reaction in reactions:
                    to_average = list()

                    for member in members:

                        if member in group_reaction_abundance_dict:

                            if str.encode(reaction) in group_reaction_abundance_dict[member]:
                                to_average.append(group_reaction_abundance_dict[member][str.encode(reaction)])
                            else:
                                to_average.append(0.0)
                        else:
                            to_average.append(0.0)

                    average_value = sum(to_average) / len(to_average)
                    output_dict[genome_group_name][group][reaction] = average_value

        return output_dict

    def aggregate_dictionary(self, reference_dict, matrix_dict):

        output_dict_mean   = dict()

        for sample, ko_abundances in matrix_dict.items():
            output_dict_mean[sample]   = dict()

            for reaction, ko_list in reference_dict.items():
                abundances = list()

                for ko in ko_list:

                    if ko in ko_abundances:
                        if ko_abundances[ko]>0:
                            abundances.append(ko_abundances[ko])

                    else:
                        logging.debug("ID not found in input matrix: %s" % ko)

                if any(abundances):
                    abundance_mean = sum(abundances)/len(abundances) # average of the abundances...

                else:
                    abundance_mean = 0

                output_dict_mean[sample][reaction] = abundance_mean

        return output_dict_mean

    def mock_metadata(self, genomes):
        genome_to_group = {genome:set([genome]) for genome in genomes}
        genome_groups = set(genomes)
        group_to_genome = dict(genome_to_group) # Make a copy here
        return genome_to_group, genome_groups, group_to_genome

    def network_pipeline(self,
           subparser_name,
           matrix, genome_metadata_path,
           transcriptome_abundances_path, transcriptome_metadata_path,
           metagenome_abundances, metagenome_metadata_path,
           metabolome,
           enrichment_output,
           depth, filter, limit, queries, output_directory):
        '''
        Parameters
        ----------
        matrix
        transcriptome_abundances_path
        metagenome_abundances
        metagenome_metadata_path
        metabolome
        enrichment_output
        depth
        filter
        limit
        queries
        output_directory
        '''
        orthology_matrix, genome_names, _ = Parser.parse_simple_matrix(matrix)
        if genome_metadata_path:
            genome_to_group, genome_groups, group_to_genome = \
                    Parser.parse_metadata_matrix(genome_metadata_path)
        else:
            genome_to_group, genome_groups, group_to_genome = \
                    self.mock_metadata(genome_names)

        reaction_matrix = self.aggregate_dictionary(self.reaction_to_ko, orthology_matrix)

        # Read in fisher results
        if enrichment_output:
            logging.info('Parsing input enrichment results')
            fisher_results = Parser.parse_enrichment_output(enrichment_output)
        else:
            logging.info('No enrichment results provided')
            fisher_results = None

        # Read in metabolome abundances
        if metabolome:
            logging.info('Parsing metabolome abundances')
            abundances_metabolome = Parser.parse_simple_matrix(metabolome)
        else:
            logging.info('No metabolome abundances provided')
            abundances_metabolome = None

        # Read in genome metagenome_abundances
        if metagenome_abundances:
            logging.info('Parsing input genome abundances')
            sample_abundance = Parser.parse_simple_matrix(metagenome_abundances)[0]
            sample_metadata = Parser.parse_metadata_matrix(metagenome_metadata_path)[2]
        else:
            # FIXME : There's always a better way than faking it.
            logging.info('No genome abundances provided')
            sample_abundance = {'MOCK': {x:1 for x in list(reaction_matrix.keys())} }
            sample_metadata = {"abundance": ['MOCK']}

        median_sample_abundances = self.median_genome_abundance(sample_abundance, sample_metadata)
        normalised_abundance_dict = self.normalise_by_abundance(median_sample_abundances, reaction_matrix, group_to_genome, genome_to_group, genome_groups)
        abundances_metagenome = self.average(normalised_abundance_dict)

        # Read in expression (TPM) values
        if transcriptome_abundances_path:
            logging.info("Parsing detectM TPM abundances")
            transcriptome_metadata = Parser.parse_metadata_matrix(transcriptome_metadata_path)[2]
            transcriptome_abundance_dict = self.average_tpm_by_sample(Parser.parse_tpm_values(transcriptome_abundances_path), transcriptome_metadata)
            transcriptome_abundances = self.average_tpm_values(transcriptome_abundance_dict, group_to_genome)
        else:
            transcriptome_abundances = None

        network_builder = NetworkBuilder(group_to_genome, abundances_metagenome,
                                         transcriptome_abundances, abundances_metabolome, 
                                         fisher_results)

        # Run the subcommand specified
        if subparser_name == self.EXPLORE:
            network_lines, node_metadata = network_builder.query_matrix(queries, depth)
        elif subparser_name == self.PATHWAY:
            network_lines, node_metadata = network_builder.pathway_matrix(limit, filter)

        # Write the outputs
        Writer.write(network_lines, os.path.join(output_directory, self.NETWORK_OUTPUT_FILE))
        Writer.write(node_metadata, os.path.join(output_directory, self.METADATA_OUTPUT_FILE))
        
        logging.info('Finished the %s pipeline' % subparser_name)