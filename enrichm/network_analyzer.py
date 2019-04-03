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

import logging
import os
import itertools
from enrichm.network_builder import NetworkBuilder
from enrichm.parser import Parser
from enrichm.databases import Databases
from enrichm.kegg_matrix import KeggMatrix

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.7"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

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

    NETWORK_OUTPUT_FILE  = 'network.tsv'
    METADATA_OUTPUT_FILE = 'metadata.tsv'    
    TRAVERSE_OUTPUT_FILE = 'traverse.tsv'    

    def __init__(self, metadata):
        self.metadata = dict()
        
        for line in open(metadata):
            
            if line.startswith('#'): 
                continue
            
            split_line = line.strip().split('\t')
            
            if len(split_line) == 1:
                raise Exception("Only one column detected in metadata file, please check that your file is tab separated")
            else:
                sample_id, group = split_line 

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
    
    def _average(self, d):
            
        for sample_group, group_dict in d.items():

            for group, reaction_dict in group_dict.items():
            
                for reaction, value in reaction_dict.items():
                    d[sample_group][group][reaction] = sum(value) / len(value)
        
        return d

    def normalise_by_abundance(self, sample_abundance_dict, sample_metadata, reaction_abundance_dict, metadata):

        reversed_metadata = dict()
        for key, items in metadata.items():

            for item in items:
                reversed_metadata[item] = key
        
        averaged_sample_abundance = dict()
        
        for group, samples in sample_metadata.items():
            averaged_sample_abundance[group] = dict()
            dicts = [sample_abundance_dict[x] for x in samples]
            genomes = set(list(itertools.chain(*[list(x.keys()) for x in dicts])))
            
            for genome in genomes:
                abundances = [d[genome] for d in dicts]
                average = sum(abundances) / len(abundances)
                averaged_sample_abundance[group][genome] = average # Median might be better?
        new_dict = {x:dict() for x in list(averaged_sample_abundance.keys())}

        for sample_group, genome_abundances in averaged_sample_abundance.items():
            
            for genome, genome_abundance in genome_abundances.items():

                if genome in reversed_metadata:

                    if genome in reaction_abundance_dict[sample_group]:

                        for reaction in list(reaction_abundance_dict[sample_group][genome].keys()):
                            
                            if reaction in reaction_abundance_dict[sample_group][genome]:
                                normalised_value = reaction_abundance_dict[sample_group][genome][reaction]*genome_abundance             

                            else:
                                normalised_value = 0.0
                            group = reversed_metadata[genome]

                            if group not in new_dict[sample_group]:
                                new_dict[sample_group][group] = dict()
                            
                            if reaction not in new_dict[sample_group][group]:
                                new_dict[sample_group][group][reaction] = [0.0]

                            new_dict[sample_group][group][reaction].append( normalised_value )
        
        new_dict = self._average(new_dict) # taking averages here again, might be better accumulated?

        return new_dict

    def _parse_enrichment_output(self, enrichment_output):
        fisher_results = dict()
        
        for file in os.listdir(enrichment_output):

            if file.endswith("fisher.tsv"):
                file = os.path.join(enrichment_output, file)
                file_io = open(file)
                file_io.readline()

                for line in file_io:
                    split_line = line.strip().split('\t')
                     
                    if len(fisher_results) == 0:
                        fisher_results[split_line[1]] = list()
                        fisher_results[split_line[2]] = list()

                    if float(split_line[-2])<0.05:
                        g1_t = float(split_line[3])
                        g1_f = float(split_line[4])
                        g2_t = float(split_line[5])
                        g2_f = float(split_line[6])

                        if g1_t == 0:
                            fisher_results[split_line[2]].append( split_line[0] )

                        elif g2_t == 0:
                            fisher_results[split_line[1]].append( split_line[0] )

                        elif ( ((g1_t/(g1_t+g1_f))) / ((g2_t/(g2_t+g2_f))) )>1:
                            fisher_results[split_line[1]].append( split_line[0] )

                        else:
                            fisher_results[split_line[2]].append( split_line[0] )
                
        if len(fisher_results.keys())>0:
            return fisher_results
        
        else:
            raise Exception("Malformatted enrichment output: %s" % enrichment_output)
    
    def average_tpm_by_sample(self, tpm_results, sample_metadata, group_metadata):
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
        
        new_output_dict = dict()
        reactions = list(Databases().r.keys())
        for key, item in output_dict.items():
            new_output_dict[key] = dict()

            for group, members in group_metadata.items():
                new_output_dict[key][group] = dict()

                for reaction in reactions:
                    to_average = list()
                    for member in members:
                        if member in item:
                            if str.encode(reaction) in item[member]:
                                to_average.append(item[member][str.encode(reaction)])
                            else:
                                to_average.append(0.0)
                    average_value = sum(to_average) / len(to_average)
                    new_output_dict[key][group][reaction] = average_value

        return new_output_dict

    def do(self,
            matrix,
            tpm_values,
            abundance,
            abundance_metadata,
            metabolome,
            enrichment_output,
            depth,
            filter,
            limit,
            queries,
            starting_compounds,
            steps,
            number_of_queries,
            output_directory):
        '''
        Parameters
        ----------
        matrix
        tpm_values
        abundance
        abundance_metadata
        metabolome
        enrichment_output
        depth
        filter
        limit
        queries
        starting_compounds
        steps
        number_of_queries
        output_directory
        '''

        km = KeggMatrix(matrix)
        nb = NetworkBuilder(self.metadata.keys())

        # Read in fisher results
        if enrichment_output:
            fisher_results = self._parse_enrichment_output(enrichment_output)
        else:
            fisher_results = None    

        # Read in genome abundance
        if abundance:
            sample_abundance = km._parse_matrix(abundance)
            sample_metadata = list(km._parse_matrix(abundance_metadata).values())[0]
            d = dict()
            
            for key, item in sample_metadata.items():
                
                if item not in d:
                    d[item] = list()
                
                d[item].append(key)
            
            sample_metadata = d
        
        else:
            sample_abundance = {'MOCK': {x:1 for x in list(km.reaction_matrix.keys())} }
            sample_metadata = {"a": ['MOCK']}
            km.reaction_matrix = {"a": km.reaction_matrix}
        
        # Read in expression (TPM) values
        if tpm_values:
            logging.info("Parsing detectM TPM values")
            normalised_abundances \
                    = self.average_tpm_by_sample(Parser.parse_tpm_values(tpm_values), sample_metadata, self.metadata)

        else:
            normalised_abundances \
                = self.normalise_by_abundance(sample_abundance,
                                              sample_metadata,
                                              km.reaction_matrix,
                                              self.metadata)
        if metabolome:
            abundances_metabolome = km._parse_matrix(metabolome)
        
        else:
            abundances_metabolome = None

        logging.info('Generating pathway network')

        network_lines, node_metadata = \
                        nb.pathway_matrix(normalised_abundances, 
                                          abundances_metabolome,
                                           fisher_results,
                                           limit,
                                           filter)

        self._write_results(os.path.join(output_directory, self.NETWORK_OUTPUT_FILE), network_lines)
        self._write_results(os.path.join(output_directory, self.METADATA_OUTPUT_FILE), node_metadata)
