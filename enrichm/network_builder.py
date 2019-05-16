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
import random
import re
from itertools import chain, product
# Local
from enrichm.databases import Databases
from enrichm.parser import Parser
###############################################################################

def nested_dict_vals(d):
    reaction_regex = '(R\d{5})$'

    for key, item in d.items():

        if isinstance(item, dict):
            yield from nested_dict_vals(item)
        
        else:
            if type(reaction_regex) == str:
                if re.match(reaction_regex, key):
                    yield key
            elif type(reaction_regex) == bytes:
                key = key.decode()
                if re.match(reaction_regex, key):
                    yield key

            else:
                raise Exception("Invalid key in in nested dict!")

class NetworkBuilder:
       
    MODULE_PREFIX   = 'M'
    MAP_PREFIX      = 'map'
    PATHWAY_PREFIX  = 'rn'
    REACTION_PREFIX = 'R'
    ZERO            = '0.0'

    def __init__(self, metadata):
        
        self.databases=Databases()
        self.metadata_keys = list(metadata.keys())
        self.matrix_header = ["compound", "reaction", 'type'] 
        self.transcriptome_header = [key + '_reaction_transcriptome' for key in self.metadata_keys] 
        self.compound_header = [key + '_compound' for key in self.metadata_keys]
        self.metadata_header = ['node', 'description', 'type', 'module', 'module_descr', 'pathway', 'pathway_descr', 'node_type']
        self.query_header =  ['query', 'step']
        self.step_header = ['step']
        self.to_omit = set(["C00828",  # Menaquinone
                    "C00534",  # Pyridoxamine
                    "C00006",  # NADP+
                    "C00003",  # NAD+
                    "C00002",  # ATP
                    "C00314",  # Pyridoxine
                    "C00864",  # Pantothenate
                    "C00504",  # Folate
                    "C00032",  # Heme
                    "C05443",  # Vitamin D3
                    "C00253",  # Nicotinate
                    "C00250",  # Pyridoxal
                    "C11378",  # Ubiquinone-10
                    "C05777",  # Coenzyme F430
                    "C00072",  # Ascorbate
                    "C00378",  # Thiamine
                    "C00101",  # Tetrahydrofolate
                    "C00029",  # UDP-glucose
                    "C00068",  # Thiamin diphosphate
                    "C00061",  # FMN
                    "C00063",  # CTP
                    "C05776",  # Vitamin B12
                    "C00113",  # PQQ
                    "C18237",  # Molybdoenzyme molybdenum cofactor
                    "C00051",  # Glutathione
                    "C00010",  # CoA
                    "C00016",  # FAD
                    "C00018",  # Pyridoxal phosphate
                    "C00019",  # S-Adenosyl-L-methionine
                    "C00153",  # Nicotinamide
                    "C04628",  # Coenzyme B
                    "C00862",  # Methanofuran
                    "C15672",  # Heme O
                    "C15670",  # Heme A
                    "C02059",  # Phylloquinone
                    "C03576",  # Coenzyme M
                    "C05441",  # Vitamin D2
                    "C00272",  # Tetrahydrobiopterin
                    "C02477",  # alpha-Tocopherol
                    "C00473",  # Retinol
                    "C00120",  # Biotin
                    "C00725",  # Lipoate
                    "C00053",  # 3'-Phosphoadenylyl sulfate
                    "C00194",  # Cobamide coenzyme
                    "C00255",  # Riboflavin
                    'C00001',  # H2O
                    'C00008',  # ADP
                    'C00013',  # Diphosphate
                    'C00004',  # NADH
                    'C00005',  # NADPH
                    'C00080',  # H+
                    'C00009',  # Orthophosphate
                    'C00008',  # ADP
                    'C00004',  # NADH
                    'C00020',  # AMP
                    'C00007',  # Oxygen
                    'C00015']) # UDP

    def _gather_module(self, key):
        
        if key in self.databases.r2m:
            len_list =  [len(self.databases.m2r[x]) for x in self.databases.r2m[key]]
            module = self.databases.r2m[key][len_list.index(max(len_list))]
            module_description = self.databases.m[module]
        else:
            module = 'NA'
            module_description='NA'
        
        return module, module_description
    
    def _gather_pathway(self, key):
        
        if key in self.databases.r2p:
            len_list =  [len(self.databases.p2r[x]) for x in self.databases.r2p[key]]
            pathway = self.databases.r2p[key][len_list.index(max(len_list))]
            pathway_description = self.databases.p[pathway]
        else:
            pathway = 'NA'
            pathway_description='NA'
        
        return pathway, pathway_description
    
    def traverse_matrix(self,
                        abundances_metagenome,
                        abundances_transcriptome,
                        possible_reactions,
                        query_list,
                        number_of_queries,
                        steps): 

        network_visitations={key:dict() for key in self.metadata_keys}

        for _ in range(number_of_queries):
            previous_reaction=None
            
            for _ in range(steps):

                for key in self.metadata_keys:
                    reactions = list()
                    possible_reaction_abundance = list()
                    possible_reaction_name = list()
                    
                    for reaction in reactions:  
                    
                        if(reaction in abundances_metagenome[key] and
                           reaction in abundances_transcriptome[key]):
                    
                            if reaction!=previous_reaction:
                    
                                if abundances_metagenome[key][reaction]>0:
                                    possible_reaction_name.append(reaction)
                                    possible_reaction_abundance.append(abundances_transcriptome[key][reaction])

                    
        output_lines = ['\t'.join(['C'] + self.metadata_keys)]
        
        for possible_reaction in set(chain(*possible_reactions.values())):
            output_line = [possible_reaction]
            
            for key in self.metadata_keys:
            
                if possible_reaction in network_visitations[key]:
                    output_line.append(str(network_visitations[key][possible_reaction]))
                else:
                    output_line.append('0')    
            
            output_lines.append('\t'.join(output_line))
        
        return output_lines


    def gather_compound_metadata(self, compound, abundances_metabolome):
        if compound in self.databases.compound_desc_dict:
            compound_type = ','.join(self.databases.compound_desc_dict[compound]['A'])
        else:
            compound_type = 'NA'
        
        if abundances_metabolome:   

            for sample in self.metadata_keys:
                
                if compound in abundances_metabolome[sample]:
                    c_ab = abundances_metabolome[sample][compound]
                else:
                    c_ab = 0    
                        
        compound_metadata_list = [compound, self.databases.c[compound], compound_type, 'NA',
                                  'NA', 'NA', 'NA', 'compound', str(c_ab)]
        
        return compound_metadata_list
    
    def gather_reaction_metadata(self, reaction, abundances_metabolome):
        
        module, module_description = self._gather_module(reaction)
        pathway, pathway_description = self._gather_pathway(reaction)
        reaction_metadata_list = [reaction, self.databases.r[reaction], 'NA', module, 
                         module_description, pathway, pathway_description, 'reaction']
        if abundances_metabolome:   
            reaction_metadata_list += ['-5' for sample in self.metadata_keys]
        return reaction_metadata_list
    
    def gather_enriched_term(self, reaction, fisher_results):
        enriched_term = list()
        
        for compared_group in list(fisher_results.keys()):

            if any(set(self.databases.r2k[reaction]).intersection(fisher_results[compared_group])):
                enriched_term.append(compared_group)
        
        if len(enriched_term)>0:
            enriched_term = '_'.join(enriched_term)
        
        else:
            enriched_term = 'NA'
        
        return enriched_term

    def get_reaction_compounds(self, compounds_list):
        
        for compound in compounds_list:
        
            if compound not in self.to_omit:
                yield compound
    
    def gather_reaction_edge_data(self, compound, reaction, fisher_results, abundances_metagenome, abundances_transcriptome):
        
        if fisher_results:
            enriched_term = self.gather_enriched_term(reaction, fisher_results)
        else:
            enriched_term = 'NA'

        reaction_line = [compound, reaction, enriched_term]
        
        for key in self.metadata_keys:
            # FIXME: What am abundances_metagenome
            for _, group_abundances in abundances_metagenome.items():
                if reaction in abundances_metagenome[key]:
                    reaction_line.append(str(abundances_metagenome[key][reaction]))
                else:
                    reaction_line.append(self.ZERO)

        if abundances_transcriptome:

            for key in self.metadata_keys:

                if reaction in abundances_transcriptome[key]:
                    reaction_line.append(str(abundances_transcriptome[key][reaction]))
                else:
                    reaction_line.append(self.ZERO)
        if sum([float(x) for x in reaction_line[3:]])>0:
            return reaction_line

    def all_matrix(self,
                   abundances_metagenome,
                   abundances_transcriptome,
                   abundances_metabolome,
                   fisher_results,
                   reference_dict):
        '''
        Parameters
        ----------
        '''
        seen_nodes = set()
        # Construct headers to network matrices
        seen_reactions = set(nested_dict_vals(abundances_metagenome))
        
        groups = list(abundances.keys())
        network_lines = [self.matrix_header + list(['_'.join([a,b]) for a,b in product(self.metadata_keys, groups)])]
        
        if abundances_metabolome:
            node_metadata_lines = [self.metadata_header + self.compound_header]
        else:     
            node_metadata_lines = [self.metadata_header]

        for reaction, entry in reference_dict.items():

            for compound in entry:
                
                if reaction in seen_reactions:
                    
                    #~#~#~#~#~#~#~#~ Fill in abundances_metagenome ~#~#~#~#~#~#~#~
                    reaction_line = self.gather_reaction_edge_data(compound, reaction, fisher_results, abundances_metagenome, abundances_transcriptome)
                    
                    if reaction_line:
                        network_lines.append(reaction_line)
                            
                    #~#~#~#~#~#~#~#~#~ Fill in node metadata ~#~#~#~#~#~#~#~#~#
                    if compound not in seen_nodes:
                        compound_metadata_list = self.gather_compound_metadata(compound, abundances_metabolome)

                        node_metadata_lines.append(compound_metadata_list)
                        seen_nodes.add(compound)

                    if reaction not in seen_nodes:
                        reaction_metadata_list = self.gather_reaction_metadata(reaction, abundances_metabolome)
                        
                        node_metadata_lines.append(reaction_metadata_list)
                        seen_nodes.add(reaction)

        return network_lines, node_metadata_lines
    
    def query_matrix(self, 
                     abundances_metagenome, 
                     abundances_transcriptome,
                     abundances_metabolome,
                     fisher_results,
                     queries, 
                     depth):
        '''
        Parameters
        ----------
        '''
        steps         = 0
        queries_list  = Parser.parse_single_column_text_file(queries)
        seen_reactions = set(nested_dict_vals(abundances_metagenome))
        seen_nodes    = set()
        level_queries = set()
        
        node_metadata_lines = [self.metadata_header + self.query_header]
        
        if abundances_metabolome:
            node_metadata_lines[0] += self.compound_header
        
        if abundances_transcriptome:
            network_lines = [self.matrix_header + self.transcriptome_header + self.step_header]
        else:
            network_lines = [self.matrix_header + self.step_header] 
        
        while depth>0:
            
            if any(level_queries):
                queries_list = set(level_queries)
                level_queries = set()
            
            for reaction, entry in self.databases.r2c.items():
                
                if reaction in seen_reactions:

                    if any(queries_list.intersection(entry)):

                        for compound in self.get_reaction_compounds(entry):
                            
                            #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                            #~#~#~#~#~#~#~#~ Fill in abundances ~#~#~#~#~#~#~#~
                            reaction_line = self.gather_reaction_edge_data(compound, reaction, fisher_results, abundances_metagenome, abundances_transcriptome)
                            
                            if reaction_line:
                                reaction_line.append(str(steps))
                                network_lines.append(reaction_line)

                            #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                            #~#~#~#~#~#~#~ Fill in node metadata  ~#~#~#~#~#~#~
                            index = str((steps if compound in queries_list else steps+1))
                            
                            if compound not in seen_nodes:
                                
                                metadata_list = self.gather_compound_metadata(compound, abundances_metabolome)
                                    
                                if compound in queries_list:
                                    is_query = 'True'
                                else:
                                    is_query = 'False'
                                    
                                node_metadata_lines.append(metadata_list + [is_query, index])
                                seen_nodes.add(compound)
                                
                                if compound not in queries_list:
                                    level_queries.add(compound)
                            
                            if reaction not in seen_nodes:
                                metadata_list = self.gather_reaction_metadata(reaction, abundances_metabolome)
                                
                                node_metadata_lines.append(metadata_list + ['False', index])
                                seen_nodes.add(reaction)

            steps+=1
            depth-=1
            
            logging.info("Step %i complete with %i queries to continue with" % (steps, len(level_queries)))
        
        return network_lines, node_metadata_lines
    
    def pathway_matrix(self, abundances_metagenome, abundances_transcriptome, abundances_metabolome, fisher_results, limit, filter):

        possible_reactions = set()

        if any(limit):
            
            for entry in limit:
            
                if(entry.startswith(self.MAP_PREFIX) or 
                   entry.startswith(self.PATHWAY_PREFIX)): 
            
                    for reaction in self.databases.p2r[entry]:
                        possible_reactions.add(reaction)
            
                elif(entry.startswith(self.MODULE_PREFIX)):
            
                    for reaction in self.databases.m2r[entry]:
                        possible_reactions.add(reaction)
            
                elif(entry.startswith(self.REACTION_PREFIX)):
                    possible_reactions.add(entry)
            possible_reactions = {reaction:self.databases.r2c[reaction] 
                                  for reaction in
                                  possible_reactions
                                  if reaction in self.databases.r2c}
        else:
            possible_reactions = self.databases.r2c
        
        for entry in filter:

            if entry in possible_reactions:
            
                del possible_reactions[entry]

        network_lines, node_metadata_lines = \
            self.all_matrix(abundances_metagenome,
                            abundances_transcriptome,
                            abundances_metabolome,
                            fisher_results,
                            possible_reactions)

        return network_lines, node_metadata_lines
  
    def traverse(self, abundances_metagenome, abundances_transcriptome, limit, filter, 
                 starting_compounds, steps, number_of_queries):

        if any(limit):

            for entry in limit:
            
                if(entry.startswith(self.MAP_PREFIX) or 
                   entry.startswith(self.PATHWAY_PREFIX)): 
            
                    for reaction in self.databases.p2r[entry]:
                        possible_reactions.add(reaction)
            
                elif(entry.startswith(self.MODULE_PREFIX)):
            
                    for reaction in self.databases.m2r[entry]:
                        possible_reactions.add(reaction)
            
                elif(entry.startswith(self.REACTION_PREFIX)):
                    possible_reactions.add(entry)

            possible_reactions = {reaction:self.databases.r2c[reaction] for reaction in possible_reactions}
        
        else:
            possible_reactions = self.databases.r2c
        
        for entry in filter:
        
            if entry in possible_reactions:
                del possible_reactions[entry]
        
        possible_compounds = set(chain(*possible_reactions.values()))

        if len(starting_compounds)==0:
            query_list=[x for x in self.databases.c.keys()
                        if x in possible_compounds]
        else:
            query_list=[x for x in starting_compounds 
                        if x in possible_compounds]

        output_lines = self.traverse_matrix(abundances_metagenome, 
                                            abundances_transcriptome, 
                                            possible_reactions, 
                                            query_list,
                                            number_of_queries, 
                                            steps)
        return output_lines
