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
from enrichm.traverse import NetworkTraverser
from enrichm.databases import Databases
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
    COMPOUND_PREFIX = 'C'    
    MAP_PREFIX      = 'map'
    PATHWAY_PREFIX  = 'rn'
    REACTION_PREFIX = 'R'
    ZERO            = '0.0'

    def __init__(self, metadata_keys):
        
        self.d=Databases()

        self.metadata_keys \
                        = list(metadata_keys)
        self.matrix_header \
                        = ["compound", "reaction", 'type'] 
        self.transcriptome_header \
                        = [key + '_reaction_transcriptome' 
                           for key in self.metadata_keys] + \
                          [key + '_reaction_expression' 
                           for key in self.metadata_keys]
        self.compound_header \
                        = [key + '_compound' 
                           for key in self.metadata_keys]
        self.metadata_header \
                        = ['node', 
                           'description',
                           'type',
                           'module', 
                           'module_descr',
                           'pathway', 
                           'pathway_descr', 
                           'node_type']
        self.query_header \
                        =  ['query',
                            'step']
        self.compound_reaction_index_header \
                        = list()
        self.step_header = ['step']

    def _gather_module(self, key):
        if key in self.d.r2m:
            len_list =  [len(self.d.m2r[x]) for x in self.d.r2m[key]]
            module = self.d.r2m[key][len_list.index(max(len_list))]
            module_description = self.d.m[module]
        else:
            module = 'NA'
            module_description='NA'
        return module, module_description
    
    def _gather_pathway(self, key):
        if key in self.d.r2p:
            len_list =  [len(self.d.p2r[x]) for x in self.d.r2p[key]]
            pathway = self.d.r2p[key][len_list.index(max(len_list))]
            pathway_description = self.d.p[pathway]
        else:
            pathway = 'NA'
            pathway_description='NA'
        return pathway, pathway_description
    
    def _parse_queries(self, queries_file):
        '''
        Parse file with one column, each entry a compound ID from kegg
        Parameters
        ----------
        queries_file    - String. Path to file containing Compound ids 
        
        Output
        ------
        Set. Compound Ids to start from in the explore network
        '''
        return set([x.strip() for x in open(queries_file)])

    def normalise(self, reaction_abundances):
        probability_list = list()
        for reaction_abundance in reaction_abundances:
            probability = float(reaction_abundance)/sum(reaction_abundances)
            probability_list.append(probability)
        return probability_list

    def get_transition(self, probabilities):
        n=random.uniform(0,1)
        x = list()
        for probability in probabilities:
            x.append(probability-n)
        return x.index(min(x, key=abs))

    def traverse_matrix(self,
                        abundances_metagenome,
                        abundances_transcriptome,
                        possible_reactions,
                        query_list,
                        number_of_queries,
                        steps): 

        network_visitations={key:dict() for key in self.metadata_keys}

        for _ in range(number_of_queries):
            previous_reaction=''
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
        for c in set(chain(*possible_reactions.values())):
            output_line = [c]
            for key in self.metadata_keys:
                if c in network_visitations[key]:
                    output_line.append(str(network_visitations[key][c]))
                else:
                    output_line.append('0')    
            output_lines.append('\t'.join(output_line))
        return output_lines

    def all_matrix(self, 
                   abundances, 
                   abundances_metabolome,
                   fisher_results,
                   reference_dict):
        '''
        Parameters
        ----------
        '''
        seen_nodes = set()
        # Construct headers to network matrices
        seen_reactions = set(nested_dict_vals(abundances))

        groups = list(abundances.keys())
        network_lines  = ['\t'.join(self.matrix_header + list(['_'.join([a,b]) for a,b in product(self.metadata_keys, groups)]))]

        if abundances_metabolome:
            node_metadata_lines = ['\t'.join(self.metadata_header + 
                                             self.compound_header)]
        else:     
            node_metadata_lines = ['\t'.join(self.metadata_header)]

        for reaction, entry in reference_dict.items():

            for compound in entry:

                if reaction in seen_reactions:

                    reaction_line = [compound, reaction]
                    if fisher_results:
                        enriched_term = list()
                        
                        for compared_group in list(fisher_results.keys()):

                            if any(set(self.d.r2k[reaction]).intersection(fisher_results[compared_group])):
                                enriched_term.append(compared_group)
                        
                        if len(enriched_term)>0:
                            enriched_term = '_'.join(enriched_term)
                        
                        else:
                            enriched_term = 'NA'
                    else:
                        enriched_term = 'NA'

                    reaction_line.append(enriched_term)

                    for key in self.metadata_keys:

                        for _, group_abundances in abundances.items():
                            if reaction in group_abundances[key]:
                                reaction_line.append(str(group_abundances[key][reaction]))
                            else:
                                reaction_line.append(self.ZERO)

                    output_line = '\t'.join(reaction_line)

                    if sum([float(x) for x in reaction_line[3:]])>0:
                        if output_line not in network_lines:
                            network_lines.append(output_line)
                    #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                    #~#~#~#~#~#~#~#~#~ Fill in node metadata ~#~#~#~#~#~#~#~#~#
                    if compound not in seen_nodes:
                        if compound in self.d.compound_desc_dict:
                            compound_type = \
                                ','.join(self.d.compound_desc_dict[compound]['A'])
                        else:
                            compound_type = 'NA'
                        
                        metadata_list = [compound, self.d.c[compound], compound_type, 
                                         'NA', 'NA', 'NA', 'NA', 'compound']
                        if abundances_metabolome:   
                            # This is fragile - compound abundances need to be aggregated
                            # by the user. TODO: Integrate compounds into kegg_matrix
                            # module.
                            for sample in self.metadata_keys:
                                if compound in abundances_metabolome[sample]:
                                    c_ab = abundances_metabolome[sample][compound]
                                else:
                                    c_ab = 0    
                                metadata_list.append(str(c_ab))
                        node_metadata_lines.append('\t'.join(metadata_list))
                        seen_nodes.add(compound)

                    if reaction not in seen_nodes:
                        module, module_description \
                                        = self._gather_module(reaction)
                        pathway, pathway_description \
                                        = self._gather_pathway(reaction)
                        metadata_list = [reaction, 
                                         self.d.r[reaction],
                                         'NA', 
                                         module, 
                                         module_description, 
                                         pathway,
                                         pathway_description,
                                         'reaction']
                        if abundances_metabolome:   
                            metadata_list += ['-5' for sample in self.metadata_keys]
                        node_metadata_lines.append('\t'.join(metadata_list))
                        seen_nodes.add(reaction)

        return network_lines, node_metadata_lines

    def query_matrix(self, 
                     abundances_metagenome, 
                     abundances_transcriptome,
                     abundances_expression,
                     queries, 
                     depth):
        '''
        Parameters
        ----------
        '''
        steps         = 0
        queries_list  = self._parse_queries(queries)
        seen_steps    = set()
        seen_nodes    = set()
        level_queries = set()

        if(abundances_transcriptome and abundances_expression):
            network_lines  \
                    = ['\t'.join(self.matrix_header + 
                                 self.transcriptome_header +
                                 self.step_header)]
        
        else:
            network_lines  \
                    = ['\t'.join(self.matrix_header+
                                 self.step_header)] 
        node_metadata_lines \
                    = ['\t'.join(self.metadata_header + 
                                 self.query_header )]
        
        while depth>0:
            if any(level_queries):
                queries_list = set(level_queries)
                level_queries = set()
            for reaction, entry in self.d.r2c.items():
                if any([reaction in x.keys() 
                        for x in abundances_metagenome.values()]):

                    if any(queries_list.intersection(entry)):
                        reaction_compounds = [x for x in entry 
                                              if x not in NetworkTraverser.to_omit]
                        for compound in reaction_compounds:
                            
                            #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                            #~#~#~#~#~#~#~#~ Fill in abundances ~#~#~#~#~#~#~#~
                            reaction_line = [compound, reaction]
                        
                            for key in self.metadata_keys:
                                if reaction in abundances_metagenome[key]:
                                    reaction_line.append(str(abundances_metagenome[key][reaction]))
                                else:
                                    reaction_line.append(self.ZERO)

                            if(abundances_transcriptome and abundances_expression):
                                for key in self.metadata_keys:
                                    if reaction in abundances_transcriptome[key]:
                                        reaction_line.append(str(abundances_transcriptome[key][reaction]))
                                    else:
                                        reaction_line.append(self.ZERO)

                                for key in self.metadata_keys:
                                    if reaction in abundances_expression[key]:
                                        reaction_line.append(str(abundances_expression[key][reaction]))
                                    else:
                                        reaction_line.append(self.ZERO)

                            output_line = '\t'.join(reaction_line)
                            if output_line not in seen_steps:
                                seen_steps.add(output_line)
                                network_lines.append(output_line+'\t%i' % steps)
                                
                            #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                            #~#~#~#~#~#~#~ Fill in node metadata  ~#~#~#~#~#~#~
                            index = str((steps if compound 
                                         in queries_list else steps+1))
                            
                            if compound not in seen_nodes:
                                #if compound in query_list:
                                #    query_ab = [str(query_list[compound][key]) 
                                #                for key in self.metadata_keys]
                                #else:
                                #    query_ab = ['NA' for key 
                                #                in self.metadata_keys]
                                
                                if compound in self.d.compound_desc_dict:
                                    compound_type = \
                                        ','.join(self.d.compound_desc_dict[compound]['A'])
                                else:
                                    compound_type = 'NA'
                                    
                                if compound in queries_list:
                                    is_query = 'True'
                                else:
                                    is_query = 'False'
                                    
                                node_metadata_lines.append(
                                    '\t'.join([compound, 
                                               self.d.c[compound],
                                               compound_type, 
                                               'NA',
                                               'NA', 
                                               'NA',
                                               'NA',  
                                               'compound'] +
                                               #query_ab +
                                               [is_query, index])
                                                           )
                                
                                seen_nodes.add(compound)
                                if compound not in queries_list:
                                    level_queries.add(compound)
                            if reaction not in seen_nodes:
                                module, module_description \
                                                = self._gather_module(reaction)
                                pathway, pathway_description \
                                                = self._gather_pathway(reaction)
                                node_metadata_lines.append(
                                    '\t'.join([reaction, 
                                               self.d.r[reaction],
                                               'NA', 
                                               module, 
                                               module_description, 
                                               pathway,
                                               pathway_description,
                                               'reaction',
                                               'False',
                                               index])
                                                           )
                                seen_nodes.add(reaction)

            steps+=1
            depth-=1
            
            logging.info("Step %i complete with %i queries to continue with" \
                                              % (steps, len(level_queries)))
        return network_lines, node_metadata_lines
    
    def pathway_matrix(self,
                       abundances_metagenome,
                       abundances_metabolome,
                       fisher_results,
                       limit,
                       filter):

        possible_reactions = set()

        if any(limit):
            
            for entry in limit:
            
                if(entry.startswith(self.MAP_PREFIX) or 
                   entry.startswith(self.PATHWAY_PREFIX)): 
            
                    for reaction in self.d.p2r[entry]:
                        possible_reactions.add(reaction)
            
                elif(entry.startswith(self.MODULE_PREFIX)):
            
                    for reaction in self.d.m2r[entry]:
                        possible_reactions.add(reaction)
            
                elif(entry.startswith(self.REACTION_PREFIX)):
                    possible_reactions.add(entry)
            possible_reactions = {reaction:self.d.r2c[reaction] 
                                  for reaction in
                                  possible_reactions
                                  if reaction in self.d.r2c}
        else:
            possible_reactions = self.d.r2c
        
        for entry in filter:

            if entry in possible_reactions:
            
                del possible_reactions[entry]

        network_lines, node_metadata_lines = \
            self.all_matrix(abundances_metagenome,
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
                    for reaction in self.d.p2r[entry]:
                        possible_reactions.add(reaction)
                elif(entry.startswith(self.MODULE_PREFIX)):
                    for reaction in self.d.m2r[entry]:
                        possible_reactions.add(reaction)
                elif(entry.startswith(self.REACTION_PREFIX)):
                    possible_reactions.add(entry)

            possible_reactions = {reaction:self.d.r2c[reaction] 
                                  for reaction in
                                  possible_reactions}                
        else:
            possible_reactions = self.d.r2c
        
        for entry in filter:
            if entry in possible_reactions:
                del possible_reactions[entry]
        
        possible_compounds = set(chain(*possible_reactions.values()))

        if len(starting_compounds)==0:
            query_list=[x for x in self.d.c.keys()
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
