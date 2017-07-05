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
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
 
###############################################################################

import logging
import pickle
import os
import random
import itertools

from traverse import NetworkTraverser

###############################################################################

class NetworkBuilder:
    
    DATA_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             '..', 
                             'data')
    
    VERSION = os.path.join(DATA_PATH, 'VERSION')
    PICKLE  = 'pickle'

    COMPOUND_DESC_PICKLE = os.path.join(DATA_PATH, 'br08001')    
    R2RPAIR = os.path.join(DATA_PATH, 'reaction_to_rpair')
    R2K     = os.path.join(DATA_PATH, 'reaction_to_orthology')
    R2C = os.path.join(DATA_PATH, 'reaction_to_compound')
    R2M = os.path.join(DATA_PATH, 'reaction_to_module')
    M2R = os.path.join(DATA_PATH, 'module_to_reaction')
    R2P = os.path.join(DATA_PATH, 'reaction_to_pathway')
    P2R = os.path.join(DATA_PATH, 'pathway_to_reaction')
    C2R = os.path.join(DATA_PATH, 'compound_to_reaction')
    C   = os.path.join(DATA_PATH, 'compound_descriptions')    
    R   = os.path.join(DATA_PATH, 'reaction_descriptions')
    P   = os.path.join(DATA_PATH, 'pathway_descriptions')
    M   = os.path.join(DATA_PATH, 'module_descriptions')
    
    MODULE_PREFIX   = 'M'
    COMPOUND_PREFIX = 'C'    
    MAP_PREFIX      = 'map'
    PATHWAY_PREFIX  = 'rn'
    REACTION_PREFIX = 'R'
    ZERO            = '0.0'

    def __init__(self, metadata_keys):
        
        self.VERSION = open(self.VERSION).readline().strip()
        
        logging.info("Loading reaction to pathway information")
        self.r2p = pickle.load(open('.'.join([self.R2P, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading pathway to reaction information")
        self.p2r = pickle.load(open('.'.join([self.P2R, 
                                                 self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading reaction to orthology information")
        self.r2k = pickle.load(open('.'.join([self.R2K, self.VERSION, 
                                              self.PICKLE])))
        logging.info("Done")
        logging.info("Loading reaction to module information")
        self.r2m = pickle.load(open('.'.join([self.R2M, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading reaction to module information")
        self.m2r = pickle.load(open('.'.join([self.M2R, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading reaction to compound information")
        self.r2c = pickle.load(open('.'.join([self.R2C, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading reaction to rpair information")
        self.r2rpair = pickle.load(open('.'.join([self.R2RPAIR, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading compound to reaction information")
        self.c2r = pickle.load(open('.'.join([self.C2R,
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading compound descriptions")
        self.c   = pickle.load(open('.'.join([self.C, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading pathway descriptions")
        self.p   = pickle.load(open('.'.join([self.P, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading reaction descriptions")
        self.r   = pickle.load(open('.'.join([self.R, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading module descriptions")
        self.m   = pickle.load(open('.'.join([self.M, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading compound classifications")
        self.compound_desc_dict \
                 = pickle.load(open('.'.join([self.COMPOUND_DESC_PICKLE, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")

        self.metadata_keys \
                        = metadata_keys
        self.matrix_header \
                        = ["compound", "reaction"] + \
                          [key + '_reaction_metagenome' 
                           for key in self.metadata_keys]
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
                        = []
    def _gather_module(self, key):
        if key in self.r2m:
            len_list =  [len(self.m2r[x]) for x in self.r2m[key]]
            module = self.r2m[key][len_list.index(max(len_list))]
            module_description = self.m[module]
        else:
            module = 'NA'
            module_description='NA'
        return module, module_description
    
    def _gather_pathway(self, key):
        if key in self.r2p:
            len_list =  [len(self.p2r[x]) for x in self.r2p[key]]
            pathway = self.r2p[key][len_list.index(max(len_list))]
            pathway_description = self.p[pathway]
        else:
            pathway = 'NA'
            pathway_description='NA'
        return pathway, pathway_description
    
    def _parse_queries(self, queries):
        output_dict = {}
        queries_io  = open(queries)
        header      = queries_io.readline().strip().split('\t')[1:]
        
        for line in queries_io:
            sline = line.strip().split()
            compound = sline[0]
            output_dict[compound] = {}
            
            if any(header):
                for idx, group in enumerate(header):
                    output_dict[compound][group] = sline[idx+1]
                
            else:
                for idx, group in enumerate(self.metadata_keys):
                    output_dict[compound][group] = 'NA'
        
        return output_dict


    def normalise(self, reaction_abundances):
        probability_list = []
        for reaction_abundance in reaction_abundances:
            probability = float(reaction_abundance)/sum(reaction_abundances)
            probability_list.append(probability)
        return probability_list

    def get_transition(self, probabilities):
        n=random.uniform(0,1)
        x = []
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

        network_visitations={key:{} for key in self.metadata_keys}
        for _ in range(number_of_queries):
            starting_compound = random.choice(query_list)
            previous_reaction=''
            next_compound=''
            for i in range(steps):
                if i==0:
                    iter_batch={key:starting_compound for key in self.metadata_keys}

                for key in self.metadata_keys:
                    compound = iter_batch[key]
                    
                    try:
                        reactions = [x for x in self.c2r[compound]
                                     if x in self.r2rpair]
                    except:
                        reactions=[]
                    possible_reaction_abundance = []
                    possible_reaction_name = []
                    for reaction in reactions:  
                        if(reaction in abundances_metagenome[key] and
                           reaction in abundances_transcriptome[key]):
                            if reaction!=previous_reaction:
                                if abundances_metagenome[key][reaction]>0:
                                    possible_reaction_name.append(reaction)
                                    possible_reaction_abundance.append(abundances_transcriptome[key][reaction])

                    if any(possible_reaction_abundance):
                        probabilities = self.normalise(possible_reaction_abundance)
                        transition_reaction = possible_reaction_name[self.get_transition(probabilities)]
                        if any(self.r2rpair[transition_reaction]):
                            for pair in self.r2rpair[transition_reaction]:
                                if compound in pair.split('_'):
                                    transition_compound = pair.replace(compound, '')\
                                                              .replace('_', '')   
                            if transition_compound not in network_visitations[key]:
                                network_visitations[key][transition_compound]=1
                            else:
                                network_visitations[key][transition_compound]+=1    
                            iter_batch[key] = transition_compound
                            previous_reaction=transition_reaction
                            
        output_lines = ['\t'.join(['C'] + self.metadata_keys)]
        for c in set(itertools.chain(*possible_reactions.values())):
            output_line = [c]
            for key in self.metadata_keys:
                if c in network_visitations[key]:
                    output_line.append(str(network_visitations[key][c]))
                else:
                    output_line.append('0')    
            output_lines.append('\t'.join(output_line))
        return output_lines

    def all_matrix(self, 
                   abundances_metagenome, 
                   abundances_transcriptome,
                   abundances_expression,
                   abundances_metabolome,
                   reference_dict):
        '''
        Parameters
        ----------
        '''
        seen_nodes = set()
        # Construct headers to network matrices
        if(abundances_transcriptome and abundances_expression):
            network_lines  = ['\t'.join(self.matrix_header + 
                                        self.transcriptome_header)]
        else:
            network_lines  = ['\t'.join(self.matrix_header)]  
        if abundances_metabolome:
            node_metadata_lines = ['\t'.join(self.metadata_header + 
                                             self.compound_header)]
        else:     
            node_metadata_lines = ['\t'.join(self.metadata_header)]
        
        for reaction, entry in reference_dict.items():
            for compound in entry:
                if any([reaction in x.keys() 
                        for x in abundances_metagenome.values()]):
                    #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                    #~#~#~#~#~#~#~#~#~#~ Fill in abundances ~#~#~#~#~#~#~#~#~#~
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
                    if output_line not in network_lines:
                        network_lines.append(output_line)
                        
                    #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                    #~#~#~#~#~#~#~#~#~ Fill in node metadata ~#~#~#~#~#~#~#~#~#                    
                    if compound not in seen_nodes:
                        if compound in self.compound_desc_dict:
                            compound_type = \
                                ','.join(self.compound_desc_dict[compound]['A'])
                        else:
                            compound_type = 'NA'
                        
                        metadata_list = [compound, self.c[compound], compound_type, 
                                         'NA', 'NA', 'NA', 'NA', 'compound']
                        if abundances_metabolome:   
                            # This is fragile - compound abundances need to be aggregated
                            # by the user. TODO: Integrate compounds into kegg_matrix
                            # module.
                            for sample in self.metadata_keys:
                                c_ab = abundances_metabolome.get_entry(sample, compound)
                                if c_ab == 0:
                                    c_ab = '-5'
                                metadata_list.append(str(c_ab))
                        node_metadata_lines.append('\t'.join(metadata_list))
                        seen_nodes.add(compound)

                    if reaction not in seen_nodes:
                        module, module_description \
                                        = self._gather_module(reaction)
                        pathway, pathway_description \
                                        = self._gather_pathway(reaction)
                        metadata_list = [reaction, 
                                         self.r[reaction],
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
        query_list    = self._parse_queries(queries)
        check_list    = set(query_list.keys())
        seen_steps    = set()
        seen_nodes    = set()
        level_queries = set()

        if(abundances_transcriptome and abundances_expression):
            network_lines  \
                    = ['\t'.join(self.matrix_header + 
                                 self.transcriptome_header)]
        else:
            network_lines  \
                    = ['\t'.join(self.matrix_header)]
        node_metadata_lines \
                    = ['\t'.join(self.metadata_header + 
                                 self.compound_header +
                                 self.query_header)]
        compound_reaction_index_lines \
                    = ['\t'.join(self.compound_reaction_index_header)]
        while depth>0:
            if any(level_queries):
                check_list = set(level_queries)
                level_queries = set()
            for reaction, entry in self.r2c.items():
                if any([reaction in x.keys() 
                        for x in abundances_metagenome.values()]):
                    if any(check_list.intersection(entry)):
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
                                         in check_list else steps+1))
                            
                            if compound not in seen_nodes:
                                if compound in query_list:
                                    query_ab = [str(query_list[compound][key]) 
                                                for key in self.metadata_keys]
                                else:
                                    query_ab = ['NA' for key 
                                                in self.metadata_keys]
                                
                                if compound in self.compound_desc_dict:
                                    compound_type = \
                                        ','.join(self.compound_desc_dict[compound]['A'])
                                else:
                                    compound_type = 'NA'
                                    
                                if compound in query_list:
                                    is_query = 'True'
                                else:
                                    is_query = 'False'
                                    
                                node_metadata_lines.append(
                                    '\t'.join([compound, 
                                               self.c[compound],
                                               compound_type, 
                                               'NA',
                                               'NA', 
                                               'NA',
                                               'NA',  
                                               'compound'] +
                                               query_ab +
                                               [is_query, index])
                                                           )
                                
                                seen_nodes.add(compound)
                                if compound not in query_list:
                                    level_queries.add(compound)
                            if reaction not in seen_nodes:
                                module, module_description \
                                                = self._gather_module(reaction)
                                pathway, pathway_description \
                                                = self._gather_pathway(reaction)
                                node_metadata_lines.append(
                                    '\t'.join([reaction, 
                                               self.r[reaction],
                                               'NA', 
                                               module, 
                                               module_description, 
                                               pathway,
                                               pathway_description,
                                               'reaction',
                                               'NA',
                                               'NA',
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
                     abundances_transcriptome,
                     abundances_expression,
                     abundances_metabolome,
                     limit,
                     filter,
                     node_from,
                     node_to,
                     catabolic,
                     anabolic,
                     bfs_shortest_path):

        
        possible_reactions=set()

        if any(limit):
            for entry in limit:
                if(entry.startswith(self.MAP_PREFIX) or 
                   entry.startswith(self.PATHWAY_PREFIX)): 
                    for reaction in self.p2r[entry]:
                        possible_reactions.add(reaction)
                elif(entry.startswith(self.MODULE_PREFIX)):
                    for reaction in self.m2r[entry]:
                        possible_reactions.add(reaction)
                elif(entry.startswith(self.REACTION_PREFIX)):
                    possible_reactions.add(entry)

            possible_reactions = {reaction:self.r2c[reaction] 
                                  for reaction in
                                  possible_reactions}                
        else:
            possible_reactions = self.r2c
        
        for entry in filter:
            if entry in possible_reactions:
                del possible_reactions[entry]
 

        if(node_to and node_from):
            if bfs_shortest_path:

                bfs_paths = NetworkTraverser\
                                    .shortest_bfs_path(self.r2rpair, 
                                                       self.c2r,
                                                       node_from,
                                                       node_to)
                shortest_path_reactions = set()
                
                for entry in bfs_paths:
                    if entry.startswith(self.REACTION_PREFIX):
                        shortest_path_reactions.add(entry)
            else:
                pass
            
            possible_reactions = {reaction:possible_reactions[reaction]
                                  for reaction in shortest_path_reactions}
        
        network_lines, node_metadata_lines = \
            self.all_matrix(abundances_metagenome, 
                            abundances_transcriptome, 
                            abundances_expression, 
                            abundances_metabolome,
                            possible_reactions)
        
        return network_lines, node_metadata_lines
  
    def traverse(self,
                 abundances_metagenome,
                 abundances_transcriptome,
                 limit,
                 filter,
                 starting_compounds,
                 steps,
                 number_of_queries):

        if any(limit):
            for entry in limit:
                if(entry.startswith(self.MAP_PREFIX) or 
                   entry.startswith(self.PATHWAY_PREFIX)): 
                    for reaction in self.p2r[entry]:
                        possible_reactions.add(reaction)
                elif(entry.startswith(self.MODULE_PREFIX)):
                    for reaction in self.m2r[entry]:
                        possible_reactions.add(reaction)
                elif(entry.startswith(self.REACTION_PREFIX)):
                    possible_reactions.add(entry)

            possible_reactions = {reaction:self.r2c[reaction] 
                                  for reaction in
                                  possible_reactions}                
        else:
            possible_reactions = self.r2c
        
        for entry in filter:
            if entry in possible_reactions:
                del possible_reactions[entry]
        
        possible_compounds = set(itertools.chain(*possible_reactions.values()))

        if len(starting_compounds)==0:
            query_list=[x for x in self.c.keys()
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