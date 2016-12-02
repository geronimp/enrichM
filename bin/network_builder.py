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
__copyright__ = "Copyright 2015"
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

###############################################################################

class NetworkBuilder:
    DATA_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             '..', 
                             'data')

    COL_BLACK        = 'Black'
    COL_RED          = 'Red'
    COL_BLUE         = 'Blue'
    COL_UNIQUE_1     = 'Unique_1'
    COL_UNIQUE_2     = 'Unique_2'
    
    VERSION = os.path.join(DATA_PATH, 'VERSION')
    COMPOUND_DESC_PICKLE = os.path.join(DATA_PATH, 'br08001')
    
    PICKLE = 'pickle'
    
    R2C = os.path.join(DATA_PATH, 'reaction_to_compound')
    R2M = os.path.join(DATA_PATH, 'reaction_to_module')
    M2R = os.path.join(DATA_PATH, 'module_to_reaction')
    R2P = os.path.join(DATA_PATH, 'reaction_to_pathway')
    P2R = os.path.join(DATA_PATH, 'pathway_to_reaction')

    C   = os.path.join(DATA_PATH, 'compound_descriptions')    
    R   = os.path.join(DATA_PATH, 'reaction_descriptions')
    P   = os.path.join(DATA_PATH, 'pathway_descriptions')
    M   = os.path.join(DATA_PATH, 'module_descriptions')

    def __init__(self):
        
        self.VERSION = open(self.VERSION).readline().strip()
        
        logging.info("Loading reaction to pathway information")
        self.r2p = pickle.load(open('.'.join([self.R2P, 
                                              self.VERSION, self.PICKLE])))
        logging.info("Done")
        logging.info("Loading pathway to reaction information")
        self.p2r = pickle.load(open('.'.join([self.P2R, 
                                                 self.VERSION, self.PICKLE])))
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
        for line in open(queries):
            sline = line.strip().split()
            output_dict[sline[0]] = [sline[1], sline[2]]
        return output_dict
    
    def all_matrix(self, abundances_1, abundances_2):
        output_lines=['\t'.join(["compound", "reaction", 'FC_col', 'FC', 
                                'ab_1', 'ab_2', 'C', 'module', 
                                'module_description', 'compound_description', 
                                'reaction_description', 'compound_type'])]
        
        for reaction, entry in self.r2c.items():
            C=self.COL_BLACK
            FC=0
            if(reaction in abundances_1 and reaction in abundances_2):
                abundance_1 = abundances_1[reaction]
                abundance_2 = abundances_2[reaction]
                 
                if abundance_1>0 and \
                   abundance_2>0:
                    FC = max(abundance_1, abundance_2) /\
                         min(abundance_1, abundance_2)
                    if abundance_2>abundance_1:
                        C=self.COL_RED
                        FC_col=FC
                    else:
                        FC_col=FC*-1
                        C=self.COL_BLUE
                else:
                    FC_col=0
                    FC=0
                    if abundance_1>0:
                        C=self.COL_UNIQUE_1
                    elif abundance_2>0 :
                        C=self.COL_UNIQUE_2   
                    else:
                        continue             
            else:
                continue
            for compound in entry:
                compound_description = self.c[compound]
                reaction_description = self.r[reaction]
                module, module_description = self._gather_module(reaction)

                if compound in self.compound_desc_dict:
                    compound_type = \
                        ','.join(self.compound_desc_dict[compound]['A'])
                else:
                    compound_type = 'NA'
                reaction_line = [compound, reaction, str(FC_col), str(FC), 
                                 str(abundance_1), str(abundance_2), C, module, 
                                 module_description, compound_description,
                                 reaction_description, compound_type] 
                output_line = '\t'.join(reaction_line)
                if output_line not in output_lines:
                    output_lines.append(output_line)
        return output_lines

    def query_matrix(self, abundances_1, abundances_2, queries, depth, 
                     abundances_1_expression, abundances_2_expression,
                     group1_transcriptome_abundances,
                     group2_transcriptome_abundances, 
                     abundances_1_name, abundances_2_name):
               
        steps=0
        query_list = self._parse_queries(queries)
        check_list = set(query_list.keys())
        seen_steps = set()
        seen_nodes = set()
        level_queries = set()
        output_lines=['\t'.join(["compound", "reaction", 'FC_col', 'FC', 
                                '%s_reaction_ab' % abundances_1_name,
                                '%s_reaction_ab' % abundances_2_name, 'C', 
                                '%s_compound_ab' % abundances_1_name, 
                                '%s_compound_ab' % abundances_2_name])]
        node_metadata = ['\t'.join(['node', 'description','type',
                                    'query', 'module', 'module_descr',
                                    'pathway', 'pathway_descr', 'step',
                                    'node_type'])]
        if(abundances_1_expression and 
           abundances_2_expression and
           group1_transcriptome_abundances and 
           group1_transcriptome_abundances):
            output_lines[0]+= '\t' + \
                    '\t'.join(["%s_reaction_transcriptome" % abundances_1_name,
                               "%s_reaction_transcriptome" % abundances_2_name,
                               "%s_reaction_expression" % abundances_1_name,
                               "%s_reaction_expression" % abundances_2_name])                         
            
        to_omit = set([x for x,y in self.compound_desc_dict.items() 
                       if "Vitamins and Cofactors" in y['A']])
        to_omit.add('C00001') # H2O
        to_omit.add('C00008') # ADP
        to_omit.add('C00013') # Diphosphate
        to_omit.add('C00004') # NADH
        to_omit.add('C00005') # NADPH
        to_omit.add('C00080') # H+
        to_omit.add('C00009') # Orthophosphate
        to_omit.add('C00008') # ADP
        to_omit.add('C00004') # NADH
        to_omit.add('C00020') # AMP
        to_omit.add('C00007') # Oxygen
        to_omit.add('C00015') # UDP
        
        while depth>0:
            if any(level_queries):
                check_list = set(level_queries)
                level_queries = set()
            for reaction, entry in self.r2c.items():
                if any(check_list.intersection(entry)):
                    C=self.COL_BLACK
                    FC=0
                    if(reaction in abundances_1 and reaction in abundances_2):
                        abundance_1 = abundances_1[reaction]
                        abundance_2 = abundances_2[reaction]
                         
                        if abundance_1>0 and \
                           abundance_2>0:
                            FC = max(abundance_1, abundance_2) /\
                                 min(abundance_1, abundance_2)
                            if abundance_2>abundance_1:
                                C=self.COL_RED
                                FC_col=FC
                            else:
                                FC_col=FC*-1
                                C=self.COL_BLUE
                        else:
                            FC_col=0
                            FC=0
                            if abundance_1>0:
                                C=self.COL_UNIQUE_1
                            elif abundance_2>0:
                                C=self.COL_UNIQUE_2
                            else:
                                continue             
                    else:
                        continue
                    reaction_compounds = [x for x in entry if x not in to_omit]
                    for compound in reaction_compounds:
                        level_queries.add(compound)
                        
                        query = ('True' if compound in query_list
                                 else 'False')
                        
                        compound_description = self.c[compound]
                        reaction_description = self.r[reaction]
                        
                        module, module_description = self._gather_module(reaction)
                        pathway, pathway_description = self._gather_pathway(reaction)
                            
                        if compound in query_list:
                            c1_ab, c2_ab = query_list[compound]
                        else:
                            c1_ab, c2_ab = ('-10', '-10')
                        if compound in self.compound_desc_dict:
                            compound_type = \
                                ','.join(self.compound_desc_dict[compound]['A'])
                        else:
                            compound_type = 'NA'
                        
                        reaction_line = [compound, reaction, str(FC_col), str(FC), 
                                         str(abundance_1), str(abundance_2), C, 
                                         c1_ab, c2_ab] 
                        
                        index = str((steps if compound 
                                     in check_list else steps+1))

                        if(abundances_1_expression and abundances_2_expression):
                            if reaction in(abundances_1_expression and abundances_2_expression):                                
                                reaction_line.append(str(group1_transcriptome_abundances[reaction]))
                                reaction_line.append(str(group2_transcriptome_abundances[reaction]))
                                reaction_line.append(str(abundances_1_expression[reaction]))
                                reaction_line.append(str(abundances_2_expression[reaction]))

                                
                        output_line = '\t'.join(reaction_line)
                        if output_line not in seen_steps:
                            seen_steps.add(output_line)
                            output_lines.append(output_line+'\t%i' % steps)
                        
                        if compound not in seen_nodes:
                            node_metadata.append('\t'.join([compound, 
                                                compound_description,
                                                compound_type, query, 'NA',
                                                'NA', 'NA', 'NA',index,  
                                                'compound']))
                        if reaction not in seen_nodes:
                            node_metadata.append('\t'.join([reaction, reaction_description,
                                                'NA', 'False', module, 
                                                module_description, pathway,
                                                pathway_description, index,
                                                'reaction']))
                        
                        
            steps+=1
            depth-=1
            logging.info("Step %i complete with %i queries to continue with" \
                                              % (steps, len(level_queries)))
        return output_lines, node_metadata
