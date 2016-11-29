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
import urllib2     
from itertools import product, chain
from bs4 import BeautifulSoup

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

    REACTIONS_PICKLE     = os.path.join(DATA_PATH, 'reactions_22-11-2016.pickle')
    MODULES_PICKLE       = os.path.join(DATA_PATH, 'modules_22-11-2016.pickle')
    COMPOUNDS_PICKLE     = os.path.join(DATA_PATH, 'compounds_22-11-2016.pickle')
    COMPOUND_DESC_PICKLE = os.path.join(DATA_PATH, 'compound_descriptions.pickle')
    C2M_PICKLE           = os.path.join(DATA_PATH, 'c2m.pickle')
    M2C_PICKLE           = os.path.join(DATA_PATH, 'm2c.pickle')
    
    MODULE           = 'module'
    REACTION         = 'Reaction'
    ORTHOLOGY        = 'orthology_def'
    REACTANT_PAIRS   = 'reactant_pairs'
    MAPS             = 'maps'
    COMPOUNDS        = 'compound_def'
    MODULE           = 'modules'
    KEGG_C_API       = 'http://rest.kegg.jp/link/compound/'
    KEGG_L_API       = 'http://rest.kegg.jp/list/'
    
    R2C = 'http://rest.kegg.jp/link/compound/reaction'
    R2M = 'http://rest.kegg.jp/link/module/reaction'
    M2R = 'http://rest.kegg.jp/link/reaction/module'
    R2P = 'http://rest.kegg.jp/link/pathway/reaction'
    P2R = 'http://rest.kegg.jp/link/reaction/pathway'

    C   = 'http://rest.kegg.jp/list/compound'    
    R   = 'http://rest.kegg.jp/list/reaction'
    P   = 'http://rest.kegg.jp/list/pathway'
    M   = 'http://rest.kegg.jp/list/module'

    def __init__(self):
        logging.info("Downloading reaction to pathway information from KEGG")
        self.r2p = self.build_dict(self.R2P)
        logging.info("Done")
        logging.info("Downloading pathway to reaction information from KEGG")
        self.p2r = self.build_dict(self.P2R)
        logging.info("Done")
        logging.info("Downloading reaction to module information from KEGG")
        self.r2m = self.build_dict(self.R2M)
        logging.info("Done")
        logging.info("Downloading reaction to module information from KEGG")
        self.m2r = self.build_dict(self.M2R)
        logging.info("Done")
        logging.info("Downloading reaction to compound information from KEGG")
        self.r2c = self.build_dict(self.R2C)
        logging.info("Done")
        logging.info("Downloading compound descriptions from KEGG")
        self.c   = self.build_name_dict(self.C)
        logging.info("Done")
        logging.info("Downloading pathway descriptions from KEGG")
        self.p   = self.build_name_dict(self.P)
        logging.info("Done")
        logging.info("Downloading reaction descriptions from KEGG")
        self.r   = self.build_name_dict(self.R)
        logging.info("Done")
        logging.info("Downloading module descriptions from KEGG")
        self.m   = self.build_name_dict(self.M)
        logging.info("Done")
        
        self.reactions_dict \
                = pickle.load(open(self.REACTIONS_PICKLE))
        self.modules_dict \
                = pickle.load(open(self.MODULES_PICKLE))
        self.compounds_dict \
                = pickle.load(open(self.COMPOUNDS_PICKLE))
        self.compound_desc_dict \
                = pickle.load(open(self.COMPOUND_DESC_PICKLE))
        self.c2m \
                = pickle.load(open(self.C2M_PICKLE))
        self.m2c \
                = pickle.load(open(self.M2C_PICKLE))
        self.maps = {x:[] 
                     for x in set(chain(*[self.modules_dict[x]['maps'] 
                                          for x in self.modules_dict.keys()]))}
        for map in self.maps.keys():
            for module, entry in self.modules_dict.items():
                if map in entry[self.MAPS]:
                    self.maps[map].append(module)
        self.compound_names={}
    
    @staticmethod
    def build_name_dict(url):
        output_dictionary = {}
        for entry in urllib2.urlopen(url).read().strip().split('\n'):
            sentry = entry.split('\t')
            key = sentry[0].split(':')[1]
            item = sentry[1].split(';')[0]
            output_dictionary[key]=item
        return output_dictionary    
    
    @staticmethod
    def build_dict(url):
        output_dictionary = {}
        for entry in urllib2.urlopen(url).read().strip().split('\n'):
            key, item = [x.split(':')[1] for x in entry.split('\t')]
            if key in output_dictionary:
                output_dictionary[key].append(item)
            else:
                output_dictionary[key] = [item]
        return output_dictionary
            
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
    
    def _gather_compound_name(self, compound_list):   

        soup = BeautifulSoup(urllib.urlopen("%s%s" % (self.KEGG_L_API,'+'.join(compound_list))), "html.parser")
        self.compound_names.update({x.split('\t')[0].split(':')[1]:x.split('\t')[1].split(';')[0] 
                         for x in soup.getText().strip().split('\n')})
    
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

    def query_matrix(self, abundances_1, abundances_2, queries, depth):
        
        query_list = self._parse_queries(queries)
        check_list = set(query_list.keys())
        
        level_queries = set()
        output_lines=['\t'.join(["compound", "reaction", 'FC_col', 'FC', 
                                'ab_1', 'ab_2', 'C', 'module', 
                                'module_description', 'compound_description', 
                                'reaction_description', 'compound_type', 
                                'query', 'pathway', 'pathway_descriptions',
                                'c1_ab', 'c2_ab'])]
        to_omit = set([x for x,y in self.compound_desc_dict.items() 
                       if "Vitamins and Cofactors" in y['A']])
        to_omit.add('C00001')
        to_omit.add('C00008')
        to_omit.add('C00013')
        to_omit.add('C00004')
        to_omit.add('C00005')
        to_omit.add('C00080')
        to_omit.add('C00009')
        to_omit.add('C00008') # O2
        to_omit.add('C00004')
        to_omit.add('C00020')
        to_omit.add('C00007')
        while depth>0:
            if any(level_queries):
                check_list = level_queries
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
                            elif abundance_2>0 :
                                C=self.COL_UNIQUE_2   
                            else:
                                continue             
                    else:
                        continue
                    reaction_compounds = [x for x in entry if x not in to_omit]

                    for compound in reaction_compounds:
                        level_queries.add(compound)
                        compound_description = self.c[compound]
                        query = ('True' if compound in query_list
                                 else 'False')
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
                                         str(abundance_1), str(abundance_2), C, module, 
                                         module_description, compound_description,
                                         reaction_description, compound_type,
                                         query,pathway,pathway_description, 
                                         c1_ab, c2_ab] 
                        output_line = '\t'.join(reaction_line)
                        if output_line not in output_lines:
                            output_lines.append(output_line)
            depth-=1
        return output_lines
