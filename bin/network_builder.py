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
from cookielib import reach
 
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
from itertools import product, chain

###############################################################################

class NetworkBuilder:
    COL_BLACK        = 'Black'
    COL_RED          = 'Red'
    COL_BLUE         = 'Blue'
    COL_UNIQUE_1     = 'Unique_1'
    COL_UNIQUE_2     = 'Unique_2'

    REACTIONS_PICKLE = 'reactions_22-11-2016.pickle'
    MODULES_PICKLE   = 'modules_22-11-2016.pickle'
    COMPOUNDS_PICKLE = 'compounds_22-11-2016.pickle'
    MODULE           = 'module'
    REACTION         = 'Reaction'
    ORTHOLOGY        = 'orthology_def'
    REACTANT_PAIRS   = 'reactant_pairs'
    MAPS             = 'maps'
    COMPOUNDS        = 'compound_def'
    MODULE           = 'modules'

    def __init__(self):
        self.reactions_dict \
                = pickle.load(open(self.REACTIONS_PICKLE))
        self.modules_dict \
                = pickle.load(open(self.MODULES_PICKLE))
        self.compounds_dict \
                = pickle.load(open(self.COMPOUNDS_PICKLE))
        self.maps = {x:[] 
                     for x in set(chain(*[self.modules_dict[x]['maps'] 
                                          for x in self.modules_dict.keys()]))}
        for map in self.maps.keys():
            for module, entry in self.modules_dict.items():
                if map in entry[self.MAPS]:
                    self.maps[map].append(module)
        
    def build(self, abundances_1, abundances_2):
                
        for map, module_list in self.maps.iteritems():
            output_lines = ['\t'.join(["Compound_1", "Compound_2", "Reaction", 
                                   "FC", "Group_1_abundance", 
                                   "Group_2_abundance", "C", "Module",
                                   "Map", "Compound_1_modules",
                                   "Compound_2_modules"])]
            
            for module in module_list:
                for reaction in self.modules_dict[module][self.REACTION]:
                    
                    C=self.COL_BLACK
                    FC=0
                    
                    if reaction in abundances_1 and \
                       reaction in abundances_2:
                        
                        abundance_1 = abundances_1[reaction]
                        abundance_2 = abundances_2[reaction]
                        
                        if abundance_1>0 and \
                           abundance_2>0:
                            FC = max(abundance_1, abundance_2) /\
                                 min(abundance_1, abundance_2)
                            if abundance_2>abundance_1:
                                C=self.COL_RED
                            else:
                                FC=FC*-1
                                C=self.COL_BLUE
                        else:
                            FC=0
                            if abundance_1>0:
                                C=self.COL_UNIQUE_1
                            elif abundance_2>0 :
                                C=self.COL_UNIQUE_2   
                            else:
                                continue             
                    else:
                        continue
                        #raise Exception("Reaction %s does not exist in input matrix" \
                        #                                % reaction)
                    # Pesky unicode
                    rps = \
                        chain(*[rp.encode('ascii', 'replace').split('??') 
                        for rp in self.reactions_dict[reaction][self.REACTANT_PAIRS]])
                    
                    for rp in rps:
                        compound_1, compound_2 = rp.split('_') 
                        
                        if compound_1 in self.compounds_dict:
                            compound_1_modules = [x for x in self.compounds_dict[compound_1] if x in module_list]
                        
                            if any(compound_1_modules):
                                compound_1_modules = ','.join(compound_1_modules)
                            else:
                                compound_1_modules = 'NA'
                        else:
                            compound_2_modules = 'NA'
                            
                        if compound_2 in self.compounds_dict:
                            compound_2_modules = [x for x in self.compounds_dict[compound_2] if x in module_list]
                            if any(compound_2_modules):
                                compound_2_modules = ','.join(compound_2_modules)
                            else:
                                compound_2_modules = 'NA'
                        else:
                            compound_2_modules = 'NA'
                            
                        modules = self.reactions_dict[reaction][self.MODULE]
                        maps = []
                        
                        for module in modules:
                            maps += self.modules_dict[module][self.MAPS]
                        
                        modules = ','.join(modules)
                        maps    = ','.join(maps)
                        reaction_line = [compound_1, compound_2, reaction, 
                                         str(FC), str(abundance_1), 
                                         str(abundance_2), C, module, maps,
                                         compound_1_modules, compound_2_modules] 
                        output_lines.append('\t'.join(reaction_line))
            yield map, output_lines
                    
