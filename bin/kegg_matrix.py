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
from math import sqrt

###############################################################################

class KeggMatrix:
    
    REACTIONS_PICKLE = 'reactions_22-11-2016.pickle'
    MODULES_PICKLE   = 'modules_22-11-2016.pickle'
    MODULE           = 'module'
    REACTION         = 'Reaction'
    ORTHOLOGY        = 'orthology_def'
    
    def __init__(self, matrix):
        self.reactions_dict \
                = pickle.load(open(self.REACTIONS_PICKLE))
        self.modules_dict \
                = pickle.load(open(self.MODULES_PICKLE))
        self.orthology_matrix \
                = self._parse_matrix(matrix)
        self.reaction_matrix \
                = self._calculate_abundances(self.reactions_dict,
                                             self.orthology_matrix)
        self.module_matrix \
                = self._calculate_abundances(self.modules_dict,
                                             self.reaction_matrix)

    def _parse_matrix(self, matrix):
        
        output_dict = {}
        
        for idx, line in enumerate(open(matrix)):
            sline = line.strip().split('\t')
            
            if idx==0:
                self.sample_names = sline[1:]
                    
                for sample in self.sample_names: 
                    output_dict[sample] = {}
            else:   
                ko_id      = sline[0]
                abundances = sline[1:]
                
                for abundance, sample in zip(abundances, self.sample_names):
                    output_dict[sample][ko_id] = float(abundance)
        return output_dict
    
    def group_abundances(self, level, samples):
        
        output_dict = {}
        
        if level == self.MODULE:
            reference_dict = self.module_matrix
        if level == self.REACTION:
            reference_dict = self.reaction_matrix
        if level == self.ORTHOLOGY:
            reference_dict = self.orthology_matrix
        
        for sample in samples:
            new_dict = {key:entry for key,entry in reference_dict.items() if 
                        key in samples}
            reference_list = new_dict[new_dict.keys()[0]].keys()
            
            for reference in reference_list:
                # If samples are missing from stored data, this will crash 
                abundances = [new_dict[sample][reference] for sample in samples]
                average    = sum(abundances)/float(len(abundances))
                output_dict[reference] = average
        return output_dict
                
                
    def _calculate_abundances(self, reference_dict, matrix_dict):
        
        output_dict_mean   = {}
        definition_key = (self.ORTHOLOGY if self.ORTHOLOGY 
                          in reference_dict[reference_dict.keys()[0]]
                          else self.REACTION)

        for sample, ko_abundances in matrix_dict.items():
            output_dict_mean[sample]   = {}
            
            for key, definition in reference_dict.items():
                abundances = []
                for id in definition[definition_key]:
                    if id in matrix_dict[sample]:
                        abundances.append(matrix_dict[sample][id])
                    else:
                        logging.debug("ID not found in input matrix: %s" % id)
                        
                if any(abundances):
                    abundance_mean = sum(abundances)/len(abundances)
                else:
                    abundance_mean = 0

                output_dict_mean[sample][key] = abundance_mean
        
        return output_dict_mean

    
    def calculate_expression(self, km_t):
        '''
        '''
        import IPython ; IPython.embed()