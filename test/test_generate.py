#!/usr/bin/env python
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
# Imports
import unittest
import os.path
import sys
import subprocess
import tempfile
import filecmp

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','enrichm')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_annotate = os.path.join(path_to_data, 'enrichm_annotate')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from enrichm.generate import GenerateModel

###############################################################################

class Tests(unittest.TestCase):
    
    ML_DATA = 'ml_data'
    sample_matrix = 'matrix.tsv'
    sample_metadata = 'metadata.tsv'
    
    sample_matrix_path = os.path.join(path_to_data, ML_DATA, 'matrix.tsv')
    sample_metadata_path = os.path.join(path_to_data, ML_DATA, 'metadata.tsv')

    def test_hello_generate(self):
        tmp = tempfile.mkdtemp()
        generateModel = GenerateModel()
        generateModel.do(self.sample_matrix_path,
                         self.sample_metadata_path,
                         generateModel.CLASSIFIER,
                         0.2, # Default testing portion
                         False, # Dont do a grid search for fine tuning
                         2, # Threads
                         tmp # Output directory
                         )
        
if __name__ == "__main__":
    unittest.main()
