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

###############################################################################

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','enrichm')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_annotate = os.path.join(path_to_data, 'enrichm_annotate')
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from enrichm.classifier import Classify
from enrichm.module_description_parser import ModuleDescription

###############################################################################

class Tests(unittest.TestCase):

    def test_classify_from_matrix(self):
        tmp = tempfile.mkdtemp()
        ko_matrix = os.path.join(path_to_annotate, 'ko_frequency_table.tsv')
        cmd = f'{path_to_script} classify \
                    --genome_and_annotation_matrix {ko_matrix} \
                    --output {tmp} \
                    --force \
                    --verbosity 1'
        subprocess.call(cmd, shell=True)

    def test_classify_with_custom_modules(self):
        tmp = tempfile.mkdtemp()
        ko_matrix = os.path.join(path_to_annotate, 'ko_frequency_table.tsv')
        cmd = f'{path_to_script} classify \
                    --genome_and_annotation_matrix {ko_matrix} \
                    --output {tmp} \
                    --force \
                    --verbosity 1'
        subprocess.call(cmd, shell=True)

    def test_curate(self):
        # M00516 
        module = ModuleDescription("K11231,K19690,K19691,K19692 K11232 K11233,K15859")
        total_steps = module.num_steps()
        fails1 = ['K11231', 'K11232']
        passes1 = ['K11231', 'K11232', 'K11233']
        self.assertEqual(False,
                         module.num_covered_steps(fails1)[0] == total_steps)
        self.assertEqual(True,
                         module.num_covered_steps(passes1)[0] == total_steps)
        # M00156
        module = ModuleDescription("((K00404+K00405,K15862)+K00407+K00406)")
        total_steps = module.num_steps()
        fails1 = ['K15862', 'K00406']
        fails2 = ['K00404', 'K00407', 'K00406']
        passes1 = ['K00404', 'K00405', 'K00407', 'K00406']
        passes2 = ['K15862', 'K00407', 'K00406']
        self.assertEqual(False,
                         module.num_covered_steps(fails1)[0] == total_steps)
        self.assertEqual(True,
                         module.num_covered_steps(passes1)[0] == total_steps)
        self.assertEqual(False,
                         module.num_covered_steps(fails2)[0] == total_steps)
        self.assertEqual(True,
                         module.num_covered_steps(passes2)[0] == total_steps)
        
        # M00006
        module = ModuleDescription("(K13937,((K00036,K19243) (K01057,K07404))) K00033")
        total_steps = module.num_steps()
        fails1 = ['K19243', 'K00033', 'K19243']
        fails2 = ['K19243', 'K00033']
        passes1 = ['K00036', 'K01057', 'K00033']
        passes2 = ['K13937', 'K00033']
        passes3 = ['K19243', 'K07404', 'K00033']
        passes4 = ['K00036', 'K07404', 'K00033']
        self.assertEqual(True,
                         module.num_covered_steps(passes1)[0] == total_steps)
        self.assertEqual(True,
                         module.num_covered_steps(passes2)[0] == total_steps)
        self.assertEqual(True,
                         module.num_covered_steps(passes3)[0] == total_steps)
        self.assertEqual(True,
                         module.num_covered_steps(passes4)[0] == total_steps)
        self.assertEqual(False,
                         module.num_covered_steps(fails1)[0] == total_steps)
        self.assertEqual(False,
                         module.num_covered_steps(fails2)[0] == total_steps)
    
    def test_update(self):
        
        with open(tempfile.mktemp(), 'w') as out_io:
            out_io.write(f"test\tcontents\n")
            out_io.flush()
            classify = Classify()
            classify.update_with_custom_modules(out_io.name)
            self.assertEqual(True, 'test' in classify.modules)
            self.assertEqual('contents', classify.m2def['test'])


if __name__ == "__main__":
    unittest.main()
