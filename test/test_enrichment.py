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

###############################################################################

path_to_script 		= os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','enrichm')
path_to_data 		= os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_annotate	= os.path.join(path_to_data, 'enrichm_annotate')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

###############################################################################

class Tests(unittest.TestCase):

    def test_enrichment_from_ko_matrix(self):

        tmp = tempfile.mkdtemp()
        bin             = os.path.join(path_to_data, 'test_nucleic_bin')
        expected_output = os.path.join(path_to_data, 'enrichm_enrichment_ko')
        metadata        = os.path.join(path_to_data, 'metadata.tsv')
        cmd             = '%s enrichment --annotate_output %s --metadata %s --output %s --force --ko' \
                            % (path_to_script, path_to_annotate, metadata, tmp)
        subprocess.call(cmd, shell=True)
        self.assertTrue(filecmp.dircmp(tmp, expected_output))
        
        for file in os.listdir(tmp):
            if file         == 'enrichment.log': continue
            output_file     = os.path.join(tmp, file)
            expected_file   = os.path.join(expected_output, file)
            self.assertTrue(filecmp.cmp(output_file, expected_file))

    def test_enrichment_from_pfam_matrix(self):
        
        tmp = tempfile.mkdtemp()
        bin             = os.path.join(path_to_data, 'test_nucleic_bin')
        expected_output = os.path.join(path_to_data, 'enrichm_enrichment_pfam')
        metadata        = os.path.join(path_to_data, 'metadata.tsv')
        genomes_to_compare = os.path.join(path_to_data, 'genomes_to_compare.tsv')
        cmd             = '%s enrichment --annotate_output %s --metadata %s --output %s --force --pfam' \
                            % (path_to_script, path_to_annotate, metadata, tmp)

        subprocess.call(cmd, shell=True)
        
        self.assertTrue(filecmp.dircmp(tmp, expected_output))
        for file in os.listdir(tmp):
            if file         == 'enrichment.log': continue
            output_file     = os.path.join(tmp, file)
            expected_file   = os.path.join(expected_output, file)
            self.assertTrue(filecmp.cmp(output_file, expected_file))

if __name__ == "__main__":
    unittest.main()
