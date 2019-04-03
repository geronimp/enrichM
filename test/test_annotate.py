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
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

###############################################################################

class Tests(unittest.TestCase):

    def test_hello_world_nucleic(self):
        tmp = tempfile.mkdtemp()
        bin = os.path.join(path_to_data, 'test_nucleic_bin')
        cmd = '%s annotate \
                        --threads 10 \
                        --ko \
                        --pfam \
                        --tigrfam \
                        --genome_directory %s \
                        --output %s \
                        --force' % (path_to_script, bin, tmp)
        subprocess.call(cmd, shell=True)

    def test_hello_world_protein(self):
        tmp = tempfile.mkdtemp()
        bin = os.path.join(path_to_data, 'test_protein_bin')
        cmd = '%s annotate \
                        --threads 10 \
                        --ko \
                        --pfam \
                        --tigrfam \
                        --protein_directory %s \
                        --output %s \
                        --force' % (path_to_script, bin, tmp)
        subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    unittest.main()
