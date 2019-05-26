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
import os
import sys
import tempfile
import filecmp
import shutil

PATH_TO_DATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from enrichm.predict import Predict

###############################################################################

class Tests(unittest.TestCase):

    ml_data = 'ml_data'
    sample_matrix = 'matrix.tsv'
    sample_metadata = 'metadata.tsv'
    sample_generate = 'generate_example'
    sample_predict = 'predict_example'

    sample_matrix_path = os.path.join(PATH_TO_DATA, ml_data, sample_matrix)
    sample_generate_path = os.path.join(PATH_TO_DATA, ml_data, sample_generate)
    sample_predict_path = os.path.join(PATH_TO_DATA, ml_data, sample_predict)

    def test_hello_predict(self):
        tmp = tempfile.mkdtemp()
        predict = Predict()
        predict.predict_pipeline(self.sample_generate_path, self.sample_matrix_path, tmp)

        expected_files = sorted(os.listdir(tmp))
        observed_files = sorted(os.listdir(self.sample_predict_path))

        self.assertEqual(len(expected_files), len(observed_files))

        for expected_file, observed_file in zip(expected_files, observed_files):
            expected_file_path = os.path.join(tmp, expected_file)
            observed_file_path = os.path.join(self.sample_predict_path, observed_file)
            # Are all files present?
            self.assertEqual(expected_file, observed_file)
        
        shutil.rmtree(tmp)
        # Note I chose not to match files exactly here. They change because ml models are 
        # estimations and will be different every time you make them. More tests that ensure the
        # internal functions of generate will be implemented as an alternative.

if __name__ == "__main__":
    unittest.main()
