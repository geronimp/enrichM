#!/usr/bin/env python3
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


import logging
import pickle
import os
import numpy as np
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from enrichm.parser import ParseGenerate
from enrichm.writer import Writer
from enrichm.parser import Parser
from enrichm.generate import GenerateModel

class Predict:

    def __init__(self):
        self.PREDICTIONS_OUTPUT_PATH = 'predictions.tsv'
        self.PREDICTIONS_HEADER = ["Sample", "Prediction", "Probability"]

    def make_predictions(self, model, sample_list, content_list, attribute_dictionary):
        '''
        Inputs
        ------

        Outputs
        -------

        '''

        sample_list 	= np.array(sample_list)
        content_list 	= np.array(content_list)
        predictions 	= model.predict(content_list)
        probabilities 	= model.predict_proba(content_list)

        output_list = [self.PREDICTIONS_HEADER]

        for sample, prediction, probability in zip(sample_list, predictions, probabilities):
            max_prob = str(round(max(list(probability)), 2))
            prediction = str(attribute_dictionary[prediction])
            output_list.append([sample, prediction, max_prob])

        return output_list

    def do(self, forester_model_directory, input_matrix_path, output_directory):
        '''
        Inputs
        ------

        Outputs
        -------

        '''
        forester_model = ParseGenerate(forester_model_directory)

        logging.info('Parsing input')
        logging.info('Loading model: %s' % (forester_model.RF_MODEL))

        logging.info('Parsing data')
        features, _, _ = Parser.parse_simple_matrix(input_matrix_path)

        sample_list = list()
        content_list = list()

        for sample, content in features.items():
            sample_list.append(sample)
            sample_content = list()

            for attribute in forester_model.attributes:

                if attribute in content:
                    sample_content.append(content[attribute])
                else:
                    sample_content.append('0')

            content_list.append(sample_content)

        logging.info('Making predictions')
        output_lines = self.make_predictions(forester_model.model,
                                              sample_list,
                                                content_list,
                                                forester_model.labels)
        Writer.write(output_lines, os.path.join(output_directory, self.PREDICTIONS_OUTPUT_PATH))
