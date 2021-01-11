#!/usr/bin/env python3
import logging
import os
import numpy as np
from enrichm.parser import ParseGenerate
from enrichm.writer import Writer
from enrichm.parser import Parser

class Predict:

    def __init__(self):
        self.predictions_output_file = 'predictions.tsv'
        self.predictions_header = ["Sample", "Prediction", "Probability"]

    def make_predictions(self, model, sample_list, content_list, attribute_dictionary):
        '''
        Inputs
        ------

        Outputs
        -------

        '''

        sample_list = np.array(sample_list)
        content_list = np.array(content_list)
        predictions = model.predict(content_list)
        probabilities = model.predict_proba(content_list)

        output_list = [self.predictions_header]

        for sample, prediction, probability in zip(sample_list, predictions, probabilities):
            max_prob = str(round(max(list(probability)), 2))
            prediction = str(attribute_dictionary[prediction])
            output_list.append([sample, prediction, max_prob])

        return output_list

    def predict_pipeline(self, forester_model_directory, input_matrix_path, output_directory):
        '''
        Inputs
        ------

        Outputs
        -------

        '''
        forester_model = ParseGenerate(forester_model_directory)

        logging.info('Parsing input')
        logging.info('Loading model: %s' % (forester_model.rf_model))

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

        Writer.write(output_lines, os.path.join(output_directory, self.predictions_output_file))
