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

__author__      = "Joel Boyd"
__copyright__   = "Copyright 2018"
__credits__     = ["Joel Boyd"]
__license__     = "GPL3"
__version__     = "0.0.7"
__maintainer__  = "Joel Boyd"
__email__       = "joel.boyd near uq.net.au"
__status__      = "Development"

###############################################################################
# Imports
import logging
import pickle
import os
import numpy as np
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
# Local
from enrichm.generate import GenerateModel
################################################################################


class Predict():

	def __init__(self):
		self.PREDICTIONS_OUTPUT_PATH = 'predictions.tsv'

		self.LABELS_DICT = "labels_dict.pickle"
		self.RF_MODEL = "rf_model.pickle"
		self.ATTRIBUTE_LIST = "attribute_list.txt"

	def _write_predictions(self, output_lines, output_directory):
		'''		
		Inputs
		------
		
		Outputs
		-------
		
		'''
		output_path = os.path.join(output_directory, self.PREDICTIONS_OUTPUT_PATH)
		logging.info('Writing predictions to output %s' % (output_path))
		with open(output_path, 'wb') as out_io:
			out_io.write(str.encode('\t'.join(["Sample", "Prediction", "Probability"]) + '\n'))
			for line in output_lines:
				out_io.write(str.encode(line + '\n'))

	def _make_predictions(self, model, sample_list, content_list, attribute_dictionary):
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

		output_list = []

		for sample, prediction, probability in zip(sample_list, predictions, probabilities):
			max_prob = str(round(max(list(probability)), 2))
			prediction = str(attribute_dictionary[prediction])
			output_list.append('\t'.join([sample, prediction, max_prob]))
		
		return output_list

	def parse_input_model_directory(self, forester_model_directory):

		output_dictionary = {
			self.LABELS_DICT: None,
			self.RF_MODEL: None,
			self.ATTRIBUTE_LIST:None
		}

		contents = os.listdir(forester_model_directory)

		for content in contents:
			content_path = os.path.join(forester_model_directory, content)
			if content in output_dictionary:
				if content.endswith("pickle"):
					output_dictionary[content] = pickle.load(open(content_path, 'rb'))
				else:
					output_dictionary[content] = [x.strip() for x in open(content_path)]

		if None in list(output_dictionary.values()):
			raise Exception("Malformatted forester model directory: %s" % (forester_model_directory))

		return output_dictionary


	def do(self, forester_model_directory, input_matrix_path, output_directory):
		'''		
		Inputs
		------
		
		Outputs
		-------
		
		'''
		forester_model = self.parse_input_model_directory(forester_model_directory)		

		logging.info('Parsing input')
		gm = GenerateModel()
		logging.info('Loading model: %s' % (self.RF_MODEL))

		logging.info('Parsing data')
		features, _ \
			= gm.parse_input_matrix(input_matrix_path)
		
		sample_list = list()
		content_list = list()
		for sample, content in features.items():
			sample_list.append(sample)
			sample_content = []
			for attribute in forester_model['attribute_list.txt']:
				if attribute in content:
					sample_content.append(content[attribute])
				else:	
					sample_content.append('0')
					
			content_list.append(sample_content)

		logging.info('Making predictions')
		output_lines = self._make_predictions(forester_model[self.RF_MODEL],
											  sample_list,
										  	  content_list,
										  	  forester_model[self.LABELS_DICT])
		self._write_predictions(output_lines, output_directory)
