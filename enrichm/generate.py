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
import pickle
import os
import logging
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import train_test_split, RandomizedSearchCV, GridSearchCV
# Local
import enrichm.generate
################################################################################

class GenerateModel():

    def __init__(self):
        '''     
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        # Subparser names
        self.REGRESSOR = "regressor"
        self.CLASSIFIER = "classifier"

        # Output file names
        self.ATTRIBUTE_IMPORTANCES = 'attribute_importances.tsv'
        self.ATTRIBUTE_LIST = "attribute_list.txt"
        self.MODEL_PICKLE = "rf_model.pickle"
        self.LABELS_DICT = "labels_dict.pickle"
    

    
    def parse_input_matrix(self, input_matrix_path):
        '''     
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        input_matrix_io = open(input_matrix_path)
        headers = input_matrix_io.readline().strip().split('\t')[1:]
        
        output_dictionary = {header:{} for header in headers}
        attribute_list = list()

        for line in input_matrix_io:
            values = line.strip().split('\t')
            attribute_value = values[0]
            attribute_list.append(attribute_value)
        
            for header, value in zip(headers, values[1:]):
                output_dictionary[header][attribute_value] = value
        logging.info('\t\tRead in matrix of %i samples associated with %i attributes' \
                % (len(output_dictionary), len(attribute_list)))
        return output_dictionary, attribute_list

    def parse_groups_matrix(self, groups_matrix_path):
        '''     
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        output_dictionary = dict()

        for line in open(groups_matrix_path):
            attribute_value, group = line.strip().split('\t')

            if group in output_dictionary:
                raise Exception("Duplicated entry in groups file: %s" % attribute_value)
            else:
                output_dictionary[attribute_value] = group

        logging.info('\t\tRead in metadata for %i samples' % (len(output_dictionary)))
        
        return output_dictionary

    def _numerify(self, input_list):
        '''
        
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        idx = 0
        output_dictionary = {}
        output_list = []
        
        for group in input_list:
            if group not in output_dictionary:
                output_dictionary[group] = idx
                idx += 1
            output_list.append(output_dictionary[group])
        output_dictionary = {item:key for key, item in output_dictionary.items()}
        return output_dictionary, output_list

    def _write_attribute_list(self, attribute_list, output_directory):
        with open(os.path.join(output_directory, self.ATTRIBUTE_LIST), 'wb') as out_io:
            out_io.write(str.encode('\n'.join(attribute_list)))
    
    def _write_importances(self, model, attribute_list, output_directory):
        '''
        
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        importances = list(model.feature_importances_)
        # List of tuples with variable and importance
        feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(attribute_list, importances)]
        # Sort the feature importances by most important first
        feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
        # Print out the feature and importances 
        logging.info('%i attributes found with an importance > 0' % (len([x for x in feature_importances if x[1]>0])))
        logging.info('Writing attribute importances')
        with open(os.path.join(output_directory, self.ATTRIBUTE_IMPORTANCES), 'wb') as out_io:
            out_io.write(str.encode('\t'.join(['Variable', 'Importance']) + '\n'))
            for pair in feature_importances:
                var, imp = pair
                out_io.write(str.encode('\t'.join([str(var), str(imp)]) + '\n'))
    
    def transpose(self, labels, features, attribute_list):
        '''
        
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        labels_list     = list()
        features_list   = list()

        for col in list(labels.keys()):

            labels_list.append(labels[col])
            col_values = list()

            for attribute in attribute_list:
                col_values.append(features[col][attribute])
            
            features_list.append(col_values)

        return labels_list, features_list

    def grid_search_cv(self, random_search_cv, threads, rf):
        '''
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        logging.info('Generating model for grid search cross validation')

        bootstrap = [random_search_cv.best_params_['bootstrap']]
        max_depth = range(random_search_cv.best_params_['max_depth']-9, random_search_cv.best_params_['max_depth']+10, 2)
        max_features = [random_search_cv.best_params_['max_features'], 2, 3]
        min_samples_leaf = range(random_search_cv.best_params_['min_samples_leaf']-1, random_search_cv.best_params_['min_samples_leaf']+2)
        min_samples_split = range(random_search_cv.best_params_['min_samples_split']-2, random_search_cv.best_params_['min_samples_split']+3, 2)
        n_estimators = range(random_search_cv.best_params_['n_estimators']-200, random_search_cv.best_params_['n_estimators'] + 200, 100)

        param_grid = {
            'bootstrap':bootstrap,
            'max_depth':[x for x in max_depth if x>0],
            'max_features':max_features,
            'min_samples_leaf':[x for x in min_samples_leaf if x>0],
            'min_samples_split':[x for x in min_samples_split if x>0],
            'n_estimators':[x for x in n_estimators if x>0]
        }
        

        grid_search = GridSearchCV(estimator = rf,
                                   param_grid = param_grid, 
                                   cv = 3,
                                   n_jobs = threads)
        return grid_search

    def random_search_cv(self, threads, rf):
        '''
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        logging.info('Generating model for random search cross validation')
        
        n_estimators            = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
        max_features            = ['auto', 'sqrt']
        max_depth               = [int(x) for x in np.linspace(10, 110, num = 11)]
        max_depth.append(None)
        min_samples_split       = [2, 5, 10]
        min_samples_leaf        = [1, 2, 4]
        bootstrap               = [True, False]
        random_grid = {'n_estimators': n_estimators,
                       'max_features': max_features,
                       'max_depth': max_depth,
                       'min_samples_split': min_samples_split,
                       'min_samples_leaf': min_samples_leaf,
                       'bootstrap': bootstrap}

        rf_random = RandomizedSearchCV(estimator = rf,
                                        param_distributions = random_grid,
                                        n_iter = 100,
                                        cv = 3,
                                        random_state = 42,
                                        n_jobs = threads)
        return rf_random

    def tune(self, features_list, labels_list_numeric, testing_portion, grid_search, threads, rf):
        '''     
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        logging.info('Splitting data into test and training datasets')
        train_features, test_features, train_labels, test_labels \
            = train_test_split(features_list,
                               labels_list_numeric,
                               test_size = testing_portion,
                               random_state = 7)


        rf_random_model = self.random_search_cv(threads, rf)
        logging.info('Fitting model')
        rf_random_trained_model = rf_random_model.fit(train_features, train_labels)

        logging.info('Best parameters from random search cross validation:')
        for x,y in rf_random_trained_model.best_params_.items():
            logging.info("\t\t%s: %s" % (x, str(y)))
        if grid_search:
            rf_grid_model = self.grid_search_cv(rf_random_model, threads, rf)   
            logging.info('Fitting model')
            rf_grid_trained_model = rf_grid_model.fit(train_features, train_labels)

            logging.info('Best parameters from grid search cross validation:')
            for x,y in rf_grid_trained_model.best_params_.items():
                logging.info("\t\t%s: %s" % (x, str(y)))

            best_model = rf_grid_trained_model.best_estimator_
        else:
            best_model = rf_random_trained_model.best_estimator_

        return best_model, test_features, test_labels

    def do(self, input_matrix_path, groups_path, model_type,
           testing_portion, grid_search, threads, output_directory):
        '''
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        logging.info('Using %f%% of the input data for testing' % (testing_portion*100))

        if model_type == self.REGRESSOR:
            model = RandomForestRegressor()
        elif model_type == self.CLASSIFIER:
            model = RandomForestClassifier()
        else:
            raise Exception("Model type not recognised: %s" % (model_type))

        logging.info('Parsing inputs:')
        labels \
            = self.parse_groups_matrix(groups_path)
        features, attribute_list \
            = self.parse_input_matrix(input_matrix_path)

        labels_list, features_list = self.transpose(labels, features, attribute_list)

        labels_dict, labels_list_numeric = self._numerify(labels_list)
        labels_list_numeric     =   np.array(labels_list_numeric)
        features_list   =   np.array(features_list)
        logging.info('Splitting data into training and testing portions')

        logging.info("Tuning hyperparameters")
        rf, test_features, test_labels = self.tune(features_list,
                                                   labels_list_numeric,
                                                   testing_portion,
                                                   grid_search,
                                                   threads,
                                                   model)

        logging.info('Making predictions on test data:')
        predictions = rf.predict(test_features)
        errors = abs(predictions - test_labels)
        logging.info('\t\tMean Absolute Error: %f degrees' % (round(np.mean(errors), 2)))
        
        correctness = []
        for prediction, label in zip(np.round(predictions), test_labels):
            if prediction==label:
                correctness.append(1)
            else:
                correctness.append(0)

        accuracy = (sum(correctness)/float(len(correctness)))*100
        logging.info('\t\tAccuracy: %f%%' % (round(accuracy, 2)))
        self._write_importances(rf, attribute_list, output_directory)

        self._write_attribute_list(attribute_list, output_directory)

        logging.info("Preserving model")
        pickle.dump(rf, open(os.path.join(output_directory, self.MODEL_PICKLE) , 'wb'))

        logging.info("Preserving group labels")
        pickle.dump(labels_dict, open(os.path.join(output_directory, self.LABELS_DICT) , 'wb'))
