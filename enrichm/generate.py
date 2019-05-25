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
from enrichm.writer import Writer
from enrichm.parser import Parser

################################################################################

class GenerateModel:

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
        self.MODEL_PICKLE = "rf_model.pickle"
        self.LABELS_DICT = "labels_dict.pickle"
        self.MODEL_ACCURACY = "accuracy.tsv"
        # Headers
        self.ATTRIBUTE_IMPORTANCES_HEADER = ['Variable', 'Importance']

    def numerify(self, input_list):
        '''

        Inputs
        ------

        Outputs
        -------

        '''
        counter = 0
        output_dictionary = dict()
        output_list = list()

        for group in input_list:
            group = group.pop()

            if group not in output_dictionary:
                output_dictionary[group] = counter
                counter += 1

            output_list.append(output_dictionary[group])

        output_dictionary = {item:key for key, item in output_dictionary.items()}

        return output_dictionary, np.array(output_list)

    def get_importances(self, model, attribute_list):
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

        output_lines = [self.ATTRIBUTE_IMPORTANCES_HEADER]

        for pair in feature_importances:
            var, imp = pair
            output_lines.append([str(var), str(imp)])

        return output_lines

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

        return labels_list, np.array(features_list)

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
        best_params_list = list()

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
                best_params_list.append([x, str(y)])

            best_model = rf_grid_trained_model.best_estimator_

        else:
            best_model = rf_random_trained_model.best_estimator_

        return best_model, test_features, test_labels, best_params_list

    def estimate_correctness(self, predictions, test_labels):
        correctness = list()

        for prediction, label in zip(np.round(predictions), test_labels):

            if prediction==label:
                correctness.append(1)
            else:
                correctness.append(0)

        accuracy = round( (sum(correctness)/float(len(correctness)))*100, 2)

        return accuracy

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
        labels, _, _ = Parser.parse_metadata_matrix(groups_path)
        features, _, attribute_list = Parser.parse_simple_matrix(input_matrix_path)
        labels_list, features_list = self.transpose(labels, features, attribute_list)
        labels_dict, labels_list_numeric = self.numerify(labels_list)

        logging.info("Tuning hyperparameters")
        rf, test_features, test_labels, best_params_list = self.tune(features_list, labels_list_numeric, testing_portion, grid_search, threads, model)

        logging.info('Making predictions on test data:')
        predictions = rf.predict(test_features)
        errors = abs(predictions - test_labels)

        logging.info('\t\tMean Absolute Error: %f degrees' % (round(np.mean(errors), 2)))
        accuracy = self.estimate_correctness(predictions, test_labels)

        logging.info('\t\tAccuracy: %f%%' % (accuracy))
        best_params_list.append( ["Accuracy", str(accuracy)] )

        logging.info("Generating attribute importances")
        output_attribute_importances = self.get_importances(rf, attribute_list)
        Writer.write(best_params_list, os.path.join(output_directory, self.MODEL_ACCURACY))

        logging.info("Generating model accuracy summary file")
        Writer.write(output_attribute_importances, os.path.join(output_directory, self.ATTRIBUTE_IMPORTANCES))

        logging.info("Preserving model")
        with open(os.path.join(output_directory, self.MODEL_PICKLE) , 'wb') as model_io:
            pickle.dump(rf, model_io)

        logging.info("Preserving group labels")
        with open(os.path.join(output_directory, self.LABELS_DICT) , 'wb') as labels_io:
            pickle.dump(labels_dict, labels_io)
