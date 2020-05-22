#!/usr/bin/env python3
# pylint: disable=line-too-long
"""
Various functions for generating machine learning models for genomes.
"""
# Imports
import pickle
import os
import logging
import numpy as np
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import train_test_split, RandomizedSearchCV, GridSearchCV
# Local
from enrichm.writer import Writer
from enrichm.parser import Parser

################################################################################

class GenerateModel:
    '''
    Functions to generate, train and optimise random forest machine learning models.
    '''
    def __init__(self):
        '''
        Inputs
        ------

        Outputs
        -------

        '''
        # Subparser names
        self.regressor = "regressor"
        self.classifier = "classifier"

        # Output file names
        self.attribute_importances = 'attribute_importances.tsv'
        self.model_pickle = "rf_model.pickle"
        self.labels_dict = "labels_dict.pickle"
        self.model_accuracy = "accuracy.tsv"

        # Headers
        self.attribute_importances_header = ['Variable', 'Importance']

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
        feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)

        # Print out the feature and importances
        logging.info('Writing attribute importances')

        output_lines = [self.attribute_importances_header]

        important_features = 0

        for pair in feature_importances:
            var, imp = pair
            output_lines.append([str(var), str(imp)])

            if imp > 0:
                important_features += 1

        logging.info('%i attributes found with an importance > 0', important_features)

        return output_lines

    def transpose(self, labels, features, attribute_list):
        '''

        Inputs
        ------

        Outputs
        -------

        '''
        labels_list = list()
        features_list = list()

        for col in list(labels.keys()):

            labels_list.append(labels[col])
            col_values = list()

            for attribute in attribute_list:
                col_values.append(features[col][attribute])

            features_list.append(col_values)

        return labels_list, np.array(features_list)

    def grid_search_cv(self, random_search_cv, threads, random_forest_model):
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

        param_grid = {'bootstrap':bootstrap,
                      'max_depth':[x for x in max_depth if x > 0],
                      'max_features':max_features,
                      'min_samples_leaf':[x for x in min_samples_leaf if x > 0],
                      'min_samples_split':[x for x in min_samples_split if x > 0],
                      'n_estimators':[x for x in n_estimators if x > 0]
                      }

        grid_search = GridSearchCV(estimator=random_forest_model,
                                   param_grid=param_grid,
                                   cv=3,
                                   n_jobs=threads)
        return grid_search

    def random_search_cv(self, threads, random_forest_model):
        '''
        Inputs
        ------

        Outputs
        -------

        '''
        logging.info('Generating model for random search cross validation')

        n_estimators = [int(x) for x in np.linspace(start=200, stop=2000, num=10)]
        max_features = ['auto', 'sqrt']
        max_depth = [int(x) for x in np.linspace(10, 110, num=11)]
        max_depth.append(None)
        min_samples_split = [2, 5, 10]
        min_samples_leaf = [1, 2, 4]
        bootstrap = [True, False]

        random_grid = {'n_estimators': n_estimators,
                       'max_features': max_features,
                       'max_depth': max_depth,
                       'min_samples_split': min_samples_split,
                       'min_samples_leaf': min_samples_leaf,
                       'bootstrap': bootstrap}

        rf_random = RandomizedSearchCV(estimator=random_forest_model,
                                       param_distributions=random_grid, n_iter=100, cv=3,
                                       random_state=42, n_jobs=threads)
        return rf_random

    def tune(self, features_list, labels_list_numeric, testing_portion, grid_search, threads,
             random_forest_model):
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
                               test_size=testing_portion,
                               random_state=7)
        best_params_list = list()

        rf_random_model = self.random_search_cv(threads, random_forest_model)
        logging.info('Fitting model')
        rf_random_trained_model = rf_random_model.fit(train_features, train_labels)

        logging.info('Best parameters from random search cross validation:')

        for parameter, value in rf_random_trained_model.best_params_.items():
            logging.info("\t\t%s: %s", parameter, str(value))

        if grid_search:
            rf_grid_model = self.grid_search_cv(rf_random_model, threads, random_forest_model)
            logging.info('Fitting model')
            rf_grid_trained_model = rf_grid_model.fit(train_features, train_labels)

            logging.info('Best parameters from grid search cross validation:')

            for parameter, value in rf_grid_trained_model.best_params_.items():
                logging.info("\t\t%s: %s", parameter, str(value))
                best_params_list.append([parameter, str(value)])

            best_model = rf_grid_trained_model.best_estimator_

        else:
            best_model = rf_random_trained_model.best_estimator_

        return best_model, test_features, test_labels, best_params_list

    def estimate_correctness(self, predictions, test_labels):
        '''
        Calculate the percentage of predictions match with their known classifications
        (test labels).

        :param predictions: List of predictions
        :type predictions: list
        :param test_labels: List of known annotations
        :type test_labels: list
        :returns: The percentage of prediction that match their known classifications
        :rtype: float
        '''
        correctness = list()

        for prediction, label in zip(np.round(predictions), test_labels):

            if prediction == label:
                correctness.append(1)
            else:
                correctness.append(0)

        accuracy = round((sum(correctness)/float(len(correctness)))*100, 2)

        return accuracy

    def generate_pipeline(self, input_matrix_path, groups_path, model_type,
                          testing_portion, grid_search, threads, output_directory):
        '''
        Inputs
        ------

        Outputs
        -------

        '''
        logging.info('Using %f%% of the input data for testing', testing_portion*100)

        if model_type == self.regressor:
            model = RandomForestRegressor()
        elif model_type == self.classifier:
            model = RandomForestClassifier()
        else:
            raise Exception("Model type not recognised: %s" % (model_type))

        logging.info('Parsing inputs:')
        labels, _, _ = Parser.parse_metadata_matrix(groups_path)
        features, _, attribute_list = Parser.parse_simple_matrix(input_matrix_path)
        labels_list, features_list = self.transpose(labels, features, attribute_list)
        labels_dict, labels_list_numeric = self.numerify(labels_list)

        logging.info("Tuning hyperparameters")
        random_forest_model, test_features, test_labels, best_params_list = self.tune(features_list, labels_list_numeric, testing_portion, grid_search, threads, model)

        logging.info('Making predictions on test data:')
        predictions = random_forest_model.predict(test_features)
        errors = abs(predictions - test_labels)

        logging.info('\t\tMean Absolute Error: %f degrees', round(np.mean(errors), 2))
        accuracy = self.estimate_correctness(predictions, test_labels)

        logging.info('\t\tAccuracy: %f%%', accuracy)
        best_params_list.append(["Accuracy", str(accuracy)])

        logging.info("Generating attribute importances")
        output_attribute_importances = self.get_importances(random_forest_model, attribute_list)
        Writer.write(best_params_list, os.path.join(output_directory, self.model_accuracy))

        logging.info("Generating model accuracy summary file")
        Writer.write(output_attribute_importances, os.path.join(output_directory, self.attribute_importances))

        logging.info("Preserving model")
        with open(os.path.join(output_directory, self.model_pickle), 'wb') as model_io:
            pickle.dump(random_forest_model, model_io)

        logging.info("Preserving group labels")
        with open(os.path.join(output_directory, self.labels_dict), 'wb') as labels_io:
            pickle.dump(labels_dict, labels_io)
