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
__copyright__   = "Copyright 2017"
__credits__     = ["Joel Boyd"]
__license__     = "GPL3"
__version__     = "0.0.7"
__maintainer__  = "Joel Boyd"
__email__       = "joel.boyd near uq.net.au"
__status__      = "Development"

###############################################################################
# Imports
import logging
import os
import numpy as np
from scipy import stats
from itertools import combinations
# Local
from enrichm.genome import Genome, AnnotationParser
###############################################################################

class Compare:

	def __init__(self):

		# Result types
		self.gc = "GC content"
		self.length = "Genome length"
		self.tests = "Tests"

	def _pan_genome(self, genome_dict, metadata):
		result = dict()
		for group, genomes in metadata.items():
			group_genomes = [genome_dict[genome] for genome in genomes]
			group_core_genome_size = len(set.intersection(*[genome.clusters for genome in group_genomes]))
			group_genome_lengths = [len(genome.sequences) for genome in group_genomes]
			group_genome_mean = float(sum(group_genome_lengths)) / float(len(group_genome_lengths))
			group_core_size = round((float(group_core_genome_size) / group_genome_mean)*100, 2)
			result[group] = [group_core_genome_size, group_core_size]
		return result

	def gather_stats(self, number_list):
		number_array = np.array(number_list)
		mean = number_array.mean()
		median = np.median(number_array)
		std = np.std(number_array)
		minimum = min(number_array)
		maximum = max(number_array)
		p90 = np.percentile(number_array, 90)
		p10 = np.percentile(number_array, 10)

		return [mean, median, std, minimum, maximum, p90, p10]

	def run_mannwhitneyu(self, feature, feature_list):
		
		test_results = dict()

		for combination in combinations(feature_list.keys(), 2):

			group_1 = np.array(feature_list[combination[0]])
			group_2 = np.array(feature_list[combination[1]])

			test_results[combination] \
				= stats.mannwhitneyu(group_1, group_2)
		return test_results

	def _compare_features(self, genome_dict, metadata):
		# GC, size, coding density
		keys = list(metadata.keys())+[self.tests]

		results = {self.gc:		{x:[] for x in keys},
				   self.length:	{x:[] for x in keys}}

		for group, genomes in metadata.items():
			group_genomes = [genome_dict[genome] for genome in genomes]
			group_lengths = list()
			group_gc = list()
			
			for genome in group_genomes:
				group_gc.append(genome.gc)
				group_lengths.append(genome.length)

			results[self.gc][group].append([group_gc, self.gather_stats(group_gc)]) 
			results[self.length][group].append([group_lengths, self.gather_stats(group_lengths)]) 

		for feature, feature_list in results.items():
			test_result = self.run_mannwhitneyu(feature, feature_list)
			results[self.tests].append([feature, test_result])
		
		return results
	
	def saturate(self, genomes_list, orthologs, output):
		result = dict()

		for group in range(1, len(genomes_list)+1):
			
			core_counts = list()
			accessory_counts = list()

			for combination in combinations(genomes_list, group):

				core_count = 0
				accessory_count = 0

				combination_length = len(combination)

				for ortholog in orthologs:
					hits = len([True for genome in combination if ortholog in genome.orthologs])
					if hits==combination_length:
						core_count+=1
					else:
						accessory_count+=1
				core_counts.append(core_count)
				accessory_counts.append(accessory_count)

			result[group] = {"core":self.gather_stats(core_counts), 'accessory':self.gather_stats(accessory_counts)}

		with open(output, 'w') as out_io:
			out_io.write('\t'.join(["Group size",
									"Core mean",
									"Core median",
									"Core standard deviation",
									"Core minimum",
									"Core maximum",
									"Core p90",
									"Core p10",
									"Accessory mean",
									"Accessory median",
									"Accessory standard deviation",
									"Accessory minimum",
									"Accessory maximum",
									"Accessory p90",
									"Accessory p10"]) +'\n')
			for group, stats in result.items():
				output_line = [str(group)]
				output_line += [str(x) for x in stats["core"]]
				output_line += [str(x) for x in stats["accessory"]]
				out_io.write('\t'.join(output_line) +'\n')

	def write_pan_genome_results(self, results, output_path):
		header = ['Genome', "Core genome size", "Percent of average genome size"]
		with open(output_path, 'w') as out_io:
			out_io.write('\t'.join(header) + '\n')
			for genome, results_list in results.items():
				output_line = [genome] + [str(x) for x in results_list]
				out_io.write('\t'.join(output_line) + '\n')

	def do(self, genome_list, metadata, output_directory):
		'''
		Parameters
		----------
		genome_objects - List.
		'''
		genome_dict = {genome.name:genome for genome in genome_list}
		ortholog_list = set.union(*[genome.orthologs for genome in genome_list])

		for group_name, genome_list in metadata.items():	
			logging.info('Generating pan genome saturation curve for group: %s' % group_name)
			group_list = [genome_dict[g] for g in genome_list]

			output = os.path.join(output_directory, "%s_saturation.tsv" % (group_name))
			self.saturate(group_list, ortholog_list, output)

		pan_genome_results = self._pan_genome(genome_dict, metadata)
		self.write_pan_genome_results(pan_genome_results, os.path.join(output_directory, "core_genome_size.tsv"))
		#if hasattr(list(genome_dict.values())[0], "gc"):
		#	compare_features_results = self._compare_features(genome_dict, metadata)
		#	self.write_compare_features_results(compare_features_results)
		
