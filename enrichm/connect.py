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
# System imports
import logging
import os
from collections import Counter
# Local imports
from enrichm.parse_annotate import ParseAnnotate
from enrichm.databases import Databases
from enrichm.module_description_parser import ModuleDescription

###############################################################################

class Connect(object):

	def __init__(self):
		d = Databases()
		
		self.m2def = d.m2def
		self.m2c = d.m2c
		self.c2m = dict()

		for module, compounds in self.m2c.items():	
			substrates = compounds[0]
			for substrate in substrates:
				if substrate in self.c2m:
					self.c2m[substrate].append(module)
				else:
					self.c2m[substrate] = [module]

		self.c = d.c
		self.m = d.m
		self.signature_modules = d.signature_modules
		self.output_file = 'linkages.tsv'
	
	def _update_with_custom_modules(self, custom_modules):
		custom_modules_dict = {line.split('\t')[0]:line.strip().split('\t')[1]
								for line in open(custom_modules)}
		self.m2def.update(custom_modules_dict)
		
		for key in custom_modules_dict.keys():
			self.m[key] = 'Custom'

	def _parse_genome_and_annotation_file_matrix(self, genome_and_annotation_file):
		genome_and_annotation_file_io = open(genome_and_annotation_file)
		headers = genome_and_annotation_file_io.readline().strip().split('\t')[1:]
		genome_to_annotation_sets = {genome_name:set() for genome_name in headers}

		for line in genome_and_annotation_file_io:
			sline = line.strip().split('\t')
			annotation, entries = sline[0], sline[1:]

			for genome_name, entry in zip(headers, entries):

				if float(entry) > 0.0:
					genome_to_annotation_sets[genome_name].add(annotation)

		return genome_to_annotation_sets

	def _write(self, lines, output_path):
		with open(output_path, 'wb') as out_io:

			for l in lines:
				out_io.write(('\t'.join(l) + '\n').encode())

	def do(self, annotate_output, metadata, custom_modules, cutoff, output_directory):
		pa = ParseAnnotate(annotate_output, 1)
		annotations = self._parse_genome_and_annotation_file_matrix(pa.ko)

		if custom_modules:
			self._update_with_custom_modules(custom_modules)

		genomes_list = list(annotations.keys())
		genomes_modules = dict()
		module_to_genome = dict()
		candidate_linkages = [["Genome 1",
					   		   "Genome 2",
					   		   "Compound",
					   		   "Compound description",
					   		   "Genome 1 module (producer)",
					   		   "Genome 2 module (consumer)"]]

		for module in self.m2def:

			for genome_name in genomes_list:
				kos = annotations[genome_name]

				if module not in self.signature_modules:
					path = ModuleDescription(self.m2def[module])
					num_all = path.num_steps()
					num_covered, _, _, _ = path.num_covered_steps(kos)
					perc_covered = num_covered / float(num_all)

					if perc_covered>= float(cutoff):
						if module in module_to_genome:
							module_to_genome[module].add(genome_name)
						else:
							module_to_genome[module] = set([genome_name])

						if genome_name in genomes_modules:
							genomes_modules[genome_name].add(module)
						else:
							genomes_modules[genome_name] = set([module])

		for genome_name, modules in genomes_modules.items():

			for module in modules:
				
				if module in self.m2c:
					produced_compounds = self.m2c[module][1]

					for compound in produced_compounds:

						if compound in self.c2m:
							possible_modules 	= self.c2m[compound]
							modules_not_covered = list()

							for possible_module in possible_modules:

								if possible_module not in modules:
									modules_not_covered.append(possible_module)
							
							for nc_module in modules_not_covered:
								
								if nc_module in module_to_genome:
									candidate_genomes = module_to_genome[nc_module]

									for genome in candidate_genomes:
										
										if genome != genome_name:
											encoded_modules = genomes_modules[genome]
											
											if module not in encoded_modules:
												candidate_linkages.append([genome_name,
																		   genome,
																		   compound,
																		   self.c[compound],
																		   module,
																		   nc_module])

		self._write(candidate_linkages, os.path.join(output_directory, self.output_file))

