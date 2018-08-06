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
import subprocess
import pickle
import tempfile
from itertools import chain

################################################################################

from genome import Genome, AnnotationParser
from annotate import Annotate

###############################################################################

class Compare:

	def __init__(self, threads):
		self.threads = threads

	def _parse_pickles(self, enrichm_annotate_output):
		'''
		Opens the pickled genome objects from a previous run of enrichm 
		annotate

		Parameters
		----------
		enrichm_annotate_output - String. Output directory of a previous run 
								  of enrichm annotate (At lease version 0.0.7)
								  or above
		Outputs
		-------
		List of Genome objects
		'''	
		output_genome_list = list()

		genome_pickle_file_path \
			= os.path.join(enrichm_annotate_output, Annotate.GENOME_OBJ)

		for pickled_genome in os.listdir(genome_pickle_file_path):
			pickled_genome_path = os.path.join(genome_pickle_file_path, pickled_genome)
			logging.info('Parsing genome: %s' % (pickled_genome_path))
			output_genome_list.append(pickle.load(open(pickled_genome_path)))

		return output_genome_list

	def _parse_blast(self, result_path):
		result = dict()
		for line in open(result_path):
			qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore	\
				= line.strip().split('\t')
			genome, sequence_id = qseqid.split('~')
			hit_genome, hit_sequence_id = sseqid.split('~')
			if genome not in result:
				result[genome] = {sequence_id:{hit_genome:[hit_genome, hit_sequence_id, evalue]}}
			else:
				if sequence_id not in result[genome]:
					result[genome][sequence_id] = {hit_genome:[hit_genome, hit_sequence_id, evalue]}
				else:
					if hit_genome not in result[genome][sequence_id]:
						result[genome][sequence_id][hit_genome] = [hit_genome, hit_sequence_id, evalue]
					else:
						p_genome, p_hit, p_eval = result[genome][sequence_id][hit_genome]
						if float(evalue)<float(p_eval):
							result[genome][sequence_id][hit_genome] = [hit_genome, hit_sequence_id, evalue]
		return result

	def _self_blast(self, genome_list):
		with tempfile.NamedTemporaryFile() as fasta:

			for genome in genome_list:
				for sequence_name, sequence in genome.sequences.items():
					fasta.write('>%s~%s\n' % (genome.name, sequence.seqname))
					fasta.write('%s\n' % (sequence.seq))
			fasta.flush()

			with tempfile.NamedTemporaryFile() as dmnd:
				cmd = 'diamond makedb --quiet --in %s --db %s ' % (fasta.name, dmnd.name)
				
				logging.debug(cmd)
				subprocess.call(cmd, shell=True)

				with tempfile.NamedTemporaryFile() as result:

					cmd = 'diamond blastp --quiet -q %s --db %s -e 1e-05 --query-cover 70 --subject-cover 70 --id 0.3 -f 6 -o %s --threads %s' \
							% (fasta.name, dmnd.name, result.name, self.threads)
					logging.debug(cmd)
					subprocess.call(cmd, shell=True)
					results = self._parse_blast(result.name)
					
		return results
	
	def calc_synt(self, tot, proteins):
		pass


	def _synteny(self, genome_list):
		best_hits = self._self_blast(genome_list)
		tot_genomes = len(genome_list)
		for genome in genome_list:
			for position, sequence_name in genome.protein_ordered_dict.items():
				
				upstream, downstream = position-1, position+1
				if upstream in genome.protein_ordered_dict:
					upstream_id = genome.protein_ordered_dict[upstream]
					if upstream_id in best_hits[genome.name]:
						up_best_hits = best_hits[genome.name][downstream_id]
						upstream_synt = self.calc_synt(tot_genomes, up_best_hits)
					else:
						upstream_synt = 1
				else:
					upstream_synt = 1
				if downstream in genome.protein_ordered_dict:
					downstream_id = genome.protein_ordered_dict[downstream]
					if downstream_id in best_hits[genome.name]:
						down_best_hits = best_hits[genome.name][downstream_id]
						downstream_synt = self.calc_synt(tot_genomes, down_best_hits)
						
					else:
						downstream_synt = 1			
				else:
					downstream_synt = 1
	def _metadata_parser(self, metadata_path):
		metadata = dict()
		seen = set()
		for line in open(metadata_path):
			genome, group = line.strip().split('\t')
			if genome in seen:
				raise Exception("Duplicate entry in metadata file: %s " % genome)
			if group in metadata:
				metadata[group].add(genome)
			else:
				metadata[group] = set([genome])
			seen.add(genome)
		return metadata


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

	def do(self, enrichm_annotate_output, metadata_path):
		'''
		Parameters
		----------
		enrichm_annotate_output - String. Output directory of a previous run 
								  of enrichm annotate (At lease version 0.0.7)
								  or above
		'''

		metadata = self._metadata_parser(metadata_path)	

		logging.info('Parsing pickled genomes from previous enrichm run: %s' \
						% (enrichm_annotate_output))
		
		genome_list = self._parse_pickles(enrichm_annotate_output)
		genome_dict = {genome.name:genome for genome in genome_list}

		pan_genome_results = self._pan_genome(genome_dict, metadata)
		import IPython ; IPython.embed()
		#self._synteny(genome_list)
		
