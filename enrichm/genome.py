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

from enrichm.databases import Databases
from enrichm.sequence_io import SequenceIO
import logging
import os
__author__ = "Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.7"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

###############################################################################

class Genome:
	'''
	A genome object which collects all the attirbutes of an imput genome,
	including protein sequences and their annotations
	'''
	def __init__(self, light, path, nucl):
		seqio = SequenceIO()
		self.clusters = set()
		self.orthologs = set()
		self.protein_ordered_dict = dict()
		self.sequences = dict()
		self.cluster_dict = dict()
		self.ortholog_dict = dict()
		self.path = path
		self.name = os.path.split(os.path.splitext(path)[0])[1]
		
		if light == False:
		
			if nucl is not None:
				self.nucl 	= nucl
				self.length = 0
				gc_list 	= 0.0

				for description, sequence in seqio.each(open(nucl)):
					self.length += len(str(sequence))
					gc_list 	+= (str(sequence).count('G') + str(sequence).count('C'))
				
				self.gc = round((gc_list/float(self.length))*100, 2)
		
			for protein_count, (description, sequence) in enumerate(seqio.each(open(path))):
				name = description.partition(' ')[0]
				sequence = Sequence(description, sequence)
				self.sequences[name] = sequence
				self.protein_ordered_dict[protein_count] = name
		
		else:
			for protein_count, (description, _) in enumerate(seqio.each(open(path))):
				name = description.partition(' ')[0]
				sequence = Sequence(description)
				self.sequences[name] = sequence
				self.protein_ordered_dict[protein_count] = name

	def add(self, annotations, evalue_cutoff, bitscore_cutoff, 
         percent_aln_query_cutoff, percent_aln_reference_cutoff, specific_cutoffs,
		    annotation_type, ref_ids):
		'''
		Adds a series of annotations to the proteins within a genome.

		Parameters
		----------
		annotations						- String. Path to file containing either blast or
										  domtblout hmmsearch, or mmseqs cluster results. 
		evalue_cutoff					- Float. E-value threshold for annotations.
		bitscore_cutoff					- Float. Bit score threshold for annotations.
		percent_aln_query_cutoff		- Float. Threshold for the percent of the query 
										  that must be aligned to consider the annotation.
		percent_aln_reference_cutoff	- Float. Threshold for the percent of the reference
										  that must be aligned to consider the annotation.
		annotation_type					- String. Either 'KO', 'TIGRFAM', 'PFAM', 	
										  'HYPOTHETICAL' or 'COG'
		'''
		# Load up annotation parser, and tell it what annotation type to expect
		ap = AnnotationParser(annotation_type)
		
		# If annotation type is a hmmsearch result
		if(annotation_type == AnnotationParser.HMMPARSER):
			# Set up an iterator to produce the results
			logging.debug("    - Parsing hmmsearch chunk")

			if(ref_ids == AnnotationParser.PFAM or
			   ref_ids == AnnotationParser.KO_HMM):
				percent_aln_query_cutoff = 0.0
				percent_aln_reference_cutoff = 0.0

			iterator = ap.from_hmmsearch_results(annotations, evalue_cutoff,
												 bitscore_cutoff, percent_aln_query_cutoff, 
												 percent_aln_reference_cutoff,
                                        specific_cutoffs,
                                        (True if ref_ids == AnnotationParser.KO_HMM else False))
			
			if ref_ids == AnnotationParser.PFAM:
				self.pfam_dict = dict()
				refdict = self.pfam_dict

			elif ref_ids == AnnotationParser.KO_HMM:
				self.ko_dict = dict()
				refdict = self.ko_dict
			
			elif ref_ids == AnnotationParser.TIGRFAM:
				self.tigrfam_dict = dict()
				refdict = self.tigrfam_dict
			
			elif ref_ids == AnnotationParser.CAZY:
				self.cazy_dict = dict()
				refdict = self.cazy_dict

		# If annotation type is a blast result
		
		elif(annotation_type == AnnotationParser.BLASTPARSER):
			# Set up an iterator to produce the results
			logging.debug("    - Parsing blast chunk")
			iterator = ap.from_blast_results(annotations, evalue_cutoff, 
											 bitscore_cutoff, percent_aln_query_cutoff)

			if ref_ids == AnnotationParser.KO:
				self.ko_dict = dict()
				refdict = self.ko_dict

			elif ref_ids == AnnotationParser.EC:
				self.ec_dict = dict()
				refdict = self.ec_dict

		for seqname, annotations, evalue, annotation_range in iterator:
			self.sequences[seqname].add(annotations, evalue, annotation_range, ref_ids)
			for annotation in annotations:

				if annotation in refdict:
					refdict[annotation].append(seqname)
				
				else:
					refdict[annotation]=[seqname]
			
	def count(self, annotation, type):
		'''
		
		Parameters
		----------
		annotation - String. An annotation ID to return a frequency for
		
		Output
		------
		The number of times this annotation was encountered in the genome
		'''

		if type==AnnotationParser.HYPOTHETICAL:
			reference_dict = self.cluster_dict
		
		elif type==AnnotationParser.KO:
			reference_dict = self.ko_dict
		
		elif type==AnnotationParser.EC:
			reference_dict = self.ec_dict
		
		elif type==AnnotationParser.PFAM:
			reference_dict = self.pfam_dict
		
		elif type==AnnotationParser.TIGRFAM:
			reference_dict = self.tigrfam_dict
		
		if annotation in reference_dict:
			return len(reference_dict[annotation])
		
		else:
			return 0

	def ordered_sequences(self):
		'''
		Iterator that yields that all protein coding Sequence objects in a genome in order.
		'''
		for sequence_id in sorted(self.protein_ordered_dict.keys()):
			
			yield self.sequences[self.protein_ordered_dict[sequence_id]]

	def add_cluster(self, sequence_id, cluster_id):
		'''
		Add hypothetical clusters to genome sequences
		
		Inputs
		------
		cluster_list - array. list where each entry is a list
					   of n = 2: sequence_name, cluster_name
		'''
	
		if cluster_id in self.cluster_dict:
			self.cluster_dict[cluster_id].append(sequence_id)
		
		else:
			self.cluster_dict[cluster_id] = [sequence_id]
		
		annotation = Annotation(cluster_id, 0, [-1,0], AnnotationParser.HYPOTHETICAL)
		
		self.sequences[sequence_id].cluster = cluster_id
		self.sequences[sequence_id].annotations.append(annotation)
		
		self.clusters.add(cluster_id)
	
	def add_ortholog(self, sequence_id, ortholog_id):
		'''
		Add hypothetical ortholog to genome sequences
		
		Inputs
		------
		cluster_list - array. list where each entry is a list
					   of n = 2: sequence_name, cluster_name
		'''
	
		if ortholog_id in self.ortholog_dict:
			self.ortholog_dict[ortholog_id].append(sequence_id)
		
		else:
			self.ortholog_dict[ortholog_id] = [sequence_id]
		
		annotation = Annotation(ortholog_id, 0, [-1,0], AnnotationParser.ORTHOLOG)
		
		self.sequences[sequence_id].ortholog = ortholog_id
		self.sequences[sequence_id].annotations.append(annotation)
		self.orthologs.add(ortholog_id)

class Sequence(Genome):
	'''
	Sequence object which collects all attributes of a sequence including its length,
	and annotations. Can compare current annotation with new annotaitons.
	'''
	def __init__(self, description, sequence=None):
		self.annotations = list()
		line_split = description.split(' # ')
		
		if sequence:
			self.seq = str(sequence)
			self.length = int(len(sequence))	

		if len(line_split)==5:
			self.seqname, self.startpos, self.finishpos, self.direction, stats \
								= line_split
			self.prod_id, self.partial, self.starttype, self.rbs_motif, self.rbs_spacer,  self.gc \
								= [x.split('=')[1] for x in stats.split(';')]
		
		elif len(line_split)==1:
			self.seqname = description
			
		else:
			raise Exception("Malformatted sequence file!")

	def all_annotations(self):
		'''
		Returns a list of all annotations assigned to this sequence
		'''
		result = []
		
		for annotation in self.annotations:
			result.append(annotation.annotation)
		
		return result

	def seqdict(self):
		'''
		Output
		------
		A dictionary where each entry is a position in the sequence,
		and values are the annotation at that position. This is important
		in particular for pfam annotations as a protein can have >1 
		domains.
		'''
		seq_dict = {x:None for x in range(self.length)}
		
		for annotation in self.annotations:
		
			for position in annotation.region:
				seq_dict[position] = annotation.annotation
		
		return seq_dict

	def what(self, query_region):
		'''
		Return annotations assigned to a list of positions within a sequence.

		Parameters
		----------
		region 		- list. List of integers specifying the positions in the 
					  sequence to return annotations for
		Outputs
		-------
		Returns a list of equal length to region, containing the annotation of
		each position.
		'''
		result = list()

		# Build reference dictionary for sequence
		seq_dict = self.seqdict()

		# Find annotation for each position in the sequence
		for position in query_region:
			result.append(seq_dict[position])
		
		return result

	def add(self, annotations, evalue, region, annotation_type):
		'''
		Return annotations assigned to a list of positions within a sequence.

		Parameters
		----------
		annotation 	- string. Annotation to assign to the sequence region
		evalue 		- float. Evalue awarded for the given annotation
		region		- list. List of integers specifying the positions in the 
					  sequence to annotate
		'''

		new_annotations = [Annotation(annotation, evalue, region, annotation_type) for annotation in annotations]
		annotation_list = [annot for annot in self.annotations if annot.type == new_annotations[0].type]
		
		if len(annotation_list) > 0:
			
			to_remove 	= list()
			to_check 	= annotation_list

			for new_annotation in new_annotations:
				overlap 	= [previous_annotation for previous_annotation in to_check
							   if len(previous_annotation.region.intersection(new_annotation.region)) > 0]
				
				if annotation_type==AnnotationParser.PFAM:
					self.annotations.append(new_annotation)
				
				else:
					
					if len(overlap)>0:
						
						for overlapping_previous_annotation in overlap:
							
							if new_annotation.compare(overlapping_previous_annotation):
								to_remove.append(overlapping_previous_annotation)

						if len(to_remove)>0:
							self.annotations = [annotation for annotation in self.annotations 
												if annotation not in to_remove]
							self.annotations.append(new_annotation)
					
					else:
							self.annotations.append(new_annotation)

		else:
			
			for new_annotation in new_annotations:
				self.annotations.append(new_annotation)
	
class Annotation(Sequence):
	'''
	Annotation object that collects all attributes associated with a given
	annotation, like where it is in the sequence, its e-value and bit score,
	and what type of annotation it is.
	'''
	def __init__(self, annotation, evalue, region, annotation_type):
		self.annotation = annotation
		self.evalue 	= float(evalue)
		self.region		= set(region)
		self.type   	= annotation_type

	def compare(self, other_annotation):
		'''
		Compares evalue of current annotation with the evaulue of another annotation object. 			

		Parameters
		----------
		other_annotation: Annotation object.
		
		Output
		------
		Returns True if self is the better annotation or False if it isn't
		'''
		if self.evalue < other_annotation.evalue:
			return True
		
		else:
			return False

class AnnotationParser:
	'''
	Annotation parser class contains fucntions to parse hmmsearch domtblouts and blast results 
	currently for: KO, PFAM and TIGRFAM. COG to come	
	'''
	KO      		= 'KO_IDS.txt'
	KO_HMM 			= 'KO_IDS.txt'
	EC				= 'EC_IDS.txt'
	PFAM    		= 'PFAM_IDS.txt'
	TIGRFAM 		= 'TIGRFAM_IDS.txt'
	CAZY      		= 'CAZY_IDS.txt'
	HYPOTHETICAL 	= 'HYPOTHETICAL.txt'
	ORTHOLOG 		= 'ORTHOLOG.txt'
	HMMPARSER 		= 'hmm'
	BLASTPARSER 	= 'blast'

	def __init__(self, annotation_type):        
		pass
	
	def from_blast_results(self,
						   blast_output,
						   evalue_cutoff,
						   bitscore_cutoff, 
						   percent_id_cutoff):
		'''
		Parse blast output in tab format.	

		Parameters
		----------
		blast_output 		- String. Path to blast output file containing results.
							  Must be in blast output format 6
		evalue_cutoff		- Float. E-value threshold for annotations.
		bitscore_cutoff		- Float. Bit score threshold for annotations.
		percent_id_cutoff 	- Float. Percent identity threshold for annotations.
		
		Yields
		------
		A sequence name, annotation, E-value and region hit for every annotation result in 
		blast_output that pass a series of specified cutoffs
		'''

		for sline in blast_output:
			# Parse out important information from each line in blast output
			evalue   = sline[10]
			bit      = sline[11]
			perc_id  = sline[2]
			seq_list = [int(x) for x in sline[6:8]]

			# If the annotation passes the specified cutoffs
			if(float(evalue) <= evalue_cutoff and
				float(bit) >= bitscore_cutoff and
				float(perc_id) >= percent_id_cutoff):
					seqname = sline[0].split('~')[1]
					annotation = sline[1].split('~')[1].split('+')
					
					yield seqname, annotation, evalue, range(min(seq_list), max(seq_list))

	def from_hmmsearch_results(self,
							   hmmsearch_output_path,
							   evalue_cutoff,
							   bitscore_cutoff, 
    						   percent_aln_query_cutoff,
    						   percent_aln_reference_cutoff,
                            	specific_cutoffs,
							   acc = False):
		'''
		Parse input hmmsearch file

		Parameters 
		----------
		hmmsearch_output_path           - String. Path to domtblout file containing 
										  hmmsearch results.
		evalue_cutoff                   - Float. E-value threshold for annotations.
		bitscore_cutoff                 - Float. Bit score threshold for annotations.
		percent_aln_query_cutoff        - Float. Threshold for the percent of the query 
										  that must be aligned to consider the annotation.
		percent_aln_reference_cutoff    - Float. Threshold for the percent of the reference 	
										  that must be aligned to consider the annotation.
		Yields
		------
		A sequence name, accession, E-value and region hit for every annottation result in 
		blast_output_path that pass a set of cutoffs
		'''
		
		# Filling in column
		for line in open(hmmsearch_output_path):
			
			# Skip headers
			if line.startswith('#'): continue
				
			# Parse HMMsearch line. '_'s represent unimportant entries. Line
			# is trimmed using [:22] to remove sequence description
			seqname, _, tlen, ko_hmm, accession, qlen, _, score, \
			_, _, _, _, i_evalue, dom_score, _, _, \
			_, seq_from, seq_to, _, _, _ = line.strip().split()[:22]				

			# Determine sequence and HMM spans
			seq_list = [int(seq_from), int(seq_to)]

			# Calculate percent of the query and reference aligned to each-other. 
			perc_seq_aln = (max(seq_list)-min(seq_list))/float(tlen)
			perc_hmm_aln = (max(seq_list)-min(seq_list))/float(qlen)
			
			if acc:
				accession = ko_hmm

			# If the annotation passes the specified cutoffs
			if specific_cutoffs:
				if ko_hmm in specific_cutoffs:
					if specific_cutoffs[ko_hmm][1]=='full':
						if float(score) < specific_cutoffs[ko_hmm][0]:
							continue
					elif specific_cutoffs[ko_hmm][1] == 'domain':
						if float(dom_score) < specific_cutoffs[ko_hmm][0]:
							continue

			if(float(i_evalue) <= evalue_cutoff and
									float(score) >= bitscore_cutoff and
									perc_seq_aln >= percent_aln_query_cutoff and
									perc_hmm_aln >= percent_aln_reference_cutoff):
				yield seqname, [accession], i_evalue, range(min(seq_list), max(seq_list))
