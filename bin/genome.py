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

import os
from Bio import SeqIO


class Genome:

	def __init__(self, path):
		self.path = path
		self.name = os.path.split(os.path.splitext(path)[0])[1]
		self.sequences = {}
		
		for protein in SeqIO.parse(path, 'fasta'):
		 	self.sequences[protein.name] = Sequence(protein.name, len(protein.seq))

class Sequence(Genome):

	def __init__(self, seqname, length):
		self.annotations = []	
		self.length = int(length)
		self.seqname = seqname

	def all_annotations(self):
		'''
		Returns a list of all annotations assigned to this sequence
		'''
		result = []
		for annotation in self.annotations:
			result.append(annotation.annotation)
		return result

	def seqdict(self):
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
		result = []

		# Build reference dictionary for sequence
		seq_dict = self.seqdict()

		# Find annotation for each position in the sequence
		for position in query_region:
			result.append(seq_dict[position])
		return result

	def add(self, annotation, evalue, region):
		'''
		Return annotations assigned to a list of positions within a sequence.

		Parameters
		----------
		annotation 	- string. Annotation to assign to the sequence region
		evalue 		- float. Evalue awarded for the given annotation
		region		- list. List of integers specifying the positions in the 
					  sequence to annotate
		'''
		new_annotation = Annotation(annotation, evalue, region)
		if len(self.annotations) > 0:
			for idx, previous_annotation in enumerate(self.annotations):
				if len(previous_annotation.region.intersection(new_annotation.region)) > 0:
					is_better = new_annotation.compare(previous_annotation)
					if is_better:
						self.annotations[idx] = new_annotation
		else:
			self.annotations.append(new_annotation)

	
class Annotation(Sequence):
	
	PFAM 	= 'PF'
	KO   	= 'K'
	TIGRFAM = 'TIGR'

	def __init__(self, annotation, evalue, region):
		self.annotation = annotation
		self.evalue 	= float(evalue)
		self.region		= set(region)

	def compare(self, other_annotation):
		'''
		Compares self with another annotation object. 

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
    KO      = 'KO_IDS.txt'
    PFAM    = 'PFAM_IDS.txt'
    TIGRFAM = 'TIGRFAM_IDS.txt'
    
	def __init__(annotation_type):        
        data_directory = os.path.join(os.path.split(os.path.realpath(__file__))[0], '../data/ids/')
        if annotation_type == self.KO:
            ids = [x.strip() for x in open(os.path.join(data_directory,self.KO))]
        elif annotation_type == self.PFAM:
            ids = [x.strip() for x in open(os.path.join(data_directory,self.PFAM))]
        elif annotation_type == self.TIGRFAM:
            ids = [x.strip() for x in open(os.path.join(data_directory,self.TIGRFAM))]
