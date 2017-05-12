#!/usr/bin/env python2
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
 
__author__ = "Joel Boyd"
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
 
###############################################################################

import logging
import os
import gzip

from genome import Genome

###############################################################################

class MatrixGenerator:
    
    COMPRESSED_SUFFIXES = set(['.gz', '.gzip'])
    REFERENCE_PATH = '/srv/db/uniprot/201607/KO.idmapping.dat.gz'
    OLD_REFERENCE_PATH='/srv/db/uniprot/uniref_20151020/idmapping.KO.dat.gz'
    MATRIX_SUFFIX  = '_matrix.tsv'
    UR100 = 'UniRef100_'
    
    KO      = 'KO_IDS.txt'
    PFAM    = 'PFAM_IDS.txt'
    TIGRFAM = 'TIGRFAM_IDS.txt'

    def __init__(self, annotation_type):        
        data_directory = os.path.join(os.path.split(os.path.realpath(__file__))[0], '../data/ids/')
        if annotation_type == self.KO:
            ids = [x.strip() for x in open(os.path.join(data_directory,self.KO))]
        elif annotation_type == self.PFAM:
            ids = [x.strip() for x in open(os.path.join(data_directory,self.PFAM))]
        elif annotation_type == self.TIGRFAM:
            ids = [x.strip() for x in open(os.path.join(data_directory,self.TIGRFAM))]

        #logging.info("Loading uniref to orthology information")
        #self.reference_dictionary = {}
        #for line in gzip.open(self.REFERENCE_PATH):
        #    protein_id, _, ko_id  = line.strip().split('\t')
        #    self.reference_dictionary[protein_id] = ko_id
        #self.possible_kos = list(set(self.reference_dictionary.values()))
    
    def _open(self, file):
        file_suffix = os.path.splitext(file)[-1]
        if file_suffix in self.COMPRESSED_SUFFIXES:
            return gzip.open(file)
        else:
            return open(file)
    
    def _clean(self, string):
        '''
        Strip off uniprot prefix if present.
        '''
        if string.startswith(self.UR100):
            return string.replace(self.UR100,'')
        else:
            return string
    
    def write_matrix(self, input_annotations, output_path):
        '''
        Writes a frequency matrix with of each annotation (rows) per sample (columns)

        Parameters
        ----------
        input_annotations   - list.
        output_path         - string. Path to file to which the results are written.
        '''
        logging.debug("    - Writing results to file: %s" % output_path)
        

        colnames = ['ID'] + [os.path.splitext(os.path.basename(filename))[0]
                             for filename in input_annotations.keys()]
        #



    def from_blast_results(self, blast_output_paths, evalue_cutoff, bitscore_cutoff, percent_id_cutoff, output_matrix_path):
        logging.info("Parsing %i blast output files" % \
                                    (len(blast_output_paths)))
        
        output_dictionary = {}
        for file in blast_output_paths:
            # Setting up sample column
            output_dictionary[file] = {}
            for ko in self.possible_kos:
                output_dictionary[file][ko]=0
            # Filling in column
            for line in self._open(file):
                sline = line.strip().split()
                
                if(float(sline[10]) <= evalue_cutoff and
                   float(sline[11]) >= bitscore_cutoff and
                   float(sline[2])  >= percent_id_cutoff and                   
                   float(sline[2])  >= percent_id_cutoff):
                        hit_id = self._clean(sline[1])
                        try:
                            hit_ko = self.reference_dictionary[hit_id]
                            output_dictionary[file][hit_ko] += 1
                        except:
                            raise Exception("UniProt ID %s not found in \
reference dictionary: %s. This probably means that an out-dated UniProt \
reference database was used. It is highly recommended that you re-run blast \
using the latest version of the UniProt database." \
                                             % (hit_id, self.REFERENCE_PATH))


        logging.info("Writing results to file: %s" % output_matrix_path)
        with open(output_matrix_path, 'w') as output_matrix_path_io:
            header = 'KO\t%s\n' % '\t'.join([os.path.basename(path) for 
                                             path in blast_output_paths])
            output_matrix_path_io.write(header)
            for ko in self.possible_kos:
                output_line = [ko]
                
                for file in blast_output_paths:
                    output_line.append( str(output_dictionary[file][ko]) )
                output_matrix_path_io.write("%s\n" % '\t'.join(output_line))
    
    def from_hmmsearch_results(self, hmmsearch_output_paths, output_matrix_path, evalue_cutoff, bitscore_cutoff, 
                               percent_aln_query_cutoff, percent_aln_reference_cutoff):
        '''
        Parse input hmmsearch files and generate a frequency matrix of annotations. 

        Parameters 
        ----------
        hmmsearch_output_paths          -
        output_matrix_path              -
        evalue_cutoff                   -
        bitscore_cutoff                 -
        percent_aln_query_cutoff        -       
        percent_aln_reference_cutoff    -

        Outputs
        -------

        '''
        logging.debug("    - Parsing %i hmmsearch output files" % len(hmmsearch_output_paths))
        
        output_dictionary = {}
        for file in hmmsearch_output_paths:
            genome = Genome(file)

            # Filling in column
            for line in open(file):
                if line.startswith('#'): continue
                
                seqname, _, tlen, query_name, accession, qlen, evalue, score, \
                bias, _, _, c_evalue, i_evalue, dom_score, dom_bias, hmm_from, \
                hmm_to, seq_from, seq_to, _, _, acc = line.strip().split()[:22]

                seq_list = [int(seq_from), int(seq_to)]
                hmm_list = [int(hmm_from), int(hmm_to)]

                perc_seq_aln = (max(seq_list)-min(seq_list))/float(tlen)
                perc_hmm_aln = (max(seq_list)-min(seq_list))/float(qlen)

                if(float(evalue)<=evalue_cutoff and
                   float(score)>=bitscore_cutoff and
                   perc_seq_aln>=percent_aln_query_cutoff and
                   perc_hmm_aln>=percent_aln_reference_cutoff):
                    if seqname in genome.sequences:
                        output_dictionary[file][seqname].add(query_name, evalue, range(min(seq_list), max(seq_list)))
                    else:
                        sequence_object = Sequence(seqname, tlen)
                        sequence_object.add(query_name, evalue, 
                                            range(min(seq_list), max(seq_list)))
                        output_dictionary[file][seqname]=sequence_object

            output_dictionary[file] = genome
        self.write_matrix(output_dictionary, output_matrix_path, )
        
        return output_dictionary