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
from network_analyzer import NetworkAnalyser

###############################################################################

class Preparer:
    
    COMPRESSED_SUFFIXES = set(['.gz', '.gzip'])
    REFERENCE_PATH = '/srv/db/uniprot/201607/KO.idmapping.dat.gz'
    OLD_REFERENCE_PATH='/srv/db/uniprot/uniref_20151020/idmapping.KO.dat.gz'
    MATRIX_SUFFIX  = '_matrix.tsv'
    UR100 = 'UniRef100_'
    
    def __init__(self):        
        logging.info("Loading uniref to orthology information")
        self.reference_dictionary = {}
        for line in gzip.open(self.REFERENCE_PATH):
            protein_id, _, ko_id  = line.strip().split('\t')
            self.reference_dictionary[protein_id] = ko_id
        logging.info('Done!')
        self.possible_kos = list(set(self.reference_dictionary.values()))
    
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
    
    def blast_output(self, 
                     blast_output_paths, 
                     evalue_cutoff, 
                     bitscore_cutoff, 
                     percent_id_cutoff, 
                     output_matrix_path):
        # TODO: select only best hit.
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
    
    def hmmsearch_output(self, 
                         hmmsearch_output_paths, 
                         evalue_cutoff, 
                         bitscore_cutoff, 
                         percent_aln_query_cutoff,
                         percent_aln_reference_cutoff,
                         output_matrix_path):
        logging.info("Parsing %i hmmsearch output files" % \
                                    (len(hmmsearch_output_paths)))
        
        output_dictionary = {}
        seq_dictionary    = {}
        
        for file in hmmsearch_output_paths:
            # Setting up sample column
            output_dictionary[file] = {}
            for ko in self.possible_kos:
                output_dictionary[file][ko]=0
            # Filling in column
            for line in open(file):
                if line.startswith('#'): continue
                seqname, _, tlen, query_name, accession, qlen, evalue, score, \
                bias, _, _, c_evalue, i_evalue, dom_score, dom_bias, hmm_from, \
                hmm_to, seq_from, seq_to, _, _, acc = line.strip().split()[:22]

                seq_list = [float(seq_from), float(seq_to)]
                hmm_list = [float(hmm_from), float(hmm_to)]
                perc_seq_aln = (max(seq_list)-min(seq_list))/float(tlen)
                perc_hmm_aln = (max(seq_list)-min(seq_list))/float(qlen)
                
                if(float(evalue)<=evalue_cutoff and
                   float(score)>=bitscore_cutoff and
                   perc_seq_aln>=percent_aln_query_cutoff and
                   perc_hmm_aln>=percent_aln_reference_cutoff):
                    if seqname in seq_dictionary:
                        if float(evalue)< seq_dictionary[seqname][1]:
                            seq_dictionary[seqname]=[query_name, float(evalue)]
                    else:
                        seq_dictionary[seqname]=[query_name, float(evalue)]
            for sequence_name, (annoatation, evalue) in \
                                                    output_dictionary.items():
                output_dictionary[file][annotation]+=1
        
        logging.info("Writing results to file: %s" % output_matrix_path)
        with open(output_matrix_path, 'w') as output_matrix_path_io:
            header = 'KO\t%s\n' % '\t'.join([os.path.basename(path) for 
                                             path in hmmsearch_output_paths])
            output_matrix_path_io.write(header)
            for ko in self.possible_kos:
                output_line = [ko]
                
                for file in hmmsearch_output_paths:
                    output_line.append( str(output_dictionary[file][ko]) )
                output_matrix_path_io.write("%s\n" % '\t'.join(output_line))
                
    def main(self, args):
        output_path = args.output_prefix + self.MATRIX_SUFFIX 
        if args.subparser_name == NetworkAnalyser.MATRIX:
            if args.blast_outputs:
                if(args.percent_aln_query_cutoff or 
                   args.percent_aln_reference_cutoff):
                    logging.warning("Cannot calculate the percent of the query \
aligned to the reference and vice versa with a blast output.")

                self.blast_output(args.blast_outputs,
                                  args.evalue_cutoff,
                                  args.bitscore_cutoff,
                                  args.percent_id_cutoff,
                                  output_path)
            elif args.hmmsearch_outputs:
                if(args.percent_id_cutoff):
                    logging.warning("Cannot calculate the percent identity \
with a     hmmsearch.")
                self.hmmsearch_output(args.hmmsearch_outputs,
                                      args.evalue_cutoff,
                                      args.bitscore_cutoff,
                                      args.percent_aln_query_cutoff,
                                      args.percent_aln_reference_cutoff,
                                      output_path)