#!/usr/bin/env python3
import logging
import os
from itertools import chain
from collections import Counter
from enrichm.databases import Databases
# Local
###############################################################################
class Writer:
    '''
    A collection of functions to write files in various formats.
    '''

    @staticmethod
    def write(output_lines_list, output_path):
        '''

        Parameters
        ----------
        output_lines_list   - List. A list of lists. Each sublist is a line,
                              where each entry is a column entry
        output_path         - String. Path to write output lines to.
        Output
        ------
        '''
        logging.info("Writing results to file: %s" % output_path)

        with open(output_path, 'w') as out_io:

            for output_line_list in output_lines_list:
                output_line_string = '\t'.join([str(column_entry) for column_entry in output_line_list]) + '\n'
                out_io.write(output_line_string)

            out_io.flush()
            out_io.close()

    @staticmethod
    def write_gff(genome, output_file):
        '''
        Writes a gff file for a Genome object

        Parameters
        ----------
        genome 			- Genome object
        output_file 	- string. file name to output results to.
        '''
        with open(output_file, 'w') as out_io:
            for sequence in genome.ordered_sequences():
                contig_id = '_'.join(sequence.seqname.split('_')[:-1])

                features = ['seq_id=%s' % sequence.seqname,
                                    'prodigal_id=%s' % sequence.prod_id,
                                    'partial=%s' % sequence.partial,
                                    'start_type=%s' % sequence.starttype,
                                    'rbs_motif=%s' % sequence.rbs_motif,
                                    'rbs_spacer=%s' % sequence.rbs_spacer,
                                    'gc=%s' % sequence.gc]

                if len(sequence.all_annotations()) > 0:
                    features.append('annotations=%s' % ','.join(sequence.all_annotations()))
                else:
                    features.append('annotations=hypothetical_protein')

                line = [contig_id,
                                    'prodigal',
                                    'CDS',
                                    sequence.startpos,
                                    sequence.finishpos,
                                    '.',
                                    ('-' if sequence.direction == '-1'
                                     else '+'),
                                    '0',
                                    ';'.join(features)]

                out_io.write('\t'.join(line) + '\n')

class MatrixGenerator:

    KO      = 'KO_IDS.txt'
    EC      = 'EC_IDS.txt'
    PFAM    = 'PFAM_IDS.txt'
    TIGRFAM = 'TIGRFAM_IDS.txt'
    CAZY = 'CAZY_IDS.txt'
    HYPOTHETICAL = 'HYPOTHETICAL'
    ORTHOLOG = 'ORTHOLOG'

    def __init__(self, annotation_type, clusters = None):
        '''
        Interpret which annotation type to write a matrix for.

        Parameters
        ----------
        annotation_type - String.
        '''
        self.annotation_type = annotation_type
        self.databases = Databases()
        if self.annotation_type == self.KO:
            self.annotation_list = [x.strip() for x in open(os.path.join(self.databases.IDS_DIR, self.KO))]

        elif self.annotation_type == self.EC:
            self.annotation_list = [x.strip() for x in open(os.path.join(self.databases.IDS_DIR, self.EC))]

        elif self.annotation_type == self.PFAM:
            self.annotation_list = [x.strip() for x in open(os.path.join(self.databases.IDS_DIR, self.PFAM))]

        elif self.annotation_type == self.TIGRFAM:
            self.annotation_list = [x.strip() for x in open(os.path.join(self.databases.IDS_DIR, self.TIGRFAM))]

        elif self.annotation_type == self.CAZY:
            self.annotation_list = [x.strip() for x in open(os.path.join(self.databases.IDS_DIR, self.CAZY))]

        elif self.annotation_type == self.HYPOTHETICAL:
            self.annotation_list = clusters

        elif self.annotation_type == self.ORTHOLOG:
            self.annotation_list = clusters

        else:
            raise Exception("Annotation type not found: %s" % (self.annotation_type))

    def write_matrix(self, genomes_list, count_domains, output_path):
        '''
        Writes a frequency matrix with of each annotation (rows) per sample (columns)

        Parameters
        ----------
        genomes_list        - list. List of Genome objects
        output_path         - string. Path to file to which the results are written.
        '''

        logging.info("    - Writing results to file: %s" % output_path)

        with open(output_path, 'w') as out_io:
            colnames = ['ID'] + [genome.name for genome in genomes_list]
            out_io.write('\t'.join(colnames) + '\n')

            if count_domains:
                genome_annotations = {genome.name:Counter(chain(*[sequence.all_annotations() for sequence in genome.sequences.values()]))
                                      for genome in genomes_list}
            else:
                genome_annotations = {genome.name:Counter(chain(*[set(sequence.all_annotations()) for sequence in genome.sequences.values()]))
                                      for genome in genomes_list}

            for annotation in self.annotation_list:
                output_line = [annotation]

                for genome in genomes_list:

                    if annotation in genome_annotations[genome.name]:
                        output_line.append(str(genome_annotations[genome.name][annotation]))

                    else:
                        output_line.append('0')

                out_io.write( '\t'.join(output_line) + '\n' )
