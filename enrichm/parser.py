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

import os
import pickle
import multiprocessing as mp
from enrichm.databases import Databases
from enrichm.annotate import Annotate

################################################################################

def parse_genomes(path):
    genome = pickle.load(open(path, 'rb'))
    return genome

################################################################################

class Parser:
    '''
    A collection of functions to parse files in various formats.
    '''
    @staticmethod
    def parse_genome_and_annotation_file_lf(genome_and_annotation_file):
        genome_to_annotation_sets = dict()

        for line in open(genome_and_annotation_file):

            try:
                genome, annotation = line.strip().split("\t")

            except:
                raise Exception("Input genomes/annotation file error on %s" % line)

            if genome not in genome_to_annotation_sets:
                genome_to_annotation_sets[genome] = set()

            genome_to_annotation_sets[genome].add(annotation)

        return genome_to_annotation_sets

    @staticmethod
    def parse_genome_and_annotation_file_matrix(genome_and_annotation_matrix):
        genome_and_annotation_matrix_io = open(genome_and_annotation_matrix)
        headers = genome_and_annotation_matrix_io.readline().strip().split('\t')[1:]
        genome_to_annotation_sets = {genome_name:set() for genome_name in headers}

        for line in genome_and_annotation_matrix_io:
            sline = line.strip().split('\t')
            annotation, entries = sline[0], sline[1:]

            for genome_name, entry in zip(headers, entries):

                if float(entry) > 0:
                    genome_to_annotation_sets[genome_name].add(annotation)

        return genome_to_annotation_sets

    @staticmethod
    def parse_taxonomy(taxonomy_path):

        output_taxonomy_dictionary = dict()

        for line in open(taxonomy_path):
            genome, taxonomy_string = line.strip().split('\t')
            output_taxonomy_dictionary[genome] = taxonomy_string.split(';')

        return output_taxonomy_dictionary

    @staticmethod
    def parse_simple_matrix(matrix):
        with open(matrix) as matrix_io:
            colnames = matrix_io.readline().strip().split('\t')[1:]
            rownames = list()
            output_dict = {colname:dict() for colname in colnames}

            for line in matrix_io:
                sline = line.strip().split('\t')
                rowname, content = sline[0], sline[1:]

                if rowname not in rownames:
                    rownames.append(rowname)

                for key, value in zip(colnames, content):
                    output_dict[key][rowname] = float(value)

        return output_dict, colnames, rownames

    @staticmethod
    def parse_metadata_matrix(matrix_path):
        '''
        Parameters
        ----------
        matrix_path : String. Path to file containing a matrix of genome rownames.
        '''
        cols_to_rows = dict()
        nr_values = set()
        attribute_dict = dict()

        with open(matrix_path) as matrix_file_io:

            for line in matrix_file_io:
                rowname, entry = line.strip().split('\t')
                nr_values.add(entry)

                if entry in attribute_dict:
                    attribute_dict[entry].add(rowname)
                else:
                    attribute_dict[entry] = set([rowname])

                if rowname not in cols_to_rows:
                    cols_to_rows[rowname] = set([entry])
                else:
                    cols_to_rows[rowname].add(entry)

        return cols_to_rows, nr_values, attribute_dict

    @staticmethod
    def parse_single_column_text_file(text_file):
        entries = set()

        for entry in open(text_file):
            entry = entry.strip()
            entries.add(entry)

        return entries

    @staticmethod
    def filter_large_matrix(columns, matrix):
        '''
        description

        Inputs
        ------

        Outputs
        -------

        '''
        columns = list(columns)
        matrix_io = open(matrix)
        header = matrix_io.readline().strip().split('\t')

        indexes = list()
        include = list()

        for column in columns:

            if column in header:
                indexes.append(header.index(column))
                include.append(column)

        columns = include

        output_dict = {column:dict() for column in columns}

        for row in matrix_io:
            srow = row.strip().split()
            annotation = srow[0]

            for column, index in zip(columns, indexes):
                count = int(srow[index])

                if count > 0:
                    output_dict[column][annotation] = int(srow[index])

        return output_dict, columns

    @staticmethod
    def parse_tpm_values(tpm_values):
        from enrichm.databases import Databases
        k2r = Databases().k2r()

        output_dict = dict()
        annotation_types = set()
        genome_types = set()

        tpm_values_io = open(tpm_values, 'rb')
        tpm_values_io.readline()

        for line in tpm_values_io:
            gene, _, _, _, _, _, _, _, _, _, tpm, \
            _, _, annotation, sample = line.strip().split(b'\t')
            annotation_list = annotation.split(b',')
            tpm = float(tpm)
            genome = '_'.join(str(gene, "utf-8").split('_')[:2]) # temporary
            genome_types.add(genome)

            if sample not in output_dict:
                output_dict[sample] = dict()

            if genome not in output_dict[sample]:
                output_dict[sample][genome] = dict()

            for annotation_type in annotation_list:

                if str(annotation_type, "utf-8") in k2r:
                    reactions = k2r[str(annotation_type, "utf-8")]

                    for reaction in reactions:
                        reaction = str.encode(reaction)

                        if reaction not in output_dict[sample][genome]:
                            output_dict[sample][genome][reaction] = 0.0
                            annotation_types.add(reaction)

                        output_dict[sample][genome][reaction] += tpm
        return output_dict, annotation_types, genome_types

    @staticmethod
    def parse_enrichment_output(enrichment_output):
        fisher_results = dict()

        for file in os.listdir(enrichment_output):

            if file.endswith("fisher.tsv"):
                file = os.path.join(enrichment_output, file)
                file_io = open(file)
                file_io.readline()

                for line in file_io:
                    split_line = line.strip().split('\t')

                    if len(fisher_results) == 0:
                        fisher_results[split_line[1]] = list()
                        fisher_results[split_line[2]] = list()

                    if float(split_line[-2])<0.05:
                        g1_t = float(split_line[3])
                        g1_f = float(split_line[4])
                        g2_t = float(split_line[5])
                        g2_f = float(split_line[6])

                        if g1_t == 0:
                            fisher_results[split_line[2]].append( split_line[0] )

                        elif g2_t == 0:
                            fisher_results[split_line[1]].append( split_line[0] )

                        elif ( ((g1_t/(g1_t+g1_f))) / ((g2_t/(g2_t+g2_f))) )>1:
                            fisher_results[split_line[1]].append( split_line[0] )

                        else:
                            fisher_results[split_line[2]].append( split_line[0] )

        if len(fisher_results.keys())>0:
            return fisher_results

        else:
            raise Exception("Malformatted enrichment output: %s" % enrichment_output)


class ParseAnnotate:

    def __init__(self, enrichm_annotate_output, processes):
        self.path = enrichm_annotate_output
        # Parse genome objects
        self.genome_pickle_file_path \
            = os.path.join(enrichm_annotate_output, Annotate.GENOME_OBJ)
        self.processes \
            = processes

        ko = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_KO)
        if os.path.isfile(ko):
            self.ko = ko
        else:
            self.ko = None
        ko_hmm = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_KO_HMM)
        if os.path.isfile(ko_hmm):
            self.ko_hmm = ko_hmm
        else:
            self.ko_hmm = None
        pfam = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_PFAM)
        if os.path.isfile(pfam):
            self.pfam = pfam
        else:
            self.pfam = None
        tigrfam = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_TIGRFAM)
        if os.path.isfile(tigrfam):
            self.tigrfam = tigrfam
        else:
            self.tigrfam = None
        cazy = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_CAZY)
        if os.path.isfile(cazy):
            self.cazy = cazy
        else:
            self.cazy = None
        ec = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_EC)
        if os.path.isfile(ec):
            self.ec = ec
        else:
            self.ec = None
        cluster = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_CLUSTER)
        if os.path.isfile(cluster):
            self.cluster = cluster
        else:
            self.cluster = None
        ortholog = os.path.join(enrichm_annotate_output, Annotate.OUTPUT_ORTHOLOG)
        if os.path.isfile(ortholog):
            self.ortholog = ortholog
        else:
            self.ortholog = None

    def parse_pickles(self, path, genome_list):
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
        paths = list()

        for pickled_genome in genome_list:
            pickled_genome_path = os.path.join(path, pickled_genome + '.pickle')
            if os.path.isfile(pickled_genome_path):
                paths.append(pickled_genome_path)

        self.pool = mp.Pool(processes=self.processes)
        output_genome_list = self.pool.map_async(parse_genomes, paths)
        output_genome_list.wait()
        genome_objects = output_genome_list.get()
        self.pool.close()

        return genome_objects

class ParseGenerate:

    def __init__(self, forester_model_directory):
        self.LABELS_DICT = "labels_dict.pickle"
        self.RF_MODEL = "rf_model.pickle"
        self.ATTRIBUTE_IMPORTANCES = "attribute_importances.tsv"
        self.forester_model_directory = forester_model_directory

        for content in os.listdir(forester_model_directory):
            content_path = os.path.join(forester_model_directory, content)

            if content == self.LABELS_DICT:
                self.labels = pickle.load(open(content_path, 'rb'))
            elif content == self.RF_MODEL:
                self.model = pickle.load(open(content_path, 'rb'))
            elif content == self.ATTRIBUTE_IMPORTANCES:
                self.attributes = list()
                content_path_io = open(content_path)
                _ = content_path_io.readline() # Junk

                for line in content_path_io:
                    attribute, _ = line.strip().split('\t')
                    self.attributes.append(attribute)

        if None in [self.labels, self.model, self.attributes]:
            raise Exception("Malformatted forester model directory: %s" \
                            % (forester_model_directory))
