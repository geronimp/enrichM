#!/usr/bin/env python3
import os
import pickle
import json
import logging
import multiprocessing as mp
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

        with open(text_file) as text_file_io:
            for entry in text_file_io:
                entries.add(entry.strip())

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
    def parse_gff(gff_file):
        feature_dict = dict()
        genome_to_annotations_dict = dict()
        
        gff_file_io = open(gff_file)
        
        for line in gff_file_io:

            if line.startswith('#'):
                continue
            contig, _, _, start_pos, finish_pos, _, strand, _, attributes = line.strip().split('\t')
            # Terminal ';' sometimes found in poorly formatted gffs. They gotta go.
            if attributes.endswith(';'):
                attributes = attributes[:-1]
            attributes_dict = dict()

            for attribute in attributes.split(';'):
                attribute_key, attribute_value = attribute.split('=')
                attribute_value_list = attribute_value.split(',')

                if attribute_key in attributes_dict:
                    raise Exception(f"Key duplicate in GFF file: {attribute_key}")
                else:
                    attributes_dict[attribute_key] = attribute_value_list 

            genome = attributes_dict['genome'][0]

            if genome not in genome_to_annotations_dict:
                genome_to_annotations_dict[genome] = dict()
                feature_dict[genome] = dict()

            for attribute in attributes_dict['annotations']:

                if contig not in feature_dict[genome]:
                    feature_dict[genome][contig] = dict()
                if attribute in feature_dict[genome][contig]:
                    feature_dict[genome][contig][attribute].append([int(start_pos), int(finish_pos), strand])
                    genome_to_annotations_dict[genome][attribute] += 1
                else:
                    feature_dict[genome][contig][attribute] = [[int(start_pos), int(finish_pos), strand]]
                    genome_to_annotations_dict[genome][attribute] = 1

        gff_file_io.close()

        return feature_dict, genome_to_annotations_dict

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
        self.labels_dict = "labels_dict.pickle"
        self.rf_model = "rf_model.pickle"
        self.attribute_importances = "attribute_importances.tsv"
        self.forester_model_directory = forester_model_directory

        for content in os.listdir(forester_model_directory):
            content_path = os.path.join(forester_model_directory, content)

            if content == self.labels_dict:
                with open(content_path, 'rb') as content_path_io:
                    self.labels = pickle.load(content_path_io)
            elif content == self.rf_model:
                with open(content_path, 'rb') as content_path_io:
                    self.model = pickle.load(content_path_io)
            elif content == self.attribute_importances:
                self.attributes = list()
                with open(content_path) as content_path_io:
                    _ = content_path_io.readline() # Junk

                    for line in content_path_io:
                        attribute, _ = line.strip().split('\t')
                        self.attributes.append(attribute)

        if None in [self.labels, self.model, self.attributes]:
            raise Exception("Malformatted forester model directory: %s" \
                            % (forester_model_directory))

class RulesJson:
    _CURRENT = "0"
    VERSION = 'version'
    
    _GENERIC_KEYS = {
        "0": [
            "modules", # Modules with rules
            "version" # Version of database
        ]
    }
    _REQUIRED_KEYS = {"0": # DB version
                         # Type of rule
                        {"synteny": # Define synteny blocks
                            # values of rule
                            ["genes", # genes in synteny block
                             "range", # range of nucleic acids within which they are found
                             "breaks", # How many breaks are allowed in the synteny block. Ignored if strict.
                             "minsubblocksize", # How many genes must be found in split  synteny blocks. Ignored if strict.
                             "min"], # Minimum number of genes that must be found in synteny to pass.
                        #"homology": {},
                        "required": {}
                        }
                     }

    
    def load(self, rules_json_filepath):
        logging.debug(f"Loading rules JSON file: {rules_json_filepath}")

        rules_json_filepath_io = open(rules_json_filepath)
        self.contents_dict = json.load(rules_json_filepath_io)
        rules_json_filepath_io.close()

        logging.debug("Finding rules version")
        current_package_version = self.contents_dict[self.VERSION]

        logging.debug(f"Loading rules version {current_package_version}")
        if current_package_version == '0':
            rj = RulesJsonVersion0(self.contents_dict, rules_json_filepath, current_package_version)
        else:
            raise Exception(f"JSON rules file version invalid: {current_package_version}")

        logging.debug("Checking package format is valid")
        self.check_version_keys(self._GENERIC_KEYS[current_package_version], self.contents_dict)

        for module_rules in self.contents_dict['modules'].values():
            self.check_version_keys(self._REQUIRED_KEYS[current_package_version], module_rules)

        return rj

    def get_version(self):
        if self.VERSION in self.contents_dict:
            return self.contents_dict[self.VERSION]
        else:
            raise Exception("No version information in rules JSON file")

    def check_version_keys(self, version_keys, contents_dict):

        if isinstance(version_keys, dict):

            for key, item in version_keys.items():
                self.check_key(key, contents_dict, item)

        elif isinstance(version_keys, list):

            for key in version_keys:
                self.check_key(key, contents_dict)

    def check_key(self, key, contents, item=None):
        
        if isinstance(contents, dict):
            if key in contents:
                if item:
                    self.check_version_keys(item, contents[key])
            else:
                raise Exception(f"Malformatted rules file. Key not found: {key}")
        
        elif isinstance(contents, list):
            for content in contents:
                if key in content:
                    if item:
                        self.check_version_keys(item, contents[key])
                else:
                    raise Exception(f"Malformatted rules file. Key not found: {key}")

class RulesJsonVersion0(RulesJson):

    def __init__(self, _contents_json, _filename, _version):
        self.VERSION = 0
        self._contents_json = _contents_json
        self._filename = _filename
        self._version = _version
        self.members = set(self._contents_json['modules'].keys())

    def get_synteny_rules(self, module_name):
        return self._contents_json['modules'][module_name]['synteny']

    def get_homology_rules(self, module_name):
        return self._contents_json['modules'][module_name]['homology']

    def get_required_rules(self, module_name):
        return self._contents_json['modules'][module_name]['required']
