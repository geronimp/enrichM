#!/usr/bin/env python3
from dataclasses import dataclass

from enrichm.writer import Writer
from enrichm.parser import Parser
from enrichm.databases import Databases
from enrichm.enrichment import mannwhitneyu_calc
from itertools import combinations

import logging
import os

@dataclass
class Uses:

    positive = 'positive'
    negative = 'negative'
    abundace = "frequency_matrix.tsv"
    enrichment = "enrichment_results.tsv"
    abundace_header = ["Compound"]
    enrichment_header = ["Compound", "Group_1", "Group_2", "group_1_mean", "group_2_mean", "score", "pvalue", "description"]
    databases: Databases()

    def gather_present_annotations(self, column_annotations):

        present_annotations = set()

        for annotation, copy_number in column_annotations.items():

            if copy_number>0:
                present_annotations.add(annotation)

        return present_annotations

    def uses(self, compound_list, annotations, column_names, count):
        output_lines_abundance = [self.abundace_header + column_names]
        enrichment_tallys = dict()

        for compound in compound_list:

            if compound in self.databases.c2r():

                enrichment_tallys[compound] = dict()
                abundance_line = [compound + '~' + self.databases.c()[compound]]

                for column_header in column_names:
                    column_positive_tally = 0
                    column_negative_tally = 0

                    # Gather all annotations present for this column (genome)
                    present_annotations = self.gather_present_annotations(annotations[column_header])
                    for reaction in self.databases.c2r()[compound]:

                        # If there are more than 0 KOs that carry out the reaction present in
                        # the genome
                        if reaction in self.databases.r2k():
                            overlapping_annotations = present_annotations.intersection(self.databases.r2k()[reaction])

                            if len(overlapping_annotations)>0:

                                if count:

                                    for annotation in overlapping_annotations:
                                        column_positive_tally+=annotations[column_header][annotation]

                                else:
                                    column_positive_tally+=1

                            else:
                                column_negative_tally+=1

                    abundance_line.append(column_positive_tally)
                    enrichment_tallys[compound][column_header] \
                        = {self.positive:column_positive_tally,
                           self.negative:column_negative_tally}

            output_lines_abundance.append(abundance_line)

        return output_lines_abundance, enrichment_tallys

    def enrichment(self, enrichment_tallys, metadata):
        output_lines = [self.enrichment_header]

        for compound, tallys in enrichment_tallys.items():

            for group_1_name, group_2_name in combinations(metadata.keys(), 2):
                group_1_tallys = {genome:tallys[genome] for genome in metadata[group_1_name]}
                group_2_tallys = {genome:tallys[genome] for genome in metadata[group_2_name]}
                group_1_values =  [tally[self.positive] for tally in group_1_tallys.values()]
                group_2_values =  [tally[self.positive] for tally in group_2_tallys.values()]

                output_line = mannwhitneyu_calc((compound, group_1_name, group_2_name, [group_1_values, None], [group_2_values, None]))
                output_lines.append(output_line + [self.databases.c()[compound]])

        return output_lines

    def uses_pipeline(self, compounds_list_path, annotation_matrix_path, metadata_path, output,
                      count):
        logging.info('Parsing input compounds list')
        compound_list = Parser.parse_single_column_text_file(compounds_list_path)
        logging.info('Parsing input annotations')
        annotations_dict, column_names, annotations = Parser.parse_simple_matrix(annotation_matrix_path)
        logging.info('Parsing input metadata')
        metadata, metadata_value_lists, attribute_dict = Parser.parse_metadata_matrix(metadata_path)
        logging.info('Tallying genes that use specified compounds')
        output_lines_abundance, enrichment_tallys = self.uses(compound_list, annotations_dict, column_names, count)
        logging.info('Writing file: %s' % self.abundace)
        Writer.write(output_lines_abundance, os.path.join(output, self.abundace))
        logging.info('Calculating enrichment between groups for each compound')
        output_lines_enrichment = self.enrichment(enrichment_tallys, attribute_dict)
        logging.info('Writing file: %s' % self.enrichment)
        Writer.write(output_lines_enrichment, os.path.join(output, self.enrichment))
        logging.info('Finished the use pipeline')