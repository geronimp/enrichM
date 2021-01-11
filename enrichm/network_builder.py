#!/usr/bin/env python3
# Imports
import logging
import re
from itertools import product
# Local
from enrichm.databases import Databases
from enrichm.parser import Parser
###############################################################################

class NetworkBuilder:

    MODULE_PREFIX   = 'M'
    MAP_PREFIX      = 'map'
    PATHWAY_PREFIX  = 'rn'
    REACTION_PREFIX = 'R'
    ZERO            = '0.0'

    def __init__(self, metadata, abundances_metagenome, abundances_transcriptome, abundances_metabolome, fisher_results):
        self.abundances_metagenome = abundances_metagenome
        self.abundances_transcriptome = abundances_transcriptome
        self.abundances_metabolome = abundances_metabolome
        self.fisher_results = fisher_results
        self.metadata_keys = list(metadata.keys())

        databases = Databases()

        self.reaction_to_module = databases.r2m()
        self.module_to_reaction = databases.m2r()
        self.module_descriptions = databases.m()
        self.reaction_to_pathway = databases.r2p()
        self.pathway_to_reaction = databases.p2r()
        self.pathway_descriptions = databases.p()
        self.compound_desc_dict = databases.compound_desc_dict()
        self.compound_descriptions = databases.c()
        self.reaction_descriptions = databases.r()
        self.reactions_to_compounds = databases.r2c()
        self.reactions_to_kos = databases.r2k()

        self.matrix_header = ["compound", "reaction", 'type']
        self.transcriptome_header = [key + '_reaction_transcriptome' for key in self.metadata_keys]
        self.compound_header = [key + '_compound' for key in self.metadata_keys]
        self.metadata_header = ['node', 'description', 'type', 'module', 'module_descr', 'pathway', 'pathway_descr', 'node_type']
        self.query_header =  ['query', 'step']
        self.step_header = ['step']
        self.to_omit = set(["C00828",  # Menaquinone
                            "C00534",  # Pyridoxamine
                            "C00006",  # NADP+
                            "C00003",  # NAD+
                            "C00002",  # ATP
                            "C00314",  # Pyridoxine
                            "C00864",  # Pantothenate
                            "C00504",  # Folate
                            "C00032",  # Heme
                            "C05443",  # Vitamin D3
                            "C00253",  # Nicotinate
                            "C00250",  # Pyridoxal
                            "C11378",  # Ubiquinone-10
                            "C05777",  # Coenzyme F430
                            "C00072",  # Ascorbate
                            "C00378",  # Thiamine
                            "C00101",  # Tetrahydrofolate
                            "C00029",  # UDP-glucose
                            "C00068",  # Thiamin diphosphate
                            "C00061",  # FMN
                            "C00063",  # CTP
                            "C05776",  # Vitamin B12
                            "C00113",  # PQQ
                            "C18237",  # Molybdoenzyme molybdenum cofactor
                            "C00051",  # Glutathione
                            "C00010",  # CoA
                            "C00016",  # FAD
                            "C00018",  # Pyridoxal phosphate
                            "C00019",  # S-Adenosyl-L-methionine
                            "C00153",  # Nicotinamide
                            "C04628",  # Coenzyme B
                            "C00862",  # Methanofuran
                            "C15672",  # Heme O
                            "C15670",  # Heme A
                            "C02059",  # Phylloquinone
                            "C03576",  # Coenzyme M
                            "C05441",  # Vitamin D2
                            "C00272",  # Tetrahydrobiopterin
                            "C02477",  # alpha-Tocopherol
                            "C00473",  # Retinol
                            "C00120",  # Biotin
                            "C00725",  # Lipoate
                            "C00053",  # 3'-Phosphoadenylyl sulfate
                            "C00194",  # Cobamide coenzyme
                            "C00255",  # Riboflavin
                            'C00001',  # H2O
                            'C00008',  # ADP
                            'C00013',  # Diphosphate
                            'C00004',  # NADH
                            'C00005',  # NADPH
                            'C00080',  # H+
                            'C00009',  # Orthophosphate
                            'C00008',  # ADP
                            'C00004',  # NADH
                            'C00020',  # AMP
                            'C00007',  # Oxygen
                            'C00015']) # UDP

    def nested_dict_vals(self, input_dictionary):
        reaction_regex = '(R\d{5})$'

        for key, item in input_dictionary.items():

            if isinstance(item, dict):
                yield from self.nested_dict_vals(item)
            else:

                if type(reaction_regex) == str:

                    if re.match(reaction_regex, key):
                        yield key

                elif type(reaction_regex) == bytes:
                    key = key.decode()

                    if re.match(reaction_regex, key):
                        yield key

                else:
                    raise Exception("Invalid key in in nested dict!")

    def gather_module(self, key):

        if key in self.reaction_to_module:
            len_list =  [len(self.module_to_reaction[x]) for x in self.reaction_to_module[key]]
            module = self.reaction_to_module[key][len_list.index(max(len_list))]
            module_description = self.module_descriptions[module]
        else:
            module = 'NA'
            module_description='NA'

        return module, module_description

    def gather_pathway(self, key):

        if key in self.reaction_to_pathway:
            len_list =  [len(self.pathway_to_reaction[x]) for x in self.reaction_to_pathway[key]]
            pathway = self.reaction_to_pathway[key][len_list.index(max(len_list))]
            pathway_description = self.pathway_descriptions[pathway]
        else:
            pathway = 'NA'
            pathway_description='NA'

        return pathway, pathway_description

    def gather_compound_metadata(self, compound):

        if compound in self.compound_desc_dict:
            compound_type = ','.join(self.compound_desc_dict[compound]['A'])
        else:
            compound_type = 'NA'

        compound_metadata_list = [compound, self.compound_descriptions[compound], compound_type, 'NA', 'NA', 'NA', 'NA', 'compound']

        if self.abundances_metabolome:

            for sample in self.metadata_keys:

                if compound in self.abundances_metabolome[sample]:
                    compound_abundance = self.abundances_metabolome[sample][compound]
                else:
                    compound_abundance = 0

                compound_metadata_list.append(compound_abundance)

        return compound_metadata_list

    def gather_reaction_metadata(self, reaction):
        module, module_description = self.gather_module(reaction)
        pathway, pathway_description = self.gather_pathway(reaction)
        reaction_metadata_list = [reaction, self.reaction_descriptions[reaction], 'NA', module, module_description, pathway, pathway_description, 'reaction']

        if self.abundances_metabolome:
            reaction_metadata_list += ['-5' for sample in self.metadata_keys]

        return reaction_metadata_list

    def gather_enriched_term(self, reaction, fisher_results):
        enriched_term = list()

        for compared_group in list(fisher_results.keys()):

            if any(set(self.reactions_to_kos[reaction]).intersection(fisher_results[compared_group])):
                enriched_term.append(compared_group)

        if len(enriched_term) > 0:
            enriched_term = '_'.join(enriched_term)
        else:
            enriched_term = 'NA'

        return enriched_term

    def get_reaction_compounds(self, compounds_list):

        for compound in compounds_list:

            if compound not in self.to_omit:
                yield compound

    def gather_reaction_edge_data(self, compound, reaction):

        if self.fisher_results:
            enriched_term = self.gather_enriched_term(reaction, self.fisher_results)
        else:
            enriched_term = 'NA'

        reaction_line = [compound, reaction, enriched_term]

        for key in self.metadata_keys:

            for _, group_abundances in self.abundances_metagenome.items():
                if reaction in group_abundances[key]:
                    reaction_line.append(str(group_abundances[key][reaction]))
                else:
                    reaction_line.append(self.ZERO)

        if self.abundances_transcriptome:

            for key in self.metadata_keys:
                for _, group_abundances in self.abundances_transcriptome.items():
                    if reaction in group_abundances[key]:
                        reaction_line.append(str(group_abundances[key][reaction]))
                    else:
                        reaction_line.append(self.ZERO)

        if sum([float(x) for x in reaction_line[3:]])>0:
            return reaction_line

    def add_to_header(self, groups):
        output_list = list()
        groups = list(groups.keys())

        for genome_group, reference_group in product(self.metadata_keys, groups):
            output_list.append(f"{genome_group}_{reference_group}")
        
        return output_list

    def all_matrix(self, reference_dict):
        '''
        Parameters
        ----------
        '''
        seen_nodes = set()
        
        # Construct headers to network matrices
        seen_reactions = set(self.nested_dict_vals(self.abundances_metagenome))

        network_lines = [self.matrix_header]

        if self.abundances_metagenome:
            network_lines[0] += self.add_to_header(self.abundances_metagenome)
        if self.abundances_transcriptome:
            network_lines[0] += self.add_to_header(self.abundances_transcriptome)
        if self.abundances_metabolome:
            node_metadata_lines = [self.metadata_header + self.compound_header]
        else:
            node_metadata_lines = [self.metadata_header]
        
        for reaction, entry in reference_dict.items():

            for compound in entry:

                if reaction in seen_reactions:

                    #~#~#~#~#~#~#~#~ Fill in abundances_metagenome ~#~#~#~#~#~#~#~
                    reaction_line = self.gather_reaction_edge_data(compound, reaction)

                    if reaction_line:
                        network_lines.append(reaction_line)

                    #~#~#~#~#~#~#~#~#~ Fill in node metadata ~#~#~#~#~#~#~#~#~#
                    if compound not in seen_nodes:
                        compound_metadata_list = self.gather_compound_metadata(compound)

                        node_metadata_lines.append(compound_metadata_list)
                        seen_nodes.add(compound)

                    if reaction not in seen_nodes:
                        reaction_metadata_list = self.gather_reaction_metadata(reaction)

                        node_metadata_lines.append(reaction_metadata_list)
                        seen_nodes.add(reaction)

        return network_lines, node_metadata_lines

    def query_matrix(self, queries, depth):
        '''
        Parameters
        ----------
        '''
        steps = 0
        queries_list = Parser.parse_single_column_text_file(queries)
        seen_reactions = set(self.nested_dict_vals(self.abundances_metagenome))
        seen_nodes = set()
        level_queries = set()

        node_metadata_lines = [self.metadata_header + self.query_header]

        if self.abundances_metabolome:
            node_metadata_lines[0] += self.compound_header

        network_lines = [self.matrix_header]
        if self.abundances_metagenome:
            network_lines[0] += self.add_to_header(self.abundances_metagenome)
        if self.abundances_transcriptome:
            network_lines[0] += self.add_to_header(self.abundances_transcriptome)
        network_lines[0] += self.step_header

        while depth>0:

            if any(level_queries):
                queries_list = set(level_queries)
                level_queries = set()

            for reaction, entry in self.reactions_to_compounds.items():

                if reaction in seen_reactions:

                    if any(queries_list.intersection(entry)):

                        for compound in self.get_reaction_compounds(entry):

                            #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                            #~#~#~#~#~#~#~#~ Fill in abundances ~#~#~#~#~#~#~#~
                            reaction_line = self.gather_reaction_edge_data(compound, reaction)

                            if reaction_line:
                                reaction_line.append(str(steps))
                                network_lines.append(reaction_line)

                            #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
                            #~#~#~#~#~#~#~ Fill in node metadata  ~#~#~#~#~#~#~
                            index = str((steps if compound in queries_list else steps+1))

                            if compound not in seen_nodes:

                                metadata_list = self.gather_compound_metadata(compound)

                                if compound in queries_list:
                                    is_query = 'True'
                                else:
                                    is_query = 'False'

                                node_metadata_lines.append(metadata_list + [is_query, index])
                                seen_nodes.add(compound)

                                if compound not in queries_list:
                                    level_queries.add(compound)

                            if reaction not in seen_nodes:
                                metadata_list = self.gather_reaction_metadata(reaction)

                                node_metadata_lines.append(metadata_list + ['False', index])
                                seen_nodes.add(reaction)

            steps+=1
            depth-=1

            logging.info("Step %i complete with %i queries to continue with" % (steps, len(level_queries)))

        return network_lines, node_metadata_lines

    def pathway_matrix(self, limit, filter):

        possible_reactions = set()

        if any(limit):

            for entry in limit:

                if(entry.startswith(self.MAP_PREFIX) or
                   entry.startswith(self.PATHWAY_PREFIX)):

                    for reaction in self.pathway_to_reaction[entry]:
                        possible_reactions.add(reaction)

                elif(entry.startswith(self.MODULE_PREFIX)):

                    for reaction in self.module_to_reaction[entry]:
                        possible_reactions.add(reaction)

                elif(entry.startswith(self.REACTION_PREFIX)):
                    possible_reactions.add(entry)

            possible_reactions = {reaction:self.reactions_to_compounds[reaction]
                                  for reaction in
                                  possible_reactions
                                  if reaction in self.reactions_to_compounds}
        else:
            possible_reactions = self.reactions_to_compounds

        for entry in filter:

            if entry in possible_reactions:

                del possible_reactions[entry]

        network_lines, node_metadata_lines = self.all_matrix(possible_reactions)
        return network_lines, node_metadata_lines

