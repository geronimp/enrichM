#!/usr/bin/env python3
# pylint: disable=line-too-long
'''
Check for blocks of enriched genes that co-occur within the genome.
'''
import logging
from itertools import product, chain, combinations
from fuzzywuzzy import fuzz, process
from collections import Counter
from enrichm.toolbox import reverse_dictionary_of_lists, cluster

class SyntenySearcher:

    def __init__(self):
        self.synteny_results_output_file = "synteny_results.tsv"
        self.synteny_results_header = ["Genome_group", "Core_gene_block", "Num_genes", "Group_genomes", "Hit_genomes", "Percent_genomes",  "Genomes"]


    def cluster_operons_by_common_elements(self, synteny_to_genome, operon_mismatch_cutoff, operon_match_score_cutoff):
        cluster_list = list()
        
        for reference_operon in synteny_to_genome:
            operon_cluster = set([reference_operon])
            reference_operon_set = set(reference_operon.split('~'))
            reference_operon_length = len(reference_operon_set)

            for query_operon in synteny_to_genome:
                  
                if query_operon not in reference_operon:
                    query_operon_set = set(query_operon.split('~'))
                    query_operon_length = len(query_operon_set)
                    overlap_length = len(reference_operon_set.intersection(query_operon_set))
                    operon_match_score = overlap_length / reference_operon_length
                    mismatch_count = query_operon_length - overlap_length
                    
                    
                    if operon_match_score >= operon_match_score_cutoff:

                        if mismatch_count <= operon_mismatch_cutoff:
                            operon_cluster.add(query_operon)
            
            # Add if not duplicate or subset of another cluster list
            if len(cluster_list) == 0:
                cluster_list.append(operon_cluster)
            else:
                new = True
                
                for previous_operon_cluster in cluster_list:
                    previous_cluster_length = len(previous_operon_cluster)
                    current_cluster_length = len(operon_cluster)
                    overlap_length = len(previous_operon_cluster.intersection(operon_cluster))

                    if(overlap_length==previous_cluster_length or overlap_length==current_cluster_length):
                        new = False
                
                if new:
                    cluster_list.append(operon_cluster)
        
        #for synteny_block_a, synteny_block_b in combinations(synteny_to_genome.keys(), 2):
        #    synteny_block_a_set = set(synteny_block_a.split('~'))
        #    synteny_block_b_set = set(synteny_block_b.split('~'))
        #    synteny_overlap = set(synteny_block_a_set).intersection(set(synteny_block_b_set))
        #    
        #    if len(synteny_overlap)>=2: 
        #        # Quick and dirtly implementation. If two synteny blocks share >=2 common elements,
        #        # cluster them together ALTHOUGH this isnt a very effective method to protect against
        #        # cluster creep, which is dangerous.
        #        found = False

        #        for cl in cluster_list:
        #            
        #            if(synteny_block_a in cl or synteny_block_b in cl):
        #                core = set.intersection(*[set(x.split('~')) for x in cl])

        #                if(len(synteny_block_a_set.intersection(core))>=2 and len(synteny_block_b_set.intersection(core))>=2):
        #                    cl.update([synteny_block_a, synteny_block_b])
        #                    grouped_clusters.update([synteny_block_a, synteny_block_b])
        #                    found = True

        #        if not found:
        #            cluster_list.append(set([synteny_block_a, synteny_block_b]))
        #            grouped_clusters.update([synteny_block_a, synteny_block_b])
        #
        ## Add in ungrouped operons
        #cluster_list += [set([cluster]) for cluster in synteny_to_genome if cluster not in grouped_clusters]
        return cluster_list


    def search_for_blocks(self, prevalence_results, gene_positions, metadata, synteny_range, min_subblock_size,
                          operon_mismatch_cutoff, operon_match_score_cutoff):
        """
        Identify groups of significantly enriched genes that are also co-located between groups of genomes
    
        :type prevalence_results: list
        :param prevalence_results:
    
        :type gene_positions: dict
        :param gene_positions:
    
        :type metadata: dict
        :param metadata:
    
        :type synteny_range: int
        :param synteny_range:
    
        :type min_subblock_size: int
        :param min_subblock_size:
    
        :raises:
    
        :rtype:
        """
        group_to_enriched_gene = dict()
        metadata_reverse_mapping = reverse_dictionary_of_lists(metadata)
        output_lines = [self.synteny_results_header]

        # TODO This for loop can be functionized
        for result in prevalence_results[1:]:
            gene = result[0]
            enriched_in = result[3]
            corrected_pvalue = float(result[-2])

            # Just ignore genes that arent enriched on one group
            # for now
            if enriched_in == 'NA':
                continue

            # Ignore genes that aren't significantly different after multiple
            # test correction
            if corrected_pvalue > 0.05:
                continue
            
            if enriched_in not in group_to_enriched_gene:
                group_to_enriched_gene[enriched_in] = list()

            group_to_enriched_gene[enriched_in].append(gene)
        
        for group, enriched_gene_list in group_to_enriched_gene.items():
            synteny_to_genomes = dict()

            for genome in metadata_reverse_mapping[group]:
                
                for contig, genome_positions in gene_positions[genome].items(): 
                    # This is iterating through a genome on a contig-by-contig basis
                    # with nothing in place to allow operons to be identified across cotnigs
                    # This is stringent, but maybe necessary?
                    
                    endings_dict = dict()

                    for gene in enriched_gene_list:

                        if gene in genome_positions:
                            for pos in genome_positions[gene]:
                                endings_dict[pos[1]] = gene 
                        else:
                            logging.debug(f"Enriched gene {gene} not found on contig {contig} in {genome}. This situation is ignored.")

                    if len(endings_dict):

                        for subblock in cluster(list(endings_dict.keys()), synteny_range):

                            if len(subblock)>=min_subblock_size:
                                block = [endings_dict[bit] for bit in subblock]
                                synteny_definition = '~'.join(block)

                                if synteny_definition not in synteny_to_genomes:
                                    synteny_to_genomes[synteny_definition] = list()

                                synteny_to_genomes[synteny_definition].append(genome)
                    else:
                        logging.debug(f"No syntenous blocks found for {genome} on contig {contig}")
            
            clustered_operons = self.cluster_operons_by_common_elements(synteny_to_genomes, operon_mismatch_cutoff, operon_match_score_cutoff)
            
            for operon_group in clustered_operons:
                core_genes = '~'.join(set.intersection(*[set(x.split('~')) for x in operon_group]))
                found_genomes = set(chain(*[synteny_to_genomes[x] for x in operon_group]))
                perc_genomes = round(len(found_genomes) / len(metadata_reverse_mapping[group])*100, 2)
                output_lines.append([group, core_genes, len(core_genes.split('~')), len(metadata_reverse_mapping[group]), len(found_genomes), perc_genomes, "~".join(found_genomes)])

        return output_lines, self.synteny_results_output_file