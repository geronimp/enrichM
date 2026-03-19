#!/usr/bin/env python3
# pylint: disable=line-too-long
'''
Check for blocks of enriched genes that co-occur within the genome,
using intergenic distance boundaries and strand/order information from GFF.
'''
import logging
from collections import Counter
from itertools import chain
from enrichm.toolbox import reverse_dictionary_of_lists


class SyntenySearcher:

    synteny_results_output_file = "synteny_results.tsv"
    synteny_results_header = [
        "Genome_group", "Core_gene_block", "Num_genes",
        "Group_genomes", "Hit_genomes", "Percent_genomes",
        "Ordered_sequence", "Strand_pattern", "Genomes",
    ]

    me_results_output_file = "mobile_element_proximity.tsv"
    me_results_header = [
        "Annotation", "Group", "N_genomes",
        "N_proximal", "Pct_proximal",
        "Mean_distance_bp", "Min_distance_bp",
    ]

    # Known mobile element annotation IDs (KO and Pfam).
    # Used to flag genes as mobile elements in flag_mobile_elements().
    MOBILE_ELEMENT_ANNOTATIONS = {
        # IS-family transposases
        'K07483', 'K07484', 'K07485', 'K07486', 'K07487',
        'K07488', 'K07489', 'K07490', 'K07491', 'K07492',
        'K07493', 'K07494', 'K07495', 'K07496', 'K07497',
        'K07498', 'K07499', 'K07500',
        # Site-specific recombinases (XerC, XerD)
        'K06400', 'K06401',
        # Phage integrase
        'K14059',
        # Pfam transposase domains
        'PF01609', 'PF00872', 'PF13586', 'PF13610', 'PF13751',
    }

    # -------------------------------------------------------------------------
    # Helpers
    # -------------------------------------------------------------------------

    @staticmethod
    def _split_by_intergenic_distance(contig_genes, threshold):
        """Split a sorted gene list into operon-candidate blocks.

        A new block starts wherever the gap between consecutive genes exceeds
        *threshold* bp (next_gene.start - prev_gene.end > threshold).
        Overlapping genes (negative gap) are always co-block.

        Parameters
        ----------
        contig_genes : list of (start, end, strand, [annotations])
            Must be sorted by start position.
        threshold : int
            Maximum intergenic distance (bp) within an operon.

        Returns
        -------
        list of lists, each inner list a run of co-operon genes.
        """
        if not contig_genes:
            return []
        blocks = []
        current = [contig_genes[0]]
        for gene in contig_genes[1:]:
            intergenic = gene[0] - current[-1][1]
            if intergenic <= threshold:
                current.append(gene)
            else:
                blocks.append(current)
                current = [gene]
        blocks.append(current)
        return blocks

    @staticmethod
    def _parse_enriched_genes(prevalence_results):
        """Return {group: set(enriched_annotation)} from a prevalence result table."""
        group_to_enriched = {}
        for result in prevalence_results[1:]:
            gene = result[0]
            enriched_in = result[3]
            corrected_pvalue = float(result[-2])
            if enriched_in == 'NA' or corrected_pvalue > 0.05:
                continue
            if enriched_in not in group_to_enriched:
                group_to_enriched[enriched_in] = set()
            group_to_enriched[enriched_in].add(gene)
        return group_to_enriched

    def cluster_operons_by_common_elements(self, synteny_to_genome,
                                           operon_mismatch_cutoff,
                                           operon_match_score_cutoff):
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

            if len(cluster_list) == 0:
                cluster_list.append(operon_cluster)
            else:
                new = True
                for previous_operon_cluster in cluster_list:
                    previous_cluster_length = len(previous_operon_cluster)
                    current_cluster_length = len(operon_cluster)
                    overlap_length = len(previous_operon_cluster.intersection(operon_cluster))
                    if (overlap_length == previous_cluster_length or
                            overlap_length == current_cluster_length):
                        new = False
                if new:
                    cluster_list.append(operon_cluster)

        return cluster_list

    # -------------------------------------------------------------------------
    # Synteny block detection
    # -------------------------------------------------------------------------

    def search_for_blocks(self, prevalence_results, gene_positions, gene_order,
                          metadata, intergenic_distance, min_subblock_size,
                          operon_mismatch_cutoff, operon_match_score_cutoff):
        """Identify enriched gene clusters sharing operon context between genome groups.

        Uses intergenic distance to define operon boundaries, and records gene
        order and strand pattern within each block.

        Parameters
        ----------
        prevalence_results : list
            Enrichment result table (Fisher/MWU output lines including header).
        gene_positions : dict
            feature_dict from parse_gff — kept for API compatibility, not used
            in the block-building loop (gene_order is used instead).
        gene_order : dict
            {genome: {contig: [(start, end, strand, [annotations])]}} sorted
            by start position, as returned by Parser.parse_gff.
        metadata : dict
            {genome: {group}} mapping from parse_metadata_matrix.
        intergenic_distance : int
            Maximum bp gap between adjacent genes within an operon.
        min_subblock_size : int
            Minimum number of enriched genes required in a block to report it.
        operon_mismatch_cutoff : int
            Allowed annotation mismatches when clustering similar blocks.
        operon_match_score_cutoff : float
            Minimum overlap fraction for clustering.
        """
        group_to_enriched_gene = self._parse_enriched_genes(prevalence_results)
        metadata_reverse_mapping = reverse_dictionary_of_lists(metadata)
        output_lines = [self.synteny_results_header]

        for group, enriched_genes in group_to_enriched_gene.items():
            # canonical_key (sorted ~-joined annotations) -> [genomes]
            synteny_to_genomes = {}
            # canonical_key -> Counter of ordered sequences across genomes
            block_sequences = {}
            # canonical_key -> Counter of strand patterns across genomes
            block_strand_patterns = {}

            for genome in metadata_reverse_mapping[group]:
                if genome not in gene_order:
                    logging.debug('Genome %s not found in gene_order — skipping', genome)
                    continue

                for contig, contig_genes in gene_order[genome].items():
                    operons = self._split_by_intergenic_distance(
                        contig_genes, intergenic_distance)

                    for operon in operons:
                        # Enriched genes present in this operon, in position order
                        enriched_in_operon = [
                            gene for gene in operon
                            if any(ann in enriched_genes for ann in gene[3])
                        ]
                        if len(enriched_in_operon) < min_subblock_size:
                            continue

                        # Canonical key: sorted annotation set for clustering
                        ann_set = sorted({
                            ann
                            for gene in enriched_in_operon
                            for ann in gene[3]
                            if ann in enriched_genes
                        })
                        canonical_key = '~'.join(ann_set)

                        # Ordered sequence: annotations in genome position order
                        ordered_seq = '~'.join(
                            ann
                            for gene in enriched_in_operon
                            for ann in gene[3]
                            if ann in enriched_genes
                        )

                        # Strand pattern: one character per enriched gene
                        strand_pat = '~'.join(gene[2] for gene in enriched_in_operon)

                        if canonical_key not in synteny_to_genomes:
                            synteny_to_genomes[canonical_key] = []
                            block_sequences[canonical_key] = Counter()
                            block_strand_patterns[canonical_key] = Counter()

                        synteny_to_genomes[canonical_key].append(genome)
                        block_sequences[canonical_key][ordered_seq] += 1
                        block_strand_patterns[canonical_key][strand_pat] += 1

            clustered_operons = self.cluster_operons_by_common_elements(
                synteny_to_genomes, operon_mismatch_cutoff, operon_match_score_cutoff)

            for operon_group in clustered_operons:
                core_genes_set = set.intersection(
                    *[set(x.split('~')) for x in operon_group])
                core_genes = '~'.join(sorted(core_genes_set))
                found_genomes = set(chain(
                    *[synteny_to_genomes[x] for x in operon_group]))
                perc_genomes = round(
                    len(found_genomes) / len(metadata_reverse_mapping[group]) * 100, 2)

                # Dominant ordered sequence and strand pattern across block variants
                merged_seqs = Counter()
                merged_strands = Counter()
                for key in operon_group:
                    merged_seqs.update(block_sequences.get(key, {}))
                    merged_strands.update(block_strand_patterns.get(key, {}))
                top_seq = merged_seqs.most_common(1)[0][0] if merged_seqs else ''
                top_strand = merged_strands.most_common(1)[0][0] if merged_strands else ''

                output_lines.append([
                    group, core_genes, len(core_genes_set),
                    len(metadata_reverse_mapping[group]),
                    len(found_genomes), perc_genomes,
                    top_seq, top_strand,
                    '~'.join(found_genomes),
                ])

        return output_lines, self.synteny_results_output_file

    # -------------------------------------------------------------------------
    # Mobile element proximity flagging
    # -------------------------------------------------------------------------

    def flag_mobile_elements(self, prevalence_results, gene_order, metadata,
                             me_distance_threshold):
        """Flag enriched annotations by proximity to mobile genetic elements.

        For each significantly enriched annotation in each group, finds the
        minimum distance (bp, by gene midpoint) to the nearest mobile element
        gene on the same contig and aggregates per (annotation, group).

        Parameters
        ----------
        prevalence_results : list
            Enrichment result table (same format as search_for_blocks input).
        gene_order : dict
            {genome: {contig: [(start, end, strand, [annotations])]}} from
            Parser.parse_gff.
        metadata : dict
            {genome: {group}} from parse_metadata_matrix.
        me_distance_threshold : int
            Distance (bp) within which a gene is considered ME-proximal.

        Returns
        -------
        (output_lines, filename)
        """
        group_to_enriched_gene = self._parse_enriched_genes(prevalence_results)
        metadata_reverse_mapping = reverse_dictionary_of_lists(metadata)

        # (annotation, group) -> {'n_genomes': int, 'distances': [float]}
        results_acc = {}

        for group, enriched_genes in group_to_enriched_gene.items():
            for genome in metadata_reverse_mapping[group]:
                if genome not in gene_order:
                    continue

                # Collect ME gene midpoints per contig for this genome
                me_midpoints = {}   # {contig: [midpoint, ...]}
                for contig, contig_genes in gene_order[genome].items():
                    for gene in contig_genes:
                        start, end, strand, annotations = gene
                        if any(ann in self.MOBILE_ELEMENT_ANNOTATIONS
                               for ann in annotations):
                            me_midpoints.setdefault(contig, []).append(
                                (start + end) / 2)

                # For each enriched annotation, find its minimum distance to an ME
                for annotation in enriched_genes:
                    min_dist = None
                    for contig, contig_genes in gene_order[genome].items():
                        contig_me = me_midpoints.get(contig)
                        if not contig_me:
                            continue
                        for gene in contig_genes:
                            start, end, strand, annotations = gene
                            if annotation not in annotations:
                                continue
                            gene_mid = (start + end) / 2
                            nearest = min(abs(gene_mid - me) for me in contig_me)
                            if min_dist is None or nearest < min_dist:
                                min_dist = nearest

                    if min_dist is None:
                        continue  # annotation not found on any contig with an ME

                    key = (annotation, group)
                    if key not in results_acc:
                        results_acc[key] = {'n_genomes': 0, 'distances': []}
                    results_acc[key]['n_genomes'] += 1
                    results_acc[key]['distances'].append(min_dist)

        output_lines = [self.me_results_header]
        for (annotation, group), data in sorted(results_acc.items()):
            n_genomes = data['n_genomes']
            distances = data['distances']
            n_proximal = sum(1 for d in distances if d <= me_distance_threshold)
            pct_proximal = round(n_proximal / n_genomes * 100, 2) if n_genomes else 0.0
            mean_dist = round(sum(distances) / len(distances), 2) if distances else 'NA'
            min_dist_val = round(min(distances), 2) if distances else 'NA'
            output_lines.append([
                annotation, group, n_genomes,
                n_proximal, pct_proximal,
                mean_dist, min_dist_val,
            ])

        return output_lines, self.me_results_output_file
