#!/usr/bin/env python3
# pylint: disable=line-too-long
import math
import os
import re
import logging
import multiprocessing as mp
from itertools import product, combinations, chain
from scipy import stats
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import pdist
import numpy as np
import statsmodels.sandbox.stats.multicomp as sm
from sklearn.decomposition import NMF
import dendropy
from enrichm.databases import Databases
from enrichm.module_description_parser import ModuleDescription
from enrichm.parser import Parser, ParseAnnotate
from enrichm.writer import Writer
from enrichm.synteny_searcher import SyntenySearcher
################################################################################

def indval_calc(x):
    '''
    Calculate Indicator Value (IndVal) for a single annotation in a focal group.
    x = [annotation, focal_group, focal_abundances, all_group_abundances, n_permutations]
    Returns [annotation, focal_group, indval, specificity, fidelity, pvalue]
    '''
    annotation, focal_group, focal_abundances, all_group_abundances, n_perms = x
    focal_arr = np.array(focal_abundances, dtype=float)
    focal_mean = np.mean(focal_arr)
    all_means = np.array([np.mean(np.array(v, dtype=float)) for v in all_group_abundances.values()])
    total_mean = np.sum(all_means)

    if total_mean == 0:
        return None

    specificity = focal_mean / total_mean
    fidelity = np.sum(focal_arr > 0) / len(focal_arr)
    indval = np.sqrt(specificity * fidelity)

    # Permutation p-value: pool all values, shuffle, recompute IndVal for focal group
    all_vals = np.array(list(chain(*[list(v) for v in all_group_abundances.values()])), dtype=float)
    focal_size = len(focal_arr)
    group_sizes = [len(v) for v in all_group_abundances.values()]
    focal_idx = list(all_group_abundances.keys()).index(focal_group)

    rng = np.random.default_rng()
    perm_indvals = np.empty(n_perms)

    for i in range(n_perms):
        shuffled = rng.permutation(all_vals)
        idx = 0
        perm_means = []
        for j, size in enumerate(group_sizes):
            chunk = shuffled[idx:idx + size]
            perm_means.append(np.mean(chunk))
            idx += size
        perm_total = np.sum(perm_means)
        if perm_total == 0:
            perm_indvals[i] = 0.0
            continue
        perm_focal = shuffled[sum(group_sizes[:focal_idx]):sum(group_sizes[:focal_idx]) + focal_size]
        perm_spec = perm_means[focal_idx] / perm_total
        perm_fid = np.sum(perm_focal > 0) / focal_size
        perm_indvals[i] = np.sqrt(perm_spec * perm_fid)

    pvalue = (np.sum(perm_indvals >= indval) + 1) / (n_perms + 1)
    return [annotation, focal_group, str(round(float(indval), 4)),
            str(round(float(specificity), 4)), str(round(float(fidelity), 4)), pvalue]


def gene_fisher_calc(x):
    annotation, group_1, group_2 = x[0], x[1], x[2]

    dat = x[3:]

    if (dat[0][0]>0 or dat[1][0]>0) and (dat[0][1]>0 or dat[1][1]>0):
        score, pval = stats.fisher_exact(dat)
        if (dat[0][0] / sum(dat[0])) > (dat[1][0] / sum(dat[1])):
            enriched_in = group_1
        else:
            enriched_in = group_2
    else:
        enriched_in = 'NA'
        score, pval = 'nan', 1.0

    return [annotation, group_1, group_2, enriched_in] + dat[0] + dat[1] + [score, pval]

def mannwhitneyu_calc(x):
    # Mann Whitney U test
    annotation, group_1, group_2, group_1_module_annotations, group_2_module_annotations = x
    group_1_module_annotations = np.array(group_1_module_annotations[0])
    group_2_module_annotations = np.array(group_2_module_annotations[0])
    group_1_mean = np.mean(group_1_module_annotations)
    group_2_mean = np.mean(group_2_module_annotations)

    if(sum(group_1_module_annotations)>0 and sum(group_2_module_annotations)>0):
        if(len(set(group_1_module_annotations)) == 1
            or
           len(set(group_2_module_annotations)) == 1):
            mw_t_stat, mw_p_value = 'NA', 1
            enriched_in = 'NA'

        else:

            if group_1_mean > group_2_mean:
                enriched_in = group_1
            else:
                enriched_in = group_2

            mw_t_stat, mw_p_value = \
                stats.mannwhitneyu(group_1_module_annotations,
                                         group_2_module_annotations)
    else:
        
        mw_t_stat, mw_p_value = 'NA', 1
        enriched_in = 'NA'

    n1 = len(group_1_module_annotations)
    n2 = len(group_2_module_annotations)

    if mw_t_stat == 'NA':
        fold_change = 'NA'
        effect_size = 'NA'
    else:
        effect_size = (2 * mw_t_stat - n1 * n2) / (n1 * n2)
        if group_2_mean == 0:
            fold_change = 'inf' if group_1_mean > 0 else 'NA'
        else:
            fold_change = group_1_mean / group_2_mean

    return [annotation, group_1, group_2, enriched_in,
            str(group_1_mean), str(group_2_mean),
            str(fold_change), str(effect_size),
            mw_t_stat, mw_p_value]

def zscore_calc(x):

    annotation, group_1, group_2, group_1_module_annotations, group_2_module_annotations = x
    group_1_module_annotations = group_1_module_annotations[0]
    group_2_module_annotations = group_2_module_annotations[0]

    if len(group_1_module_annotations)>1:
        reference = group_1_module_annotations
        reference_name = group_1
        genome = group_2_module_annotations[0]
        genome_name = group_2

    else:
        reference_name = group_2
        reference = group_2_module_annotations
        genome = group_1_module_annotations[0]
        genome_name = group_1

    if genome>0:
        reference_group_comp_sd = np.std(reference, axis=0)
        reference_group_comp_mean = np.mean(reference, axis=0)

        if (genome-reference_group_comp_mean)>0:

            if reference_group_comp_sd==0:
                z_score = np.inf
                p_value = 0.0
                enriched_in = 'NA'

            else:
                z_score = (genome-reference_group_comp_mean) / reference_group_comp_sd
                p_value = 2-2*stats.norm.cdf(z_score)
                enriched_in = genome_name

            return [annotation,
                    reference_name,
                    genome_name,
                    enriched_in,
                    str(reference_group_comp_mean),
                    str(reference_group_comp_sd),
                    str(genome),
                    str(z_score),
                    p_value]

def kruskal_wallis_calc(x):
    annotation, group_names, group_values_list = x
    group_arrays = [np.array(v) for v in group_values_list]
    non_empty = [(n, v) for n, v in zip(group_names, group_arrays) if sum(v) > 0]
    if len(non_empty) < 2:
        return None
    try:
        h_stat, pvalue = stats.kruskal(*[v for _, v in non_empty])
    except ValueError:
        return None
    return [annotation, ','.join(group_names), h_stat, pvalue]

def phylo_pairs_calc(x):
    '''
    Scoary1-style phylogenetic paired comparison for one annotation.
    x = [annotation, pairs, genome_annotations, n_permutations]
    pairs = list of (g1_genome, g2_genome) tuples derived from tree splits.
    Returns [annotation, concordant, discordant, pvalue] or None.
    '''
    annotation, pairs, genome_annotations, n_perms = x
    C, D = 0, 0
    for g1_g, g2_g in pairs:
        g1_has = genome_annotations.get(g1_g, {}).get(annotation, 0) > 0
        g2_has = genome_annotations.get(g2_g, {}).get(annotation, 0) > 0
        if g1_has and not g2_has:
            C += 1
        elif not g1_has and g2_has:
            D += 1
    if C + D == 0:
        return None
    # Within-pair permutation: randomly swap which genome is treated as g1
    rng = np.random.default_rng()
    perm_C = np.zeros(n_perms, dtype=int)
    for i in range(n_perms):
        c = 0
        for g1_g, g2_g in pairs:
            if rng.random() < 0.5:
                g1_g, g2_g = g2_g, g1_g
            g1_has = genome_annotations.get(g1_g, {}).get(annotation, 0) > 0
            g2_has = genome_annotations.get(g2_g, {}).get(annotation, 0) > 0
            if g1_has and not g2_has:
                c += 1
        perm_C[i] = c
    pvalue = (np.sum(perm_C >= C) + 1) / (n_perms + 1)
    return [annotation, C, D, pvalue]

################################################################################

class Enrichment:

    TIGRFAM = "tigrfam"
    PFAM = "pfam"
    KEGG = "kegg"
    CAZY = "cazy"
    EC = "ec"
    CLUSTER = "cluster"
    ORTHOLOG = "ortholog"
    OTHER = "other"
    COG = "cog"
    GO = "go"
    EGGNOG = "eggnog"
    TIGRFAM_PREFIX = 'TIGR'
    PFAM_PREFIX = 'PF'
    KEGG_PREFIX = 'K'
    GO_PREFIX = 'GO:'
    CAZY_PREFIX = ["GH", "AA", "GT", "PL", "CE", "CBM", "SLH", "dockerin", "cohesin", "GTCellulosesynt"]
    EC_PREFIX = ["1", "2", "3","4","5","6", "7"]
    KO_PATTERN = re.compile(r'^K\d{5}$')
    PROPORTIONS = 'proportions.tsv'
    MODULE_COMPLETENESS = 'modules.tsv'

    def _annotation_type_of(self, annotation):
        cazy_prefix = ''.join(c for c in annotation if not c.isdigit() and c != '_')
        if annotation.startswith(self.TIGRFAM_PREFIX):
            return self.TIGRFAM
        elif self.KO_PATTERN.match(annotation):
            return self.KEGG
        elif annotation.startswith(self.PFAM_PREFIX):
            return self.PFAM
        elif cazy_prefix in self.CAZY_PREFIX:
            return self.CAZY
        elif annotation.split('.')[0] in self.EC_PREFIX:
            return self.EC
        elif annotation.startswith(self.GO_PREFIX):
            return self.GO
        elif len(annotation) <= 2 and annotation.isalpha() and annotation.isupper():
            return self.COG
        elif '@' in annotation:
            return self.EGGNOG
        else:
            return self.OTHER

    def check_annotation_type(self, annotations):
        '''
        Inspects all rownames from the input matrix and determines the annotation
        type. Raises ValueError if mixed annotation types are detected.

        Parameters
        ----------
        annotations     - List. A list of strings, each a rowname
                          from the original input matrix

        Output
        ------
        The annotation type
        '''
        types = {self._annotation_type_of(a) for a in annotations}
        if len(types) > 1:
            raise ValueError(
                f"Mixed annotation types detected: {types}. "
                "All annotations must be the same type."
            )
        return types.pop()

    def weight_annotation_matrix(self,
                                 sample_abundance,
                                 annotation_abundance,
                                 sample_dict,
                                 annotations):

        output_dict = {sample_group: dict() for sample_group in sample_dict.keys()}
        
        logging.info('Aggregating abundances across samples')
        for group, samples in sample_dict.items():
            output_dict[group] = dict()

            for annotation in annotations:
                output_dict[group][annotation] = list()

                for sample in samples:
                    sample_annotation_abundance = 0.0

                    for genome, genome_annotation_dict in annotation_abundance.items():

                        if annotation in genome_annotation_dict:
                            value = genome_annotation_dict[annotation]
                            if genome in sample_abundance[sample]:
                                sample_annotation_abundance += sample_abundance[sample][genome]*value

                    output_dict[group][annotation].append(sample_annotation_abundance)

        return output_dict

    def calculate_portions(self,
                           annotations,
                           combination_dict,
                           annotations_dict,
                           genome_list,
                           proportions_cutoff):
        '''
        Calculates the portions of genome

        Parameters
        ----------
        modules             - List. List of all possible annotations for the given annotation
                              type (eg, all ko ids, or all pfam ids).
        combination_dict    - Dictionary. Metadata dictionary, with the groups as keys, and lists
                              of genomes as entries
        annotations_dict    - Dictionary. Annotation dictionary, with the the genome ids as keys,
                              and a list of annotations as the entry for each.
        genome_list         - List. List of strings, each one a genome name
        proportions_cutoff  - Float. Value with which to cutoff
        '''
        raw_proportions_output_lines        = [['Annotation'] + list(combination_dict.keys())]

        for annotation in annotations:

            annotation_values               = dict()
            raw_proportions_output_line = [annotation]

            for group_name, genome_list in combination_dict.items():
                if len(genome_list)>0:
                    coverage = len([genome for genome in genome_list
                                    if annotation in annotations_dict[genome]])
                    total    = float(len(genome_list))
                    entry    = coverage/total
                    annotation_values[group_name] = entry
                    raw_proportions_output_line.append(str(entry))
                else:
                    raw_proportions_output_line.append('0.0')

            if annotation_values and any(v >= proportions_cutoff for v in annotation_values.values()):
                raw_proportions_output_lines.append(raw_proportions_output_line)

        return raw_proportions_output_lines

    @staticmethod
    def filter_by_prevalence(annotations_dict, min_prevalence):
        """Remove annotations present in fewer than min_prevalence fraction of genomes."""
        if min_prevalence <= 0:
            return annotations_dict
        n_genomes = len(annotations_dict)
        if n_genomes == 0:
            return annotations_dict
        counts = {}
        for genome_anns in annotations_dict.values():
            for ann, val in genome_anns.items():
                if val > 0:
                    counts[ann] = counts.get(ann, 0) + 1
        min_count = min_prevalence * n_genomes
        keep = {ann for ann, c in counts.items() if c >= min_count}
        logging.info(
            f'Prevalence filter: keeping {len(keep)} of {len(counts)} annotations '
            f'(>= {min_prevalence * 100:.1f}% of {n_genomes} genomes)'
        )
        return {
            genome: {ann: val for ann, val in anns.items() if ann in keep}
            for genome, anns in annotations_dict.items()
        }

    def module_completeness(self, database, result_file_path, pval_cutoff):
        module_output = [["Module", "Lineage", "Total steps", "Steps covered", "Percentage covered", "Module description"]]
        module_descriptions = database.m()

        g1_sig_annotations = set()
        g2_sig_annotations = set()

        result_file_io = open(result_file_path)
        result_file_io.readline()

        for line in result_file_io:
            sline = line.strip().split('\t')
            if float(sline[-2])<pval_cutoff:
                if result_file_path.endswith("fisher.tsv"):
                    g1 = float(sline[4]) / (int(sline[4]) + int(sline[5]))
                    g2 = float(sline[6]) / (int(sline[6]) + int(sline[7]))
                elif result_file_path.endswith("cdf.tsv"):
                    g1 = float(sline[4])
                    g2 = float(sline[6])
                if g1>g2:
                    g1_sig_annotations.add(sline[0])
                else:
                    g2_sig_annotations.add(sline[0])

        for module, definition in database.m2def().items():

            if module not in database.signature_modules:
                pathway = ModuleDescription(definition)
                num_all = pathway.num_steps()
                g1_num_covered, _, _, _ = pathway.num_covered_steps(g1_sig_annotations)
                g1_perc_covered = g1_num_covered / float(num_all)

                g2_num_covered, _, _, _ = pathway.num_covered_steps(g2_sig_annotations)
                g2_perc_covered = g2_num_covered / float(num_all)

                if g1_perc_covered>0:
                    output_line = [module, sline[1], num_all, g1_num_covered, g1_perc_covered, module_descriptions[module]]
                    module_output.append(output_line)

                if g2_perc_covered>0:
                    output_line = [module, sline[2], num_all, g2_num_covered, g2_perc_covered, module_descriptions[module]]
                    module_output.append(output_line)

        prefix = '_vs_'.join([sline[1], sline[2]]).replace(' ', '_')

        return module_output, prefix

    def enrichment_pipeline(# Input options
           self, annotate_output, annotation_matrix, gff_files, dram_output, emapper_output,
           metadata_path, abundances_path, abundance_metadata_path, pval_cutoff,
           proportions_cutoff, min_prevalence, threshold, multi_test_correction, processes,
           ko, pfam, tigrfam, cluster, ortholog, cazy, ec, ko_hmm, cog, go, eggnog,
           intergenic_distance, subblock_size, operon_mismatch_cutoff, operon_match_score_cutoff,
           me_distance, output_directory,
           decompose=False, n_components=None, select_components=False,
           tree_path=None):
        
        database = Databases()
        syntenysearcher = SyntenySearcher()

        if gff_files:
            logging.info("Parsing .gff file input(s)")
            annotations_dict = dict()
            gene_positions = dict()
            gene_order = dict()

            for gff_file in gff_files:
                position, counts, order = Parser.parse_gff(gff_file)
                annotations_dict.update(counts)
                gene_positions.update(position)
                gene_order.update(order)

            annotations = set(chain(*[list(x.keys()) for x in annotations_dict.values()]))

        elif annotate_output or annotation_matrix:
            if annotate_output:
                logging.info('Parsing annotate output: %s' % (annotate_output))
                pa = ParseAnnotate(annotate_output, processes)

                if ko:
                    annotation_matrix = pa.ko
                elif ko_hmm:
                    annotation_matrix = pa.ko_hmm
                elif pfam:
                    annotation_matrix = pa.pfam
                elif tigrfam:
                    annotation_matrix = pa.tigrfam
                elif cluster:
                    annotation_matrix = pa.cluster
                elif ortholog:
                    annotation_matrix = pa.ortholog
                elif cazy:
                    annotation_matrix = pa.cazy
                elif ec:
                    annotation_matrix = pa.ec

            logging.info('Parsing annotation matrix')
            annotations_dict, _, annotations = Parser.parse_simple_matrix(annotation_matrix)

        elif dram_output:
            logging.info('Parsing DRAM output')
            if ko:
                parse_key = "ko_id"
            elif pfam:
                parse_key = "pfam_hits"
            elif cazy:
                parse_key = "cazy_ids"

            headers, tables = Parser.parse_dram_output(dram_output)
            long = Parser.merge_counts_long(headers, tables, key=parse_key).to_pandas()
            annotations = long[parse_key].unique().tolist()

            annotations_dict = (
                long.groupby("sample")
                    .apply(lambda g: dict(zip(g[parse_key], g["count"].astype(float))))
                    .to_dict()
            )

        elif emapper_output:
            logging.info('Parsing emapper output')
            if ko:
                parse_key = self.KEGG
            elif cog:
                parse_key = self.COG
            elif go:
                parse_key = self.GO
            elif eggnog:
                parse_key = self.EGGNOG
            elif pfam:
                parse_key = self.PFAM
            elif ec:
                parse_key = self.EC

            genome_ids, tables = Parser.parse_emapper_output(emapper_output, parse_key)
            long = Parser.merge_counts_long(genome_ids, tables, key="annotation").to_pandas()
            annotations = long["annotation"].unique().tolist()
            annotations_dict = (
                long.groupby("sample")
                    .apply(lambda g: dict(zip(g["annotation"], g["count"].astype(float))))
                    .to_dict()
            )

        if min_prevalence > 0:
            annotations_dict = self.filter_by_prevalence(annotations_dict, min_prevalence)
            annotations = list(set(chain(*[list(x.keys()) for x in annotations_dict.values()])))

        annotation_type = (
            parse_key if emapper_output
            else self.check_annotation_type(annotations)
        )
        
        if abundances_path:
            logging.info('Running abundances pipeline')
            logging.info('Parsing sample abundance')
            abundances_dict, _, _ = Parser.parse_simple_matrix(abundances_path)

            logging.info('Parsing sample metadata')
            _, _, ab_attribute_dict = Parser.parse_metadata_matrix(abundance_metadata_path)

            test = Test(annotations_dict,
                        None,
                        annotation_type,
                        threshold,
                        multi_test_correction,
                        processes,
                        database)

            weighted_abundance \
                = self.weight_annotation_matrix(abundances_dict,
                                                annotations_dict,
                                                ab_attribute_dict,
                                                annotations)
            results = test.test_weighted_abundances(weighted_abundance, annotations)
            genome_list = list(annotations_dict.keys())
            
            for result in results:
                test_result_lines, test_result_output_file = result
                test_result_output_path = os.path.join(output_directory, test_result_output_file)
                Writer.write(test_result_lines, test_result_output_path)

        else:
            
            logging.info('Parsing metadata: %s' % metadata_path)
            metadata, metadata_value_lists, attribute_dict \
                = Parser.parse_metadata_matrix(metadata_path)
                
            logging.info("Comparing sets of genomes")
            combination_dict = dict()
            
            for combination in product(*list([metadata_value_lists])):
                genome_list = list()

                for genome, attributes in metadata.items():

                    for feature in combination:

                        if feature in attributes:
                            genome_list.append(genome)

                combination_dict['_'.join(combination)] = genome_list

            test = Test(annotations_dict, combination_dict, annotation_type, threshold, multi_test_correction, processes, database)
            results = test.test_pipeline(attribute_dict, tree_path=tree_path)

            if decompose:
                logging.info('Running NMF decomposition')
                nmf_results = test.nmf_decompose(n_components=n_components,
                                                 select_components=select_components)
                results.extend(nmf_results)

            for result in results:
                test_result_lines, test_result_output_file = result
                test_result_output_path = os.path.join(output_directory, test_result_output_file)
                Writer.write(test_result_lines, test_result_output_path)

            raw_proportions_output_lines = \
                self.calculate_portions(annotations,
                                        combination_dict,
                                        annotations_dict,
                                        genome_list,
                                        proportions_cutoff)
            Writer.write(raw_proportions_output_lines, os.path.join(output_directory, self.PROPORTIONS))

        if gff_files:
            logging.info("Searching for co-located clusters of genes")
            synteny_results_output_lines, synteny_results_path \
                = syntenysearcher.search_for_blocks(results[0][0],
                                                    gene_positions,
                                                    gene_order,
                                                    metadata,
                                                    intergenic_distance,
                                                    subblock_size,
                                                    operon_mismatch_cutoff,
                                                    operon_match_score_cutoff)
            Writer.write(synteny_results_output_lines, os.path.join(output_directory, synteny_results_path))

            logging.info("Flagging mobile element proximity")
            me_output_lines, me_output_path \
                = syntenysearcher.flag_mobile_elements(results[0][0],
                                                       gene_order,
                                                       metadata,
                                                       me_distance)
            Writer.write(me_output_lines, os.path.join(output_directory, me_output_path))

        if annotation_type == self.KEGG:
            logging.info('Finding module completeness in differentially abundant KOs')

            for result_file in os.listdir(output_directory):

                if(result_file.endswith("fisher.tsv") or result_file.endswith("cdf.tsv")):
                    module_output, prefix = self.module_completeness(database, os.path.join(output_directory, result_file), pval_cutoff)
                    Writer.write(module_output, os.path.join(output_directory, prefix +'_'+ self.MODULE_COMPLETENESS))

class Test(Enrichment):
    __test__ = False

    FISHER_HEADER = [['annotation', 'group_1', 'group_2', 'enriched_in', 'group_1_true', 'group_1_false',
                      'group_2_true', 'group_2_false', 'odds_ratio', 'pvalue', 'corrected_pvalue', 'description']]

    MANNWHITNEYU_HEADER = [['annotation', 'group_1', 'group_2', 'enriched_in', 'group_1_mean', 'group_2_mean',
                            'fold_change', 'effect_size', 'U_statistic', 'pvalue', 'corrected_pvalue', 'description']]

    ZSCORE_HEADER = [['annotation', 'group_1', 'group_2', 'enriched_in', 'group_1_mean', 'group_1_sd',
                      'group_2_count', 'z_score', 'pvalue', 'corrected_pvalue', 'description']]

    KW_HEADER = [['annotation', 'groups', 'H_statistic', 'pvalue', 'corrected_pvalue', 'description']]
    KW_OUTPUT = 'kruskal_wallis.tsv'

    INDVAL_HEADER = [['annotation', 'indicator_group', 'indval', 'specificity', 'fidelity', 'pvalue', 'corrected_pvalue', 'description']]
    INDVAL_OUTPUT = 'indval_results.tsv'

    NMF_LOADINGS_OUTPUT = 'nmf_loadings.tsv'
    NMF_SCORES_OUTPUT = 'nmf_scores.tsv'
    NMF_MWU_HEADER = [['component', 'group_1', 'group_2', 'enriched_in', 'group_1_mean', 'group_2_mean',
                        'fold_change', 'effect_size', 'U_statistic', 'pvalue', 'corrected_pvalue']]
    NMF_MWU_OUTPUT = 'nmf_component_mwu.tsv'

    PHYLO_HEADER = [['annotation', 'concordant_pairs', 'discordant_pairs',
                     'pvalue', 'corrected_pvalue', 'description']]
    PHYLO_OUTPUT = 'phylo_corrected.tsv'

    PA = 'presence_absence'
    IVG_OUTPUT = 'ivg_results.cdf.tsv'
    GENE_FISHER_OUTPUT = 'gvg_results.fisher.tsv'
    GVG_OUTPUT = 'gvg_results.mannwhitneyu.tsv'

    mtc_dict = {'b': 'Bonferroni',
                's': 'Sidak',
                'h': 'Holm',
                'hs': 'Holm-Sidak',
                'sh': 'Simes-Hochberg',
                'ho': 'Hommel',
                'fdr_bh': 'FDR Benjamini-Hochberg',
                'fdr_by': 'FDR Benjamini-Yekutieli',
                'fdr_tsbh': 'FDR 2-stage Benjamini-Hochberg',
                'fdr_tsbky': 'FDR 2-stage Benjamini-Krieger-Yekutieli',
                'fdr_gbs': 'FDR adaptive Gavrilov-Benjamini-Sarkar'}

    def __init__(self, genome_annotations, groups,
                 annotation_type, threshold, multi_test_correction,
                 processes, database):
        '''
        Collects functions to count and test differential abundance among groups of genomes.
        '''

        self.threshold = threshold
        self.multi_test_correction = multi_test_correction
        self.annotation_type = annotation_type
        self.groups = groups
        self.pool = None
        if processes and processes > 1:
            self.pool = mp.Pool(processes=processes)
        self._database = database
        self._descriptions_loaded = False

        if annotation_type==self.PFAM:
            self.genome_annotations = dict()
            for key, item in genome_annotations.items():
                self.genome_annotations[key] = {key.split('.')[0]:entry for key,entry in item.items()}
        else:
            self.genome_annotations = genome_annotations

    def _map(self, func, items):
        if self.pool is None:
            return [func(item) for item in items]
        return self.pool.map(func, items)

    def test_chooser(self, groups):
        groups = [len(x) for x in groups]

        # Enrichment
        if any(group == 1 for group in groups):
            enrichment_test = self.PA
            overrepresentation_test = stats.norm.cdf

        else:
            enrichment_test = stats.fisher_exact
            overrepresentation_test = stats.mannwhitneyu

        return enrichment_test, overrepresentation_test

    def correct_multi_test(self, pvalues):
        logging.info('Applying multi-test correction using the %s method' % (self.mtc_dict[self.multi_test_correction]) )
        corrected_pvals \
            = sm.multipletests(pvalues,
                               alpha        = self.threshold,
                               method       = self.multi_test_correction,
                               returnsorted = False,
                               is_sorted    = False)[1]

        return corrected_pvals

    def count(self, annotation, group, freq):

        if freq:
            group_true = list()

        else:
            group_true = 0

        group_false = 0

        for genome in self.groups[group]:

            if annotation in self.genome_annotations[genome]:
                if freq:
                    group_true.append(self.genome_annotations[genome][annotation])
                else:
                    if self.genome_annotations[genome][annotation]>0.0:
                        group_true+=1
                    else:
                        group_false+=1

            else:

                if freq:
                    group_true.append(0)

                else:
                    group_false+=1

        return group_true, group_false

    def gene_frequencies(self, group_1, group_2, freq=False):

        res_list    = list()
        annotations = sorted(set(chain(*self.genome_annotations.values())))

        for annotation in annotations:
            passed = True
            group_1_true, group_1_false \
                = self.count(annotation,
                              group_1,
                              freq)
            group_2_true, group_2_false \
                = self.count(annotation,
                              group_2,
                              freq)

            if freq:
                if(len([x for x in group_1_true if x!='0'])==0 and
                    len([x for x in group_2_true if x!='0'])==0 ):
                    passed = False
            else:
                if(group_1_true==0 and group_2_true==0):
                    passed = False

            if passed:
                res_list.append([annotation, group_1, group_2, [group_1_true, group_1_false], [group_2_true, group_2_false]])

        return res_list

    def corrected_pvals(self, output_lines):
        pvalues = [output_line[-1] for output_line in output_lines]
        corrected_pvalues = self.correct_multi_test(pvalues)

        return corrected_pvalues

    def _load_descriptions(self):
        if not self._descriptions_loaded:
            self.k = self._database.k()
            self.tigrfamdescription = self._database.tigrfamdescription()
            self.pfam2description = self._database.pfam2description()
            self.ec2description = self._database.ec2description()
            self._descriptions_loaded = True

    def add_descriptions(self, output_lines):
        self._load_descriptions()

        if self.annotation_type == self.KEGG:
            desc = self.k

        if self.annotation_type == self.CAZY:
            desc = None

        if self.annotation_type == self.TIGRFAM:
            desc = self.tigrfamdescription

        if self.annotation_type == self.PFAM:
            desc = self.pfam2description

        if self.annotation_type == self.EC:
            desc = self.ec2description

        if self.annotation_type == self.OTHER:
            desc = None

        if self.annotation_type in (self.COG, self.GO, self.EGGNOG):
            desc = None

        for line in output_lines:
            annotation = line[0]

            if desc:

                if annotation in desc:
                    line.append(desc[annotation])

                else:
                    line.append("NA")

            else:
                line.append("NA")

        return output_lines

    def test_weighted_abundances(self,
                                 weighted_abundance,
                                 annotations):

        logging.info('Calculating enrichment across samples using Mann-Whitney U test')
        results = list()

        for combination in combinations(weighted_abundance, 2):
            prefix = '_vs_'.join(
                [sorted(combination)[0], sorted(combination)[1]]).replace(' ', '_')
            res_list = list()

            for annotation in annotations:
                group_1 = weighted_abundance[combination[0]][annotation]
                group_2 = weighted_abundance[combination[1]][annotation]
                gene_count = [annotation, combination[0],
                              combination[1], [group_1], [group_2]]
                res_list.append(gene_count)
            output_lines = self._map(mannwhitneyu_calc, res_list)

            for idx, corrected_pval in enumerate(self.corrected_pvals(output_lines)):
                output_lines[idx].append(str(corrected_pval))
            output_lines = self.add_descriptions(output_lines)
            output_lines = self.MANNWHITNEYU_HEADER + output_lines
            results.append([output_lines, prefix + '_' + self.GVG_OUTPUT])

        return results

    def kruskal_wallis_frequencies(self):
        group_names = list(self.groups.keys())
        annotations = sorted(set(chain(*self.genome_annotations.values())))
        res_list = []
        for annotation in annotations:
            group_values = [self.count(annotation, g, freq=True)[0] for g in group_names]
            if all(sum(v) == 0 for v in group_values):
                continue
            res_list.append([annotation, group_names, group_values])
        return res_list

    def _build_annotation_matrix(self):
        '''Build genomes × annotations matrix. Returns (all_genomes, annotations, X).'''
        all_genomes = list(self.genome_annotations.keys())
        annotations = sorted(set(chain(*self.genome_annotations.values())))
        X = np.array([[self.genome_annotations[g].get(a, 0.0) for a in annotations]
                      for g in all_genomes], dtype=float)
        return all_genomes, annotations, X

    def _select_k_cophenetic(self, X, k_max=30, n_runs=20):
        '''Select NMF k via cophenetic correlation coefficient.'''
        n_genomes = X.shape[0]
        k_range = range(2, min(k_max + 1, n_genomes))
        best_k = 2
        best_score = -1.0

        for k in k_range:
            consensus = np.zeros((n_genomes, n_genomes))
            for _ in range(n_runs):
                model = NMF(n_components=k, max_iter=500)
                W = model.fit_transform(X)
                labels = np.argmax(W, axis=1)
                for i in range(n_genomes):
                    for j in range(n_genomes):
                        if labels[i] == labels[j]:
                            consensus[i, j] += 1
            consensus /= n_runs
            dist = pdist(consensus, metric='euclidean')
            Z = linkage(dist, method='average')
            c, _ = cophenet(Z, dist)
            logging.info(f'k={k}: cophenetic correlation = {c:.4f}')
            if c > best_score:
                best_score = c
                best_k = k

        logging.info(f'Selected k={best_k} (cophenetic r={best_score:.4f})')
        return best_k

    def nmf_decompose(self, n_components=None, select_components=False):
        '''
        Decompose annotation matrix via NMF. Returns list of [lines, filename] results.

        Parameters
        ----------
        n_components        - int or None. Fixed k; overrides heuristic and cophenetic selection.
        select_components   - bool. Use cophenetic correlation to select k automatically.
        '''
        all_genomes, annotations, X = self._build_annotation_matrix()
        n_genomes = X.shape[0]

        if n_components is not None:
            k = n_components
        elif select_components:
            logging.info('Selecting NMF k via cophenetic correlation (this may take a while)')
            k = self._select_k_cophenetic(X)
        else:
            k = max(2, int(math.sqrt(n_genomes / 2)))

        logging.info(f'Running NMF decomposition with k={k} components on '
                     f'{n_genomes} genomes × {len(annotations)} annotations')

        model = NMF(n_components=k, random_state=42, max_iter=500)
        W = model.fit_transform(X)   # n_genomes × k
        H = model.components_        # k × n_annotations

        results = []

        # Component loadings: rows = components, cols = annotations
        loading_lines = [['component'] + list(annotations)]
        for i, row in enumerate(H):
            loading_lines.append([f'component_{i + 1}'] + [str(round(float(v), 6)) for v in row])
        results.append([loading_lines, self.NMF_LOADINGS_OUTPUT])

        # Genome component scores: rows = genomes, cols = components
        score_lines = [['genome'] + [f'component_{i + 1}' for i in range(k)]]
        for genome, row in zip(all_genomes, W):
            score_lines.append([genome] + [str(round(float(v), 6)) for v in row])
        results.append([score_lines, self.NMF_SCORES_OUTPUT])

        # MWU on component scores between all group pairs
        genome_idx = {g: i for i, g in enumerate(all_genomes)}
        group_names = list(self.groups.keys())
        mwu_lines = []
        for component_idx in range(k):
            comp_name = f'component_{component_idx + 1}'
            for g1, g2 in combinations(group_names, 2):
                g1_scores = [W[genome_idx[g], component_idx] for g in self.groups[g1] if g in genome_idx]
                g2_scores = [W[genome_idx[g], component_idx] for g in self.groups[g2] if g in genome_idx]
                row = mannwhitneyu_calc([comp_name, g1, g2, [g1_scores], [g2_scores]])
                mwu_lines.append(row)

        if mwu_lines:
            for idx, cpval in enumerate(self.corrected_pvals(mwu_lines)):
                mwu_lines[idx].append(str(cpval))
            mwu_lines = self.NMF_MWU_HEADER + mwu_lines
            results.append([mwu_lines, self.NMF_MWU_OUTPUT])

        return results

    def _get_phylo_pairs(self, tree, g1_genomes, g2_genomes):
        '''
        Traverse tree internal nodes and collect cross-clade pairs where one
        genome is in g1 and the other is in g2. Returns list of (g1_genome,
        g2_genome) tuples — one per informative cross-clade pairing.
        '''
        all_genomes = set(g1_genomes) | set(g2_genomes)
        g1_set = set(g1_genomes)
        g2_set = set(g2_genomes)
        pairs = set()
        for node in tree.internal_nodes():
            children = node.child_nodes()
            if len(children) < 2:
                continue
            child_leaves = []
            for child in children:
                leaves = {n.taxon.label.strip() for n in child.leaf_iter()
                          if n.taxon and n.taxon.label.strip() in all_genomes}
                child_leaves.append(leaves)
            for i in range(len(child_leaves)):
                for j in range(i + 1, len(child_leaves)):
                    for l in child_leaves[i]:
                        for r in child_leaves[j]:
                            if l in g1_set and r in g2_set:
                                pairs.add((l, r))
                            elif r in g1_set and l in g2_set:
                                pairs.add((r, l))
        return list(pairs)

    def phylo_correction(self, tree_path, group_1, group_2, n_permutations=999):
        '''
        Run Scoary1-style phylogenetic paired comparison for all annotations,
        testing whether group associations survive correction for shared ancestry.

        Parameters
        ----------
        tree_path       - str. Path to Newick tree file.
        group_1, group_2 - str. Group names in self.groups.
        n_permutations  - int. Number of within-pair permutations (default 999).

        Returns [lines, filename] or None if no informative pairs found.
        '''
        tree = dendropy.Tree.get(path=tree_path, schema='newick',
                                 preserve_underscores=True)
        g1_genomes = self.groups[group_1]
        g2_genomes = self.groups[group_2]
        pairs = self._get_phylo_pairs(tree, g1_genomes, g2_genomes)

        if not pairs:
            logging.warning(
                f'No phylogenetically informative pairs for {group_1} vs {group_2}. '
                'Check that tree tip labels match genome names.')
            return None

        logging.info(f'Phylo correction: {len(pairs)} paired comparisons '
                     f'({group_1} vs {group_2})')

        annotations = sorted(set(chain(*self.genome_annotations.values())))
        calc_input = [[ann, pairs, self.genome_annotations, n_permutations]
                      for ann in annotations]

        lines = [x for x in self._map(phylo_pairs_calc, calc_input) if x is not None]
        if not lines:
            return None

        for idx, cpval in enumerate(self.corrected_pvals(lines)):
            lines[idx].append(str(cpval))
        lines = self.add_descriptions(lines)
        return [self.PHYLO_HEADER + lines, self.PHYLO_OUTPUT]

    def indval_frequencies(self, n_permutations=999):
        annotations = sorted(set(chain(*self.genome_annotations.values())))
        res_list = []
        for annotation in annotations:
            all_group_abundances = {g: self.count(annotation, g, freq=True)[0] for g in self.groups}
            if all(sum(float(v) for v in vals) == 0 for vals in all_group_abundances.values()):
                continue
            for focal_group in self.groups:
                res_list.append([annotation, focal_group,
                                 all_group_abundances[focal_group],
                                 all_group_abundances,
                                 n_permutations])
        return res_list

    def test_pipeline(self, group_dict, tree_path=None):
        results = list()

        if len(group_dict) > 2:
            logging.info('Testing multi-group over-representation using Kruskal-Wallis test')
            kw_input = self.kruskal_wallis_frequencies()
            kw_lines = [x for x in self._map(kruskal_wallis_calc, kw_input) if x is not None]
            if kw_lines:
                for idx, cpval in enumerate(self.corrected_pvals(kw_lines)):
                    kw_lines[idx].append(str(cpval))
                kw_lines = self.add_descriptions(kw_lines)
                kw_lines = self.KW_HEADER + kw_lines
                results.append([kw_lines, self.KW_OUTPUT])

        logging.info('Calculating IndVal indicator scores')
        indval_input = self.indval_frequencies()
        indval_lines = [x for x in self._map(indval_calc, indval_input) if x is not None]
        if indval_lines:
            for idx, cpval in enumerate(self.corrected_pvals(indval_lines)):
                indval_lines[idx].append(str(cpval))
            indval_lines = self.add_descriptions(indval_lines)
            indval_lines = self.INDVAL_HEADER + indval_lines
            results.append([indval_lines, self.INDVAL_OUTPUT])

        for combination in combinations(group_dict, 2):
            enrichment_test, overrepresentation_test = self.test_chooser( [group_dict[member] for member in combination] )
            prefix = '_vs_'.join([sorted(combination)[0], sorted(combination)[1]]).replace(' ', '_')

            logging.info('Comparing gene frequency among groups: %s', ', '.join(combination))

            if enrichment_test == stats.fisher_exact:
                logging.info('Testing gene enrichment using Fisher\'s exact test')
                gene_count = self.gene_frequencies(*combination)
                output_lines = self._map(gene_fisher_calc, gene_count)

                for idx, corrected_pval in enumerate(self.corrected_pvals(output_lines)):
                    output_lines[idx].append(str(corrected_pval))

                header = self.FISHER_HEADER
                output = self.GENE_FISHER_OUTPUT
                output_lines = self.add_descriptions(output_lines)
                output_lines = header + output_lines
                results.append([output_lines, prefix + '_' + output])

            elif enrichment_test == self.PA:
                logging.info('enrichment statistics not possible with only one genome to compare')
                logging.info('See prevalence matrix for unique genes in groups')

            logging.info('Comparing gene over-representation among genomes')

            if(overrepresentation_test == stats.mannwhitneyu):
                logging.info('Testing over-representation using Mann-Whitney U test')
                gene_count = self.gene_frequencies(*combination, True)
                output_lines = self._map(mannwhitneyu_calc, gene_count)

                for idx, corrected_pval in enumerate(self.corrected_pvals(output_lines)):
                    output_lines[idx].append(str(corrected_pval))

                header = self.MANNWHITNEYU_HEADER
                output = self.GVG_OUTPUT

            elif overrepresentation_test == stats.norm.cdf:
                logging.info('Testing over-representation using Z score test')
                gene_count = self.gene_frequencies(*combination, True)
                output_lines = self._map(zscore_calc, gene_count)
                output_lines = [x for x in output_lines if x]

                for idx, corrected_pval in enumerate(self.corrected_pvals(output_lines)):
                    output_lines[idx].append(str(corrected_pval))

                header = self.ZSCORE_HEADER
                output = self.IVG_OUTPUT

            output_lines = self.add_descriptions(output_lines)
            output_lines = header + output_lines
            results.append([output_lines, prefix +'_'+ output])

            if tree_path:
                logging.info(f'Running phylogenetic correction for {combination[0]} vs {combination[1]}')
                phylo_result = self.phylo_correction(tree_path, combination[0], combination[1])
                if phylo_result:
                    phylo_lines, phylo_output = phylo_result
                    results.append([phylo_lines, prefix + '_' + phylo_output])

        return results
