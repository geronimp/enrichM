#!/usr/bin/env python
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

__author__      = "Joel Boyd"
__copyright__   = "Copyright 2017"
__credits__     = ["Joel Boyd"]
__license__     = "GPL3"
__version__     = "0.0.7"
__maintainer__  = "Joel Boyd"
__email__       = "joel.boyd near uq.net.au"
__status__      = "Development"

###############################################################################

# System
import logging
import os
import random
import re

import numpy as np
import multiprocessing as mp
import statsmodels.sandbox.stats.multicomp as sm

from scipy import stats
from itertools import product, combinations, chain

# Local
from enrichm.parse_annotate import ParseAnnotate
from enrichm.databases import Databases
from enrichm.draw_plots import Plot
from enrichm.comparer import Compare

################################################################################

def gene_fisher_calc(x):
    
    ko, group_1, group_2 = x[0], x[1], x[2]
    
    dat = x[3:]
    
    if (dat[0][0]>0 or dat[1][0]>0) and (dat[0][1]>0 or dat[1][1]>0):
        score, pval = stats.fisher_exact(dat)
    else:
        score, pval = 'nan', 1.0

    return [group_1, group_2, ko] + dat[0] + dat[1] + [score, pval]

def group_mannwhitneyu_calc(x):
    # Mann Whitney U test
    module, group_1, group_2, group_1_module_kos, group_2_module_kos, module_list = x

    if(len(group_1_module_kos)>1 and len(group_2_module_kos)>1):
        if(sum(group_1_module_kos)>0 and sum(group_2_module_kos)>0):
            if(len(set(group_1_module_kos)) == 1
                or
               len(set(group_2_module_kos)) == 1):
                mw_t_stat, mw_p_value = 'NA', 1
            else:
                group_1_module_kos = np.array(group_1_module_kos)       
                group_2_module_kos = np.array(group_2_module_kos)           
                mw_t_stat, mw_p_value = \
                    stats.mannwhitneyu(group_1_module_kos,
                                             group_2_module_kos)
            
            return [module, group_1, group_2, str(np.mean(group_1_module_kos)),
                    str(np.mean(group_2_module_kos)), str(len(module_list)), 
                    mw_t_stat, mw_p_value]

def zscore_calc(x):
    
    module, group, genome, genome_comp, reference_group_comp_list, module_list, description = x
    
    if genome_comp>0:
        reference_group_comp_sd \
            = np.std(reference_group_comp_list, axis=0)
        reference_group_comp_mean \
            = np.mean(reference_group_comp_list, axis=0)
        
        if (genome_comp-reference_group_comp_mean)>0:
            z_score = (genome_comp-reference_group_comp_mean) / reference_group_comp_sd
            p_value = 2-2*stats.norm.cdf(z_score)
            
            return [module, group, genome, str(reference_group_comp_mean), 
                    str(genome_comp), module_list, str(z_score), 
                    str(p_value), description]

################################################################################

class Enrichment:
    '''
    Class contianing various functions to calculate the enrichment of funcitons
    and genes between two groups of genomes...
    '''
    
    TIGRFAM                 = "tigrfam"
    PFAM                    = "pfam"
    KEGG                    = "kegg"

    def __init__(self):

        self.BACKGROUND              = 'background'
        self.MATRIX_SUFFIX           = '_enrichment_matrix.tsv'
        self.SAMPLE_MATRIX_SUFFIX    = '_sample_enrichment_matrix.tsv'
        self.COMPARE_SUFFIX          = '_compare_matrix.tsv'    
        self.TIGRFAM_PREFIX          = 'TIGR'
        self.PFAM_PREFIX             = 'PF'
        self.KEGG_PREFIX             = 'K'
        self.PROPORTIONS             = 'proportions.tsv'
        self.UNIQUE_TO_GROUPS        = 'unique_to_groups.tsv'
        
        self.taxonomy_index_dictionary = {"d__":0, "p__":1, "c__":2,
                                          "o__":3, "f__":4, "g__":5, "s__":6}    

    def _parse_matrix(self, matrix_file_io, colnames):
        for line in matrix_file_io:
            sline = line.strip().split('\t')
            rowname, entries = sline[0], sline[1:]
            for colname, entry in zip(colnames, entries):
                yield colname, entry, rowname

    def parse_metadata_matrix(self, matrix_path):
        '''        
        Parameters
        ----------
        matrix_path : String. Path to file containing a matrix of genome rownames.        
        '''

        matrix_file_io  = open(matrix_path)
        cols_to_rows    = dict()
        nr_values       = set()
        attribute_dict  = dict()
        
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

    def _parse_annotation_matrix(self, matrix_path):
        '''        
        Parameters
        ----------
        matrix_path : String. Path to file containing a matrix of genome rownames.        
        '''

        matrix_file_io  = open(matrix_path)
        colnames        = matrix_file_io.readline().strip().split('\t')[1:]
        cols_to_rows    = {genome_name:set() for genome_name in colnames}
        rownames        = set()

        for genome_name, entry, rowname \
                        in self._parse_matrix(matrix_file_io, colnames):
            
            rownames.add(rowname)
            if float(entry) > 0:
                cols_to_rows[genome_name].add(rowname)

        return cols_to_rows, rownames, colnames

    def check_annotation_type(self, annotations):
        '''
        Takes a random sample of the rownames from the input matrix
        and based on the characters they start with, determines the 
        type of annotation being used.
        
        Parameters
        ----------
        annotaitons     - List. A list of strings, each a rowname 
                          from the original input matrix

        Output
        ------
        The annotation type
        '''
        sample = random.sample(annotations, 1)[0]
        if sample.startswith(self.TIGRFAM_PREFIX):
            return self.TIGRFAM
        elif sample.startswith(self.KEGG_PREFIX):
            return self.KEGG
        else:
            return self.PFAM

    def calculate_portions(self, modules, combination_dict, annotations_dict, genome_list, proportions_cutoff):
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
        ### ~ TODO: Add a portion of genome column too?

        raw_proportions_output_lines        = [['Module'] + combination_dict.keys()]
        enriched_proportions_output_lines   = [['Module'] + combination_dict.keys()]

        for module in modules:
            
            module_values               = {}
            raw_proportions_output_line = [module]
            
            for group_name, genome_list in combination_dict.items():
                if len(genome_list)>0:
                    coverage = len([genome for genome in genome_list 
                                    if module in annotations_dict[genome]])
                    total    = float(len(genome_list))
                    entry    = coverage/total
                    module_values[group_name] = entry
                    raw_proportions_output_line.append(str(entry))
                else:
                    raw_proportions_output_line.append('0.0')
            
            if max(module_values.values())>0:
                groups = combination_dict.keys()
                for group in groups:
                    group_val = module_values[group]
                    compare_groups = []
                    if group_val>0:
                        for other_group in groups:
                            if other_group!=group:
                                other_group_val = module_values[other_group]
                                if other_group_val > 0:
                                    compare_value = group_val/other_group_val
                                else:
                                    compare_value = float("inf")
                                compare_groups.append(compare_value)
                    
                    if all([x>proportions_cutoff for x in compare_groups]):                    
                        pass 

            raw_proportions_output_lines.append(raw_proportions_output_line)

        return raw_proportions_output_lines

    def _write(self, output_lines_list, output_path):
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
            for line in output_lines_list:
                string = '\t'.join([str(x) for x in line]) + '\n'
                out_io.write(string)
    
    def parse_taxonomy(self, taxonomy, taxonomy_dict):
        genomes_set = set()
        rank = taxonomy[:3]
        
        if rank in self.taxonomy_index_dictionary:
            rank_index = self.taxonomy_index_dictionary[rank]
        else:
            raise Exception("Rank doesnt exist (%s) Does your taxonomy have a GTDB rank prefix?" % (taxonomy))
        
        for genome_id, taxonomy_list in taxonomy_dict.items():

            if taxonomy_list[rank_index] == taxonomy:
                genomes_set.add(genome_id)

        return genomes_set 

    def do(# Input options
           self, annotate_output, metadata_path, modules, abundances, 
           # Runtime options
           do_all, do_ivi, do_gvg, do_ivg, pval_cutoff, proportions_cutoff, threshold, 
           multi_test_correction, taxonomy, batchfile, processes, ko, pfam, tigrfam, 
           hypothetical, cazy,
           # Output options
           output_directory):

        p  = Plot()
        c  = Compare()
        d  = Databases()
        
        logging.info('Parsing inputs')
        pa = ParseAnnotate(annotate_output, processes)

        if ko:
            annotation_matrix = pa.ko
        elif pfam:
            annotation_matrix = pa.pfam
        elif tigrfam:
            annotation_matrix = pa.tigrfam
        elif hypothetical:
            annotation_matrix = pa.hypothetical_cluster
        elif cazy:
            annotation_matrix = pa.cazy

        logging.info('Parsing annotations')

        annotations_dict, modules, genomes \
                    = self._parse_annotation_matrix(annotation_matrix)

        if (taxonomy or batchfile):
            genomes_set = set()
            if taxonomy:
                genomes_set = genomes_set.union(self.parse_taxonomy(taxonomy, d.taxonomy))
            if batchfile:
                batchfile_metadata, batchfile_metadata_value_lists, batchfile_attribute_dict \
                            = self.parse_metadata_matrix(batchfile)
                genomes_set = genomes_set.union(set(batchfile_metadata.keys()))

        if modules:
            logging.info('Limiting to %i modules' % len(modules))
            modules = modules

        logging.info('Parsing metadata')
        metadata, metadata_value_lists, attribute_dict \
                    = self.parse_metadata_matrix(metadata_path)

        # Load pickles
        pa.genome_objects = pa.parse_pickles(pa.genome_pickle_file_path, metadata.keys())
        reference_genomes = pa.parse_pickles(d.GTDB_DIR, genomes_set)

        reference_genome_annotations = {genome.name.replace('_gene', ''):set(genome.ko_dict.keys()) for genome in reference_genomes}
        annotations_dict.update(reference_genome_annotations)

        if taxonomy:
            attribute_dict[taxonomy] = set()
            for genome in reference_genome_annotations.keys():
                genomes.append(genome)
                metadata[genome] = set([taxonomy])
                attribute_dict[taxonomy].add(genome)
            metadata_value_lists.add(taxonomy)
        elif batchfile:
            found_gtdb_genomes = set([x.name.replace('_gene', '') for x in reference_genomes])
            for group in batchfile_metadata_value_lists:
                attribute_dict[group] = set()
            for genome in found_gtdb_genomes:
                genome_group = batchfile_metadata[genome]
                for g in attribute_dict:
                    attribute_dict[g].add(genome)
                metadata[genome] = genome_group
                genomes.append(genome)
                
            metadata_value_lists = metadata_value_lists.union(batchfile_metadata_value_lists)

        logging.info("Comparing sets of genomes")        
        combination_dict = dict()
        for combination in product(*list([metadata_value_lists])):
            genome_list = list()
            for genome, attributes in metadata.items():
                for feature in combination:
                    if feature in attributes:
                        genome_list.append(genome)  
            combination_dict['_'.join(combination)] = genome_list

        annotation_type = self.check_annotation_type(modules)
        t = Test(annotations_dict,
                 modules,
                 genomes,
                 combination_dict,
                 annotation_type,
                 threshold,
                 multi_test_correction,
                 pval_cutoff, 
                 processes,
                 d)

        if do_all:
            do_gvg = True
            do_ivi = True
            do_ivg = True
        
        if do_gvg:
            gvg = list(combinations(combination_dict.keys(), 2)) # For Fisher's exact test
            if len(gvg)==0:
                logging.info('No gvg comparison possible')
                gvg=None
        else: 
            gvg = None
        
        if do_ivi:
            ivi = list(combinations(chain(*combination_dict.values()), 2)) # For T-test
            if len(ivi)==0:
                logging.info('No ivi comparison possible')
                ivi=None
        else:
            ivi = None

        if do_ivg:
            ivg = dict() # For Z score test
            for group_name, genome_list in combination_dict.items():
                if len(genome_list)>1:
                    ivg[group_name] = list()
                    for idx, genome in enumerate(genome_list):
                        other=genome_list[:idx] + genome_list[idx+1:]
                        ivg[group_name].append((genome, other))
            if len(ivg)==0:
                logging.info('No ivg comparison possible')
                ivg=None
        else:
            ivg = None

        results = t.do(ivi, ivg, gvg)

        for test_result_lines, test_result_output_file in results:
            test_result_output_path = os.path.join(output_directory,
                                                   test_result_output_file)
            self._write(test_result_lines, test_result_output_path)
        
        raw_portions_path \
            = os.path.join(output_directory, self.PROPORTIONS)
        unique_to_groups_path \
            = os.path.join(output_directory, self.UNIQUE_TO_GROUPS)
        raw_proportions_output_lines \
            = self.calculate_portions(modules, combination_dict, annotations_dict, genome_list, proportions_cutoff)

        self._write(raw_proportions_output_lines, raw_portions_path)

        logging.info('Generating summary plots')

        if annotation_type==self.KEGG:
            p.draw_barplots(os.path.join(output_directory, t.GENE_FISHER_OUTPUT), pval_cutoff, output_directory)

        p.draw_pca_plot(annotation_matrix, metadata_path, output_directory)

        #c.do(pa.genome_objects,
        #     attribute_dict,
        #     output_directory)

class Test(Enrichment):

    def __init__(self, genome_annotations, modules, genomes, groups, 
                 annotation_type, threshold, multi_test_correction, 
                 pval_cutoff, processes, d):
        '''
        
        Parameters
        ----------
        
        Output
        ------
        '''
               
        self.IVI_OUTPUT             = 'ivi_results.fisher.tsv'
        self.IVG_OUTPUT             = 'ivg_results.tsv'
        self.GENE_FISHER_OUTPUT     = 'gvg_results.fisher.tsv'
        self.GVG_OUTPUT             = 'gvg_results.mannwhitneyu.tsv'

        self.threshold              = threshold
        self.multi_test_correction  = multi_test_correction
        self.m2def                  = d.m2def
        self.m                      = d.m
        self.clan2pfam              = d.clan2pfam
        self.clan_to_description    = d.clan2name
        self.modules                = modules
        self.genomes                = genomes
        self.annotation_type        = annotation_type
        self.groups                 = groups
        self.pval_cutoff            = pval_cutoff
        self.pool                   = mp.Pool(processes = processes)

        self.mtc_dict = {'b': 'Bonferroni',
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

        if annotation_type==Enrichment.PFAM:
            self.genome_annotations = {}
            for key, item in genome_annotations.items():
                self.genome_annotations[key] = set([x.split('.')[0] for x in item])
        else:
            self.genome_annotations = genome_annotations

        self.pfam                   = dict()
        self.id2pfam                = dict()
        
        for line in open(d.PFAM_CLAN_DB):
            sline       = line.strip().split()
            pfam        = sline[0]
            clan        = sline[1]
            id          = sline[2]
            description = "%s; %s" % (id, ' '.join(sline[2:]))
            self.pfam[pfam]  = description
            self.id2pfam[id] = pfam
    


    def correct_multi_test(self, pvalues):
        logging.info('Applying multi-test correction using the %s method' % (self.mtc_dict[self.multi_test_correction]) )
        corrected_pvals \
            = sm.multipletests(pvalues,
                               alpha        = self.threshold,
                               method       = self.multi_test_correction,
                               returnsorted = False,
                               is_sorted    = False)[1]        

        return corrected_pvals

    def _strip_kegg_definitions(self, definition):
        kos_list = [ko for ko in re.split("[^\w]", definition) if ko]
        return set(kos_list)

    def gather_genome_annotations(self, group, target_annotations):
        
        genome_annotation_list = list()
        
        for genome in self.groups[group]:
            genome_target_list = self.genome_annotations[genome].intersection(target_annotations)
            genome_annotation_list.append(len(genome_target_list))
        
        return genome_annotation_list

    def get_annotations(self):

        if self.annotation_type==Enrichment.TIGRFAM:
            logging.warning('Comparisons are not possible with TIGRFAM because heirarchical classifications are needed (like in KEGG, COG or PFAM).')
        elif self.annotation_type==Enrichment.KEGG:
            iterator = self.m2def
            annotation_description = self.m
        elif self.annotation_type==Enrichment.PFAM:
            iterator = self.clan2pfam
            annotation_description = self.clan_to_description
        return iterator, annotation_description

    def gene_fisher(self):
        """"""
        header      = [['group_1', 'group_2', 'ko', 'group_1_true', 'group_1_false', 'group_2_true',
                        'group_2_false', 'score', 'pvalue', 'corrected_pvalue']]
        kos         = set(chain(*self.genome_annotations.values()))
        res_list    = list()
        
        for group_1, group_2 in combinations(self.groups.keys(), 2):
            
            for ko in kos:

                group_1_true, group_1_false, group_2_true, group_2_false = 0, 0, 0, 0

                if (len(self.groups[group_1])>=5 and len(self.groups[group_2])>=5):
                    for genome_1 in self.groups[group_1]:
                        if ko in self.genome_annotations[genome_1]:
                            group_1_true+=1
                        else:
                            group_1_false+=1
                    for genome_2 in self.groups[group_2]:
                        
                        if ko in self.genome_annotations[genome_2]:
                            group_2_true+=1
                        else:
                            group_2_false+=1

                res_list.append([ko, group_1, group_2, [group_1_true, group_1_false], [group_2_true, group_2_false]])

        output_lines = self.pool.map(gene_fisher_calc, res_list)
        pvalues = [x[-1] for x in output_lines]

        for idx, corrected_pval in enumerate(self.correct_multi_test(pvalues)):
            output_lines[idx].append(str(corrected_pval))
                
        return header + output_lines

    def ivi_fisher(self, ivi):
        ''''''
        header       = [['module', 'genome_1', 'genome_2', 'group_1_count', 'group_2_count',
                         'count', 'stat', 'p_value', 'corrected_p_value', 'description']]
        res_list    = list()
        
        iterator, annotation_description = self.get_annotations()

        for genome_1, genome_2 in ivi:

            for module, module_definition in iterator.items():
                module_kos    = self._strip_kegg_definitions(module_definition)
                
                genome_1_annotations = self.genome_annotations[genome_1]
                genome_2_annotations = self.genome_annotations[genome_2]

                genome_1_comp = module_kos.intersection(genome_1_annotations)
                genome_2_comp = module_kos.intersection(genome_2_annotations)

                row_1 = [len(genome_1_comp), 
                         len(genome_2_comp)]

                row_2 = [(len(genome_1_annotations) - len(genome_1_comp)),
                         (len(genome_2_annotations) - len(genome_2_comp))]

                res_list.append([module, genome_1, genome_2, row_1, row_2])

        output_lines = self.pool.map(gene_fisher_calc, res_list)
        pvalues      = [x[-1] for x in output_lines]

        if len(pvalues)>0:
            corrected_pvalues = self.correct_multi_test(pvalues)
        else:
            corrected_pvalues = ['nan' for x in output_lines]

        output_lines_list = [x[:len(x)-1]+[str(y)]+x[len(x)-1:] for x, y in zip(output_lines, list(corrected_pvalues))]

        return header + output_lines
        
    def ttest(self, gvg):
        
        header      = [['module', 'group_1', 'group_2', 'group_1_mean', 'group_2_mean', 'count', 
                               'mann_whitney_t_stat', 'mann_whitney_p_value', 'mann_whitney_corrected_p_value',
                               'description']]
        res_list    = list()

        iterator, annotation_description = self.get_annotations()
                
        for group_1, group_2 in gvg:

            for module, definition in iterator.items():
                module_list = self._strip_kegg_definitions(definition)
                
                group_1_module_kos = self.gather_genome_annotations(group_1, module_list)
                group_2_module_kos = self.gather_genome_annotations(group_2, module_list)
                
                res_list.append([module, group_1, group_2, group_1_module_kos, group_2_module_kos, module_list])

        output_lines = [x for x in self.pool.map(group_mannwhitneyu_calc, res_list) if x]
        pvalues      = [x[-1] for x in output_lines]      

        for idx, corrected_p in enumerate(self.correct_multi_test(pvalues)):
            output_lines[idx] = [str(x) for x in output_lines[idx]] + [str(corrected_p), annotation_description[output_lines[idx][0]]]

        return header + output_lines

    def zscore(self, ivg):
        header       = [['module', 'group', 'genome', 'group_mean', 'genome_count',
                         'count', 'z_score', 'p_value','corrected_p_value', 'description']]
        res_list     = list()

        iterator, annotation_description = self.get_annotations()

        for group, comparisons in ivg.items():
            for genome, reference_group in comparisons:
                for module, definition in iterator.items():
                    module_list = self._strip_kegg_definitions(definition)
                    
                    genome_comp \
                        = len(self.genome_annotations[genome].intersection(module_list))
                    reference_group_comp_list \
                        = np.array([len(set(self.genome_annotations[reference_genome]).intersection(module_list))
                                    for reference_genome in reference_group])

                    res_list.append([module,
                                     group,
                                     genome,
                                     genome_comp,
                                     reference_group_comp_list,
                                     str(len(module_list)),
                                     annotation_description[module]])
        
        output_lines = [x for x in self.pool.map(zscore_calc, res_list) if x]
        pvalues      = [float(x[-2]) for x in output_lines]

        if len(pvalues)>0:
            corrected_pvalues = self.correct_multi_test(pvalues)
        else:
            corrected_pvalues = ['nan' for x in output_lines]

        output_lines = [x[:len(x)-1]+[str(y)]+x[len(x)-1:] for x, y in zip(output_lines, list(corrected_pvalues))]

        return header + output_lines


    def do(self, ivi, ivg, gvg):

        logging.info('Running enrichent tests')
        results=[]
        
        logging.info('Comparing gene frequency among genomes (presence/absence)')
        results.append( (self.gene_fisher(), self.GENE_FISHER_OUTPUT) )

        if ivi:
            logging.info('Running individual vs individual comparisons')
            results.append( (self.ivi_fisher(ivi), self.IVI_OUTPUT) )
        if ivg:
            logging.info('Running individual vs group comparisons')
            results.append( (self.zscore(ivg), self.IVG_OUTPUT) )
        if gvg:
            logging.info('Running group vs group comparisons')
            results.append( (self.ttest(gvg), self.GVG_OUTPUT) )
        
        return results



