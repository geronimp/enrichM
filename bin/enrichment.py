    #!usr/bin/env python
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
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

###############################################################################

import logging
import os
import random
import re
import scipy.stats
import numpy as np
import statsmodels.sandbox.stats.multicomp as sm
from itertools import product, combinations, chain

# Local
from databases import Databases

###############################################################################

class Matrix:
    '''
    A custom class that parses the matrices the way I want it to. Kinda 
    redundant work actually. Whoops I guess.
    '''
    
    def __init__(self, matrix):
        self.suffix_to_delim = {'.tsv':'\t', '.csv':','}
        self.rownames        = []
        self.colnames        = []
        self.matrix          = {}

        self._parse_matrix(open(matrix))
        
    def _parse_matrix(self, matrix_io):
        suffix = os.path.splitext(matrix_io.name)[-1]
        try:
            delim = self.suffix_to_delim[suffix]
        except:
            raise Exception('Unrecognised file suffix: %s\nPlease provide either a .csv or .tsv file' % suffix)
        headers = matrix_io.readline().strip().split(delim)[1:]
        for header in headers:
            self.matrix[header]={}
            if header in self.colnames:
                raise Exception("Duplicate column names found in %s" % \
                                                         (matrix_io.name))
            self.colnames.append(header)
        for row in matrix_io:
            split_row         = row.strip().split(delim)
            row_name, entries = split_row[0], split_row[1:]
            for column_entry in split_row:
                for column_name, entry in zip(self.colnames, entries):
                    self.matrix[column_name][row_name]=entry
            if row_name in self.rownames:
                raise Exception("Duplicate row names found in %s" % \
                                                         (matrix_io.name))
            self.rownames.append(row_name)

    def get_row(self, rowname):
        if rowname in self.rownames:
            row = [rowname]
            for colname in self.colnames: 
                row.append(self.matrix[colname][rowname])
        return row
    
    def get_col(self, colname):      
        if colname in self.colnames:
            
            return self.matrix[colname].values()
        else:
            return 0 # raise Exception("Colname %s not in matrix" % colname)

    def get_entry(self, colname, rowname):
        if colname in self.colnames:
            if rowname in self.rownames:
                return self.matrix[colname][rowname]
            else:
                return 0 # raise Exception("Rowname %s not in matrix" % rowname)
        else:
            return 0 # raise Exception("Colname %s not in matrix" % colname)
    
    def filter_by_cols(self, set):
        current_list = []
        for row in self.rownames:
            if tuple(self.get_row(row)[1:]) == set:
                current_list.append(row)
        return current_list
        
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

    def _parse_enrichm_annotations(self, annotations_path):
        modules              = set()
        genomes              = set()
        genome_to_annotation = {}
        annotations_io = open(annotations_path)
        annotations_io.next() # Skip header
        for line in annotations_io:
            sline = line.strip().split('\t')
            genome_id, module_id = sline[:2]
            modules.add(module_id)
            genomes.add(genome_id)
            if genome_id in genome_to_annotation:
                genome_to_annotation[genome_id].add(module_id)
            else:
                genome_to_annotation[genome_id] = set([module_id])
        return genome_to_annotation, modules, genomes
    
    def _parse_annotation_matrix(self, annotation_matrix_path):

        '''        
        Parameters
        ----------
        annotation_matrix : String. Path to file containing a matrix of genome annotations.        
        '''

        annotation_matrix = Matrix(annotation_matrix_path)
        modules = set(annotation_matrix.rownames)
        genomes = set(annotation_matrix.colnames)
        genome_to_annotation = {}

        for genome in genomes:
            genome_to_annotation[genome] = set()
            for module in modules:
                if int(annotation_matrix.get_entry(genome, module))>0:
                    genome_to_annotation[genome].add(module)
        
        return genome_to_annotation, modules, genomes
    

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

    def calculate_portions(self, modules, combination_dict, annotations_dict, genome_list):
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
        '''
        ### ~ TODO: Add a portion of genome column too?

        output_lines = [['Module'] + combination_dict.keys()]
        
        for module in modules:
            output_line = [module]
            for group_name, genome_list in combination_dict.items():
                if len(genome_list)>0:
                    coverage = len([genome for genome in genome_list 
                                    if module in annotations_dict[genome]])
                    all      = float(len(genome_list))
                    output_line.append(str(coverage/all))
                else:
                    output_line.append('0.0')
            output_lines.append(output_line)

        return output_lines

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
                string = '\t'.join(line) + '\n'
                out_io.write(string)


    def do(self, enrichm_annotations, annotation_matrix, metadata,
           subset_modules, abundances, no_ivi, no_gvg, no_ivg, cutoff, 
           threshold, multi_test_correction, output_directory):
        '''
        
        Parameters
        ----------
        
        Output
        ------
        '''
        logging.info('Parsing input annotations')
        if enrichm_annotations:
            annotations_dict, modules, genomes \
                        = self._parse_enrichm_annotations(enrichm_annotations)
        elif annotation_matrix:
            annotations_dict, modules, genomes \
                        = self._parse_annotation_matrix(annotation_matrix)
        else:
            raise Exception("Programming error")      
        if subset_modules:
            modules = subset_modules

        logging.info("Comparing sets of genomes")
        metadata = Matrix(metadata)
        metadata_value_lists = \
            [set(metadata.get_col(col)) for col in metadata.colnames]
        combination_dict = {}

        for combination in product(*metadata_value_lists):
            genome_list = metadata.filter_by_cols(combination)
            combination_dict['_'.join(combination)]=genome_list
        
        t = Test(annotations_dict,
                 modules,
                 genomes ,
                 combination_dict,
                 self.check_annotation_type(modules),
                 threshold,
                 multi_test_correction)

        if no_gvg:
            gvg = None
        else:
            gvg = list(combinations(combination_dict.keys(), 2)) # For fisher's exact
        if no_ivi:
            ivi = None
        else:
            ivi = list(combinations(chain(*combination_dict.values()), 2)) # For T-test
        if no_ivg:
            ivg = None
        else:
            ivg = {} # For Z score test
            for group_name, genome_list in combination_dict.items():
                if len(genome_list)>1:
                    ivg[group_name] = []
                    for idx, genome in enumerate(genome_list):
                        other=genome_list[:idx] + genome_list[idx+1:]
                        ivg[group_name].append((genome, other))

        test_results = t.do(ivi, ivg, gvg)

        for test_result_lines, test_result_output_file in test_results:
            test_result_output_path = os.path.join(output_directory,
                                                   test_result_output_file)
            self._write(test_result_lines, test_result_output_path)

        proportions     = self.calculate_portions(modules, combination_dict, annotations_dict, genome_list)
        portions_path   = os.path.join(output_directory, self.PROPORTIONS)
        self._write(proportions, portions_path)

        ###############################################################################
        ###############################################################################
        # code below needs work...

        if abundances:
            logging.info("Generating enrichment matrix")
            abundances = Matrix(abundances)
            
            output_sample_enrichment_matrix = \
                                        output_directory + self.SAMPLE_MATRIX_SUFFIX
            output_lines = ['\t'.join(['Module'] + abundances.colnames)]
            for module in modules:
                output_line = [module]
                for sample in abundances.colnames:
                    
                    module_prevalence = 0.0
                    for genome in genomes:
                        if module in annotations_dict[genome]:
                            genome_abundance = abundances.get_entry(sample, genome)
                            module_prevalence+=float(genome_abundance)
                    output_line.append(str(module_prevalence))
                output_lines.append('\t'.join(output_line))
            
            logging.info("Writing results to file: %s" \
                            % output_sample_enrichment_matrix)
            with open(output_sample_enrichment_matrix, 'w') as output_enrichment_matrix_io:
                output_enrichment_matrix_io.write('\n'.join(output_lines))
            logging.info('Done!')
            
            output_enrichment_matrix = output_directory + self.MATRIX_SUFFIX
            output_lines = ['\t'.join(['Module', 'Sample', 
                                       'Group', 'Abundance'])]
            for module in modules:
                for sample in abundances.colnames:
                    for group in product(*metadata_value_lists):
                        group_name='_'.join(group)
                        genomes_in_group_list = metadata.filter_by_cols(group)                        
                        module_prevalence = 0.0
                        for genome in genomes_in_group_list:
                            if module in annotations_dict[genome]:
                                genome_abundance = abundances.get_entry(sample, genome)
                                module_prevalence+=float(genome_abundance)
                        output_line = [module, sample, 
                                       group_name, str(module_prevalence)]
                        
                        output_lines.append('\t'.join(output_line))
                
            logging.info("Writing results to file: %s" % output_enrichment_matrix)
            with open(output_enrichment_matrix, 'w') as output_enrichment_matrix_io:
                output_enrichment_matrix_io.write('\n'.join(output_lines))
            logging.info('Done!')

class Test(Enrichment):

    def __init__(self, genome_annotations, modules, genomes, groups, annotation_type, threshold, multi_test_correction):

        d=Databases()
        
        self.IVI_OUTPUT             = 'ivi_results.tsv'
        self.IVG_OUTPUT             = 'ivg_results.tsv'
        self.GVG_OUTPUT             = 'gvg_results.tsv'

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
        
        if annotation_type==Enrichment.PFAM:
            self.genome_annotations = {}
            for key, item in genome_annotations.items():
                self.genome_annotations[key] = set([x.split('.')[0] for x in item])
        else:
            self.genome_annotations = genome_annotations

        self.pfam                   = {}
        self.id2pfam                = {}
        
        for line in open(d.PFAM_CLAN_DB):
            sline       = line.strip().split()
            pfam        = sline[0]
            clan        = sline[1]
            id          = sline[2]
            description = "%s; %s" % (id, ' '.join(sline[2:]))
            self.pfam[pfam]  = description
            self.id2pfam[id] = pfam

    def correct_multi_test(self, pvalues):

        mtc_dict = {'fdr_bh':'Benjamini/Hochberg (non-negative)'}

        logging.info('Applying multi-test correction using the %s method' % (mtc_dict[self.multi_test_correction]) )

        corrected_pvals \
            = sm.multipletests(pvalues,
                               alpha        = self.threshold,
                               method       = self.multi_test_correction,
                               returnsorted = False,
                               is_sorted    = False)[1]        
        return corrected_pvals

    def _strip_kegg_definitions(self, definition):
        kos_list = [ko for ko in re.split("[^\w]", definition) if ko]
        return kos_list

    def gather_genome_annotations(self, group, target_kos):
        genome_ko_list = []
        for genome in self.groups[group]:
            genome_target_list = set(self.genome_annotations[genome])\
                                             .intersection(target_kos)
            genome_ko_list.append(len(genome_target_list))
        return genome_ko_list

    def fisher(self, ivi):
        '''
        Parameters
        ----------
        
        Output
        ------
        '''

        output_lines = [['module', 'genome_1', 'genome_2', 'count', 'stat', 'p_value', 'corrected_p_value', 'description']]
        
        pvals=[]
        
        if self.annotation_type==Enrichment.TIGRFAM:
            logging.warning('Comparisons are not possible with TIGRFAM because heirarchical classifications are needed (like in KEGG, COG or PFAM).')
        else:
            for genome_1, genome_2 in ivi:

                if self.annotation_type==Enrichment.KEGG:
                    iterator = self.m2def
                    annotation_description = self.m

                elif self.annotation_type==Enrichment.PFAM:
                    iterator = self.clan2pfam
                    annotation_description = self.clan_to_description

                for module, module_definition in iterator.items():
                    module_kos    = self._strip_kegg_definitions(module_definition)
                    
                    genome_1_comp = set(module_kos).intersection(self.genome_annotations[genome_1])
                    genome_2_comp = set(module_kos).intersection(self.genome_annotations[genome_2])

                    row_1 = [len(genome_1_comp), 
                             len(genome_2_comp)]
                    row_2 = [(len(self.genome_annotations[genome_1]) - len(genome_1_comp)),
                             (len(self.genome_annotations[genome_2]) - len(genome_2_comp))]
                    stat, p_value = \
                            scipy.stats.fisher_exact([row_1,
                                                      row_2])
                    pvals.append(p_value)
                    output_lines.append([module,
                                         genome_1,
                                         genome_2,
                                         str(len(module_kos)),
                                         str(stat),
                                         str(p_value),
                                         annotation_description[module]])

        corrected_pvals = self.correct_multi_test(pvals)
        output_lines = [x[:len(x)-1]+[str(y)]+x[len(x)-1:] for x, y in zip(output_lines, list(corrected_pvals))]

        return output_lines
        
    def ttest(self, gvg):
        output_lines = [['module', 'group_1', 'group_2', 'count', 't_stat', 'p_value','corrected_p_value', 'description']]
        pvals        = []

        if self.annotation_type==Enrichment.TIGRFAM:
            logging.warning('Comparisons are not possible with TIGRFAM because heirarchical classifications are needed (like in KEGG, COG or PFAM).')
        else:
            for group_1, group_2 in gvg:
                
                if self.annotation_type == Enrichment.KEGG:
                    annotation_dict = self.m2def.items()
                    annotation_description = self.m

                elif self.annotation_type == Enrichment.PFAM:
                    annotation_dict = self.clan2pfam.items()
                    annotation_description = self.clan_to_description
                else:
                    raise Exception("Programming error")

                for module, definition in annotation_dict:

                    module_list = self._strip_kegg_definitions(definition)
                    
                    group_1_module_kos = self.gather_genome_annotations(group_1, module_list)
                    group_2_module_kos = self.gather_genome_annotations(group_2, module_list)

                    if(len(group_1_module_kos)>1 and len(group_2_module_kos)>1):
                        if(sum(group_1_module_kos)>0 and sum(group_2_module_kos)>0):
                            t_stat, p_value = \
                                scipy.stats.ttest_ind(np.array(group_1_module_kos),
                                                      np.array(group_2_module_kos),
                                                      equal_var  = False,
                                                      nan_policy = 'raise')
                            pvals.append(p_value)
                            output_lines.append([module,
                                                 group_1,
                                                 group_2,
                                                 str(len(module_list)),
                                                 str(t_stat),
                                                 str(p_value),
                                                 annotation_description[module]])
        corrected_pvals = self.correct_multi_test(pvals)
        output_lines = [x[:len(x)-1]+[str(y)]+x[len(x)-1:] for x, y in zip(output_lines, list(corrected_pvals))]

        return output_lines

    def zscore(self, ivg):
        output_lines = [['module', 'group', 'genome', 'count', 'z_score', 'p_value','corrected_p_value', 'description']]
        pvals        = []

        if self.annotation_type==Enrichment.TIGRFAM:
            logging.warning('Comparisons are not possible with TIGRFAM because heirarchical classifications are needed (like in KEGG, COG or PFAM).')
        else:
            for group, comparisons in ivg.items():
                for genome, reference_group in comparisons:
                    

                    if self.annotation_type == Enrichment.KEGG:
                        annotation_dict = self.m2def.items()
                        annotation_description = self.m
                    elif self.annotation_type == Enrichment.PFAM:
                        annotation_dict = self.clan2pfam.items()
                        annotation_description = self.clan_to_description
                    else:
                        raise Exception("Programming error")

                    for module, definition in annotation_dict:
                        
                        module_list = self._strip_kegg_definitions(definition)
                        genome_comp \
                            = len(set(self.genome_annotations[genome]).intersection(module_list))
                        reference_group_comp_list \
                            = np.array([len(set(self.genome_annotations[reference_genome]).intersection(module_list))
                                        for reference_genome in reference_group])
                        reference_group_comp_sd \
                            = np.std(reference_group_comp_list, axis=0)
                        reference_group_comp_mean \
                            = np.mean(reference_group_comp_list, axis=0)

                        z_score = (genome_comp-reference_group_comp_mean) / reference_group_comp_sd
                        p_value = 2-2*scipy.stats.norm.cdf(z_score)
                        pvals.append(p_value)
                        output_lines.append([module, 
                                              group, 
                                              genome, 
                                              str(len(module_list)), 
                                              str(z_score), 
                                              str(p_value), 
                                              annotation_description[module]])

        corrected_pvals = self.correct_multi_test(pvals)
        output_lines = [x[:len(x)-1]+[str(y)]+x[len(x)-1:] for x, y in zip(output_lines, list(corrected_pvals))]

        return output_lines

    def do(self, ivi, ivg, gvg):
        
        logging.info('Running genome comparisons')
        results=[]

        if ivi:
            logging.info('Running individual vs individual comparisons')
            results.append( (self.fisher(ivi), self.IVI_OUTPUT) )
        if ivg:
            logging.info('Running individual vs group comparisons')
            results.append( (self.zscore(ivg), self.IVG_OUTPUT) )
        if gvg:
            logging.info('Running group vs group comparisons')
            results.append( (self.ttest(gvg), self.GVG_OUTPUT) )
            
        return results
