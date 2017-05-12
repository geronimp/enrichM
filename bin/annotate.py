#!/usr/bin/env python
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
# System imports
import logging
import subprocess
import os 


from databases import KO_DB, PFAM_DB, TIGRFAM_DB
from matrix_generator import MatrixGenerator
from gff_generator import GffGenerator
from genome import Genome
###############################################################################
################################ - Classes - ##################################

class Annotate:

    GENOME_BIN = 'genome_bin'
    GENOME_PROTEINS = 'genome_proteins'    
    GENOME_KO = 'annotations_ko'
    GENOME_COG = 'annotations_cog'
    GENOME_PFAM = 'annotations_pfam'
    GENOME_TIGRFAM = 'annotations_tigrfam'
    OUTPUT_PFAM = 'pfam.tsv'

    PROTEINS_SUFFIX = '.faa'
    OUTPUT_SUFFIX = '.tsv'
    
    def __init__(self, genome_files, output_directory, ko, pfam, 
                 tigrfam, cog, evalue, bit, id, aln_query, aln_reference, 
                 threads):
        # Define inputs and outputs
        self.genome_file_list = genome_files         
        self.output_directory = output_directory

        # Define type of annotation to be carried out
        self.ko               = ko 
        self.pfam             = pfam 
        self.tigrfam          = tigrfam 
        self.cog              = cog 
        # Cutoffs
        self.evalue           = evalue 
        self.bit              = bit 
        self.id               = id 
        self.aln_query        = aln_query
        self.aln_reference    = aln_reference
        # Parameters
        self.threads          = threads

    def prep(self):
        '''
        Do any preparation specific to the genome annotation pipeline. 

        Outputs
        -------
        returns the directory with all genome ids sym-linked into it.
        '''
        # link all the genomes into one file
        
        genome_directory = os.path.join(self.output_directory, 
                                       self.GENOME_BIN)
        os.mkdir(os.path.join(self.output_directory, self.GENOME_BIN))
        for genome_path in self.genome_file_list:
            os.symlink(os.path.join(os.getcwd(),
                                    genome_path), 
                       os.path.join(genome_directory, 
                                    os.path.basename(genome_path)
                                    )
                       )    
        return genome_directory

    def call_proteins(self, genome_directory):
        '''
        Use prodigal (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119) to call proteins within the genomes 
        
        Parameters
        ----------
        genome_directory  - string. Directory containing .fna files for each  
                            input genome

        Outputs
        -------
        returns the directory containing an .faa file for each input genomes
        '''   
        output_directory_path = os.path.join(self.output_directory, 
                                             self.GENOME_PROTEINS) 
        os.mkdir(output_directory_path)
        genomes_list = []

        for genome in os.listdir(genome_directory):
            
            genome_path = os.path.join(genome_directory,
                                       genome)
            output_genome = os.path.splitext(genome)[0] + self.PROTEINS_SUFFIX 
            output_genome_path = os.path.join(output_directory_path,
                                              output_genome)
            logging.info("    - Calling proteins for genome: %s" % (genome))
            cmd = 'prodigal -q -p meta -o /dev/null -a %s -i %s' % (output_genome_path,
                                                    genome_path)
            logging.debug(cmd)
            subprocess.call(cmd, shell = True)
            
            genomes_list.append(Genome(output_genome_path))
        
        return genomes_list
    
    def annotate_ko(self, genome_faa_directory):
        '''
        Annotate the proteins encoded by each genome with KO ids using either BLAST or using HMM searches (no implemented yet).

        Parameters
        ----------        
        genome_faa_directory  - string. Directory containing .faa files for 
                                each input genome
        
        Outputs
        -------
        returns a directory containing the search results for each of the input population genomes, and a frequency matrix contining with the KOs as rows, and the genomes as columns.
        '''        
        output_directory_path = os.path.join(self.output_directory, 
                                             self.GENOME_KO)
        os.mkdir(output_directory_path)
        self._diamond_search(KO_DB, genome_faa_directory, output_directory_path)
        results_files = [os.path.join(output_directory_path, genome_result) 
                        for genome_result in os.listdir(output_directory_path)]
        return results_files

    def annotate_cog(self): 
        pass
        # TODO haven't prepared the reference database yet

    def _diamond_search(self, database, genome_faa_directory, output_directory_path):
        '''
        Carry out a diamond searches on each genome within a given directory. 

        Parameters
        ----------   
        databases             - string. Path to .dmnd database to use for searching     
        genome_faa_directory  - string. Directory containing .faa files for 
                                each input genome
        output_directory_path - string. Path to directory to output results into        
        '''  
        for genome in os.listdir(genome_faa_directory):
            genome_path = os.path.join(genome_faa_directory,
                                       genome)
            output_genome = os.path.splitext(genome)[0] + self.OUTPUT_SUFFIX 
            output_genome_path = os.path.join(output_directory_path,
                                              output_genome)            
            cmd = 'diamond blastp --outfmt 6 --max-target-seqs 1 --query %s --out %s --db %s --threads %s ' % (genome_path,  output_genome_path, database, self.threads)
            if self.evalue:
                cmd += '--evalue %f ' % (self.evalue) 
            if self.bit:
                cmd += '--min-score %f ' % (self.bit)
            if self.id:
                cmd += '--id %f ' % (self.id)
            logging.debug(cmd)
            subprocess.call(cmd, shell = True)

    def annotate_pfam(self, genomes_list):
        '''
        Annotate the proteins encoded by each genome with pfam ids using HMM searches.

        Parameters
        ----------        
        genomes_list  - object.        

        Outputs
        -------
        returns a directory containing the search results for each of the input population genomes,
        and a frequency matrix contining with the pfam ids as rows, and the genomes as columns.
        '''    
        output_directory_path = os.path.join(self.output_directory, 
                                             self.GENOME_PFAM)
        os.mkdir(output_directory_path)
        for genome in genomes_list:
            self._hmm_search(PFAM_DB, genome.path, output_directory_path)
            results_files = [os.path.join(output_directory_path, genome_result) 
                        for genome_result in os.listdir(output_directory_path)]
        return results_files

    def annotate_tigrfam(self, genome_faa_directory):
        '''
        Annotate the proteins encoded by each genome with tigrfam ids using HMM searches.

        Parameters
        ----------        
        genome_faa_directory  - string. Directory containing .faa files for 
                                each input genome
        
        Outputs
        -------
        returns a directory containing the search results for each of the input population genomes,
        and a frequency matrix contining with the tigrfam ids as rows, and the genomes as columns.
        '''    
        output_directory_path = os.path.join(self.output_directory, 
                                             self.GENOME_TIGRFAM)
        os.mkdir(output_directory_path)      
        self._hmm_search(TIGRFAM_DB, genome_faa_directory, output_directory_path)

        results_files = [os.path.join(output_directory_path, genome_result) 
                        for genome_result in os.listdir(output_directory_path)]
        return results_files

    def _hmm_search(self, database, genome_path, output_directory_path):
        '''
        Carry out a hmmsearches on each genome within a given directory. 

        Parameters
        ----------   
        databases             - string. Path to HMM to use for searching     
        genome_faa_directory  - string. Directory containing .faa files for 
                                each input genome
        output_directory_path - string. Path to directory to output results into        
        '''  
        output_genome = os.path.basename(os.path.splitext(genome_path)[0]) + self.OUTPUT_SUFFIX 
        output_genome_path = os.path.join(output_directory_path, output_genome)            
        cmd = "hmmsearch --cpu %s -o /dev/null --noali --domtblout %s " \
                          % (self.threads,output_genome_path)
        if self.evalue:
            cmd += '-E %f ' % (self.evalue) 
        if self.bit:
            cmd += '-T %f ' % (self.bit)
        if self.id:
            logging.warning("--id flag not used for hmmsearch")
        cmd += "%s %s " % (database, genome_path)
        logging.debug(cmd)
        subprocess.call(cmd, shell = True)        

    def do(self, genome_directory, proteins_directory):
        '''
        Run Annotate pipeline for enrichM

        Parameters
        ----------
        genome_directory  - string. Path to directory containing genomes
        proteins_directory  - string. Path to directory containing proteins (.faa files) for genomes
        '''
        logging.info("Running pipeline: annotate")
        results = {}
        logging.info("Setting up for genome annotation")
        if genome_directory:
            self.genome_directory = genome_directory
        else:
            self.genome_directory = self.prep()

        if proteins_directory:
            pass
            # TODO: Implement me
            # genomes_dictionary = self.parse_proteins(proteins_directory)
        else:
            logging.info("Calling proteins for annotation")
            genomes_list = self.call_proteins(self.genome_directory)

        logging.info("Starting annotation:")
        if self.ko:
            logging.info('    - Annotating genomes with ko ids')
            ko_result_paths = self.annotate_ko(self.genome_faa_directory)
            freq_table = os.path.join(self.output_directory, self.OUTPUT_PFAM)
            mg.from_blast_results(ko_result_paths)

        if self.cog:
            logging.info('    - Annotating genomes with cog ids')
            logging.info('    - cog annotation currently not implemented')

        if self.pfam:
            logging.info('    - Annotating genomes with pfam ids')
            pfam_result_paths = self.annotate_pfam(genomes_list)
            logging.info('    - Generating pfam frequency table')
            mg = MatrixGenerator(MatrixGenerator.PFAM)
            freq_table = os.path.join(self.output_directory, self.OUTPUT_PFAM)
            pfam_annotations = mg.from_hmmsearch_results(pfam_result_paths, 
                                                         freq_table, 
                                                         self.evalue, 
                                                         self.bit, 
                                                         self.aln_query, 
                                                         self.aln_reference)

        if self.tigrfam:
            logging.info('    - Annotating genomes with tigrfam ids')
            freq_table = os.path.join(self.output_directory, self.OUTPUT_PFAM)
            tigrfam_result_paths = self.annotate_tigrfam(self.genome_faa_directory)
            mg.from_hmmsearch_results(tigrfam_result_paths)
        # TODO: Generate GFF files
        # TODO: Generate frequency matrix files

