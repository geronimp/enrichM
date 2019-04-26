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

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.7"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

###############################################################################
# System imports
from enrichm.genome import Genome, AnnotationParser
from enrichm.gff_generator import GffGenerator
from enrichm.matrix_generator import MatrixGenerator
from enrichm.databases import Databases
from enrichm.sequence_io import SequenceIO
import numpy as np
import statsmodels.sandbox.stats.multicomp as sm
import multiprocessing as mp
import pickle
import shutil
import tempfile
import logging
import subprocess
import os

# Local

###############################################################################
###############################################################################

def parse_genomes(params):
    genome = Genome(*params)
    return genome

###############################################################################
################################ - Classes - ##################################

class Annotate:

    GENOME_BIN                      = 'genome_bin'
    GENOME_PROTEINS                 = 'genome_proteins'    
    GENOME_KO                       = 'annotations_ko'
    GENOME_KO_HMM                   = 'annotations_ko_hmm'
    GENOME_EC                       = 'annotations_ec'
    GENOME_PFAM                     = 'annotations_pfam'
    GENOME_TIGRFAM                  = 'annotations_tigrfam'
    GENOME_HYPOTHETICAL             = 'annotations_hypothetical'
    GENOME_CAZY                     = 'annotations_cazy'
    GENOME_GFF                      = 'annotations_gff'
    GENOME_OBJ                      = 'annotations_genomes'
    OUTPUT_KO                       = 'ko_frequency_table.tsv'
    OUTPUT_KO_HMM                   = 'ko_hmm_frequency_table.tsv'
    OUTPUT_EC                       = 'ec_frequency_table.tsv'
    OUTPUT_PFAM                     = 'pfam_frequency_table.tsv'
    OUTPUT_TIGRFAM                  = 'tigrfam_frequency_table.tsv'
    OUTPUT_CAZY                     = 'cazy_frequency_table.tsv'
    OUTPUT_HYPOTHETICAL_CLUSTER     = 'cluster_frequency_table.tsv'
    OUTPUT_HYPOTHETICAL_ORTHOLOG    = 'ortholog_frequency_table.tsv'
    OUTPUT_HYPOTHETICAL_ANNOTATIONS = 'hypothetical_annotations.tsv'
    OUTPUT_DIAMOND                  = "DIAMOND_search"
    GFF_SUFFIX                      = '.gff'
    PROTEINS_SUFFIX                 = '.faa'
    ANNOTATION_SUFFIX               = '.tsv'
    PICKLE_SUFFIX                   = '.pickle'
    
    def __init__(self,
                 output_directory,
                 ko, ko_hmm, pfam, tigrfam, cluster, ortholog, cazy, ec,
                 evalue, bit, id, aln_query, aln_reference, c, cut_ga, cut_nc, cut_tc, cut_hmm, inflation, chunk_number, chunk_max, count_domains,
                 threads, parallel, suffix, light):

        # Define inputs and outputs
        self.output_directory = output_directory

        # Define type of annotation to be carried out
        self.ko               = ko 
        self.ko_hmm           = ko_hmm
        self.pfam             = pfam 
        self.tigrfam          = tigrfam 
        self.cluster          = cluster 
        self.ortholog         = ortholog

        self.cazy             = cazy
        self.ec               = ec
        
        # Cutoffs
        self.evalue           = evalue 
        self.bit              = bit 
        self.id               = id 
        self.aln_query        = aln_query
        self.aln_reference    = aln_reference
        self.c                = c
        self.cut_ga           = cut_ga
        self.cut_nc           = cut_nc
        self.cut_tc           = cut_tc
        self.cut_hmm          = cut_hmm
        self.inflation        = inflation
        self.chunk_number     = chunk_number
        self.chunk_max        = chunk_max
        self.count_domains    = count_domains

        # Parameters
        self.threads          = threads
        self.parallel         = parallel
        self.suffix           = suffix
        self.light            = light

        # Set up multiprocesses pool
        self.pool             = mp.Pool(processes = int(self.parallel))

        # Load databases
        self.databases        = Databases()

    def prep_genome(self, genome_file_list, genome_directory):
        '''
        Do any preparation specific to the genome annotation pipeline. 
    
        Inputs
        ------
        genome_file_list - List. list of strings, each a path to a file
        containing a genome 

        Outputs
        -------
        returns the directory with all genome ids sym-linked into it.
        '''
        # link all the genomes into one file
        logging.info('Preparing genomes for annotation')
        if genome_file_list:
            os.mkdir(genome_directory)
            genome_paths = list()
            for genome_path in genome_file_list:
                if genome_path.endswith(self.suffix):
                    genome_paths.append(genome_path) 
            cmd = "xargs --arg-file=/dev/stdin cp --target-directory=%s" % genome_directory
            logging.debug(cmd)
            process = subprocess.Popen(["bash", "-c", cmd], 
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       universal_newlines=True)
            process.communicate(input=str('\n'.join(genome_paths)))
        return genome_directory

    def call_proteins(self, genome_directory):
        '''
        Use prodigal (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
        to call proteins within the genomes 
        
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
        genome_list     = list()
        genome_paths    = list()

        for genome in os.listdir(genome_directory):
            if genome.endswith(self.suffix):
                genome_paths.append(os.path.splitext(genome)[0])

        logging.info("    - Calling proteins for %i genomes" % (len(genome_paths)))
        cmd = "ls %s/*%s | sed 's/%s//g' | grep -o '[^/]*$' | parallel -j %s prodigal -q -p meta -o /dev/null -a %s/{}%s -i %s/{}%s  > /dev/null 2>&1" \
                % (genome_directory,
                   self.suffix,
                   self.suffix,
                   self.parallel,
                   output_directory_path,
                   self.PROTEINS_SUFFIX,
                   genome_directory,
                   self.suffix)

        logging.debug(cmd)
        subprocess.call(cmd, shell = True)
        logging.debug('Finished')
        for genome_protein, genome_nucl in zip(os.listdir(output_directory_path), os.listdir(genome_directory)):
            output_genome_protein_path = os.path.join(output_directory_path, genome_protein)
            output_genome_nucl_path = os.path.join(genome_directory, genome_nucl)
            genome = (self.light, output_genome_protein_path, output_genome_nucl_path)
            genome_list.append(genome)
        
        return genome_list
    
    def annotate_diamond(self, genomes_list, database, parser_type, ids_type, output_subdirectory):
        '''
        Annotate the proteins encoded by each genome with KO ids using either BLAST or using HMM
        searches (no implemented yet).

        Parameters
        ----------        
        genome_faa_directory  - string. Directory containing .faa files for 
                                each input genome
        
        Outputs
        -------
        returns a directory containing the search results for each of the input population genomes, 
        and a frequency matrix contining with the KOs as rows, and the genomes as columns.
        '''        

        output_directory_path = os.path.join(self.output_directory, 
                                             output_subdirectory)
        genome_dict = {genome.name:genome for genome in genomes_list}
        os.mkdir(output_directory_path)
        specific_cutoffs = None
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(str.encode('\n'.join(["sed \"s/>/>%s~/g\" %s" % (genome.name, genome.path) for genome in genomes_list])))
            temp.flush()
            output_annotation_path = os.path.join(output_directory_path, self.OUTPUT_DIAMOND) + self.ANNOTATION_SUFFIX
            logging.info('    - BLASTing genomes' )
            self._diamond_search(temp.name, output_annotation_path, database)

            for genome_name, batch in self.get_batches(output_annotation_path):
                
                if batch:
                    genome = genome_dict[genome_name]
                    genome.add(batch, 
                                 self.evalue, 
                                 self.bit, 
                                 self.aln_query, 
                                 self.aln_reference,
                                 specific_cutoffs,
                                 parser_type,
                               ids_type)
    
    def get_batches(self, input_file):
        last = None
        input_file_io = open(input_file)

        for line in input_file_io:
            split_line = line.strip().split('\t')
            genome_id, _ = split_line[0].split('~')
            
            if last is None:
                last = genome_id
                batch = [split_line]
            
            else:
            
                if last==genome_id:
                    batch.append(split_line)
            
                else:
            
                    yield last, batch
                    batch = [split_line]
                    last = genome_id
        
        if last is None:
        
            yield None, None
        
        else:

            yield last, batch

    def _diamond_search(self, tmp_name, output_path, database):
        '''
        Carry out a diamond blastp search. 

        Parameters
        ----------   
        input_genome_path     - string. Path to file containing .faa file for 
                                an input genome
        output_path           - string. Path to file to output results into    
        databases             - string. Path to HMM to use for searching         
        '''  
         
        cmd = 'bash %s | diamond blastp --quiet --outfmt 6 --max-target-seqs 1 --query /dev/stdin --out %s --db %s --threads %s ' \
                            % (tmp_name, output_path, database, self.threads)
        
        if self.evalue:
            cmd += '--evalue %s ' % (str(self.evalue)) 
        
        if self.bit:
            cmd += '--min-score %s ' % (str(self.bit))
        
        if self.id:
            cmd += '--id %s ' % (str(self.id*100))
        
        if self.aln_query:
            cmd += "--query-cover %s " % (str(self.aln_query * 100))
        
        if self.aln_reference:
            cmd += "--subject-cover %s " % (str(self.aln_reference * 100))

        logging.debug(cmd)
        subprocess.call(cmd, shell = True)
        logging.debug('Finished')

    def hmmsearch_annotation(self, genomes_list, output_directory_path, database, ids_type, parser):
        '''
        Annotate the proteins encoded by each genome with pfam ids using HMM searches.

        Parameters
        ----------        
        genomes_list  - list. list of Genome objects        

        '''    
        os.mkdir(output_directory_path)
        genome_dict = {genome.name: genome for genome in genomes_list}

        if (ids_type == AnnotationParser.TIGRFAM or 
            ids_type == AnnotationParser.PFAM):
            hmmcutoff=True
        else:
            hmmcutoff=False
        
        if ids_type == AnnotationParser.KO_HMM:
            specific_cutoffs = self.databases.parse_ko_cutoffs()
        else:
            specific_cutoffs = None

        self._hmm_search(output_directory_path, database, hmmcutoff)

        for genome_annotation in os.listdir(output_directory_path):
            genome_id = os.path.splitext(genome_annotation)[0]
            genome = genome_dict[genome_id]
            output_annotation_path = os.path.join(output_directory_path, genome_annotation)
            genome.add(output_annotation_path, 
                         self.evalue, 
                         self.bit, 
                         self.aln_query, 
                         self.aln_reference,
                         specific_cutoffs,
                         parser,
                         ids_type)

    def annotate_hypothetical(self, genomes_list, ortholog):
        '''
        Sort proteins coded by each genome into homologous clusters.  
        
        Inputs
        ------
        genomes_list  - list. list of Genome objects        

        '''
        output_directory_path = os.path.join(self.output_directory, 
                                             self.GENOME_HYPOTHETICAL)
        os.mkdir(output_directory_path)      

        with tempfile.NamedTemporaryFile() as temp:

            temp.write(str.encode('\n'.join(["sed \"s/>/>%s~/g\" %s" % (genome.name, genome.path) for genome in genomes_list])))
            temp.flush()  

            tmp_dir = tempfile.mkdtemp()  
                
            db_path = os.path.join(output_directory_path, "db")
            clu_path = os.path.join(output_directory_path, "clu")
            align_path = os.path.join(output_directory_path, "alignDb")
            blast_output_path = os.path.join(output_directory_path, "alignDb.m8")
            formatted_blast_output_path = os.path.join(output_directory_path, "alignDb.formatted.m8")

            clu_tsv_path = os.path.join(output_directory_path, "hypothetical_clusters.tsv")
            logging.info('    - Generating MMSeqs2 database')
            cmd = "bash %s | sponge | mmseqs createdb /dev/stdin %s -v 0 > /dev/null 2>&1" % (
                temp.name, db_path)
            logging.debug(cmd)
            subprocess.call(cmd, shell = True)
            logging.debug('Finished')
            logging.info('    - Clustering genome proteins')
            cmd = "mmseqs cluster %s %s %s --max-seqs 1000 --threads %s --min-seq-id %s -e %f -c %s -v 0 " % (
                db_path, clu_path, tmp_dir, self.threads, self.id, self.evalue, self.c)
            logging.debug(cmd)
            subprocess.call(cmd, shell = True)
            logging.debug('Finished')
            logging.info('    - Extracting clusters')

            cmd = 'mmseqs createtsv %s %s %s %s -v 0 > /dev/null 2>&1' % (
                db_path, db_path, clu_path, clu_tsv_path)
            logging.debug(cmd)
            subprocess.call(cmd, shell = True)
            logging.debug('Finished')
            
            logging.info('    - Computing Smith-Waterman alignments for clustering results')
            cmd = "mmseqs alignall %s %s %s --alignment-mode 3 -v 0  " % (
                db_path, clu_path, align_path)
            logging.debug(cmd)
            subprocess.call(cmd, shell = True)
            logging.debug('Finished')

            logging.info('    - Converting to BLAST-like output')
            cmd = "mmseqs createtsv %s %s %s %s -v 0 > /dev/null 2>&1   " % (
                db_path, db_path, align_path, blast_output_path)
            # --format-output query,target,bits
            logging.debug(cmd)
            subprocess.call(cmd, shell = True)
            logging.debug('Finished')
            logging.info('    - Reformatting BLAST output')
            cmd = "OFS=\"\t\" awk 'FNR==NR{a[$1]=$2;next}{$3=a[$3]; $1=\"\"; for(i=2;i<NF;i++){printf(\"%s\t\",$i)} printf(\"\\n\")}' %s %s | cut -f1,2,5 > %s" % ("%s", db_path + '.lookup', blast_output_path, formatted_blast_output_path)
            logging.debug(cmd)
            subprocess.call(cmd, shell=True)
            logging.debug('Finished')

        ortholog_dict = self.run_mcl(formatted_blast_output_path,
                                     output_directory_path)
        ortholog_ids = ortholog_dict.keys()
        cluster_ids = self.parse_cluster_results(clu_tsv_path,
                                        genomes_list,
                                                ortholog_dict,
                                        output_directory_path)
        return cluster_ids, ortholog_ids

    def run_mcl(self, blast_abc, output_directory_path):
        dict_path = os.path.join(output_directory_path, "alignDb.dict")
        mci_path = os.path.join(output_directory_path, "alignDb.mci")
        cluster_path = os.path.join(output_directory_path, "mcl_clusters.tsv")
        output_path = os.path.join(output_directory_path, "mcl_clusters.convert.tsv")

        logging.info('    - Preparing network')
        ortholog_dict = dict()
        cmd = "mcxload -abc %s -write-tab %s -o %s --stream-mirror --stream-neg-log10 > /dev/null 2>&1" % (
            blast_abc, dict_path, mci_path)
        logging.debug(cmd)
        subprocess.call(cmd, shell=True)
        logging.debug('Finished')

        logging.info('    - Finding orthologs')
        ortholog_dict = dict()
        cmd = 'mcl %s -te %s -I %s -o %s  > /dev/null 2>&1' % (
            mci_path, self.threads, self.inflation, cluster_path)
        logging.debug(cmd)
        subprocess.call(cmd, shell = True)
        logging.debug('Finished')
        
        logging.info('    - Reformatting output')
        ortholog_dict = dict()
        cmd = 'mcxdump -icl %s -o %s -tabr %s > /dev/null 2>&1' % (cluster_path,  output_path, dict_path)
        logging.debug(cmd)
        subprocess.call(cmd, shell=True)
        logging.debug('Finished')

        ortholog = 1
        for line in open(output_path):
            ortholog_idx = "ortholog_%i" % ortholog
            ortholog_dict[ortholog_idx] = set()
        
            for protein in line.strip().split('\t'):
                ortholog_dict[ortholog_idx].add(protein)
        
            ortholog += 1
        
        return ortholog_dict

    def parse_cluster_results(self, 
                              cluster_output_path,
                              genomes_list,
                              ortholog_dict,
                              output_directory_path):
        '''
        Parse cluster output in tab format.
        
        Inputs
        ------
        from_cluster_results    - String. Path to mmseqs2 clustering output file
        
        Yields
        -------
        A cluster name, and a list of sequences in that cluster.
        
        '''
        logging.info('    - Parsing input cluster file: %s' % cluster_output_path)
        
        cluster_ids             = set()
        previous_cluster_name   = None
        counter                 = 0
        genome_dictionary       = {genome.name:genome for genome in genomes_list}

        with open(os.path.join(output_directory_path, self.OUTPUT_HYPOTHETICAL_ANNOTATIONS), 'w') as out_io:

            for line in open(cluster_output_path):
    
                cluster_id, member      = line.strip().split('\t')
                genome_id, sequence_id  = member.split('~')
                
                if cluster_id == previous_cluster_name:
                    genome_dictionary[genome_id].add_cluster(sequence_id, "cluster_%i" % counter)
                else:
                    counter += 1
                    previous_cluster_name = cluster_id 
                    cluster_ids.add("cluster_%i" % counter)
                    genome_dictionary[genome_id].add_cluster(sequence_id, "cluster_%i" % counter)
                out_io.write('\t'.join([genome_id, sequence_id, "cluster_%i" % counter]) + '\n')
        
        for ortholog, group in ortholog_dict.items():
            for member in group:
                genome, protein = member.split('~')
                genome_dictionary[genome].add_ortholog(protein, ortholog)
        
        return cluster_ids

    def _default_hmmsearch_options(self):
        cmd = ''
        if self.bit:
            cmd += '-T %s ' % (str(self.bit))    
        else:
            cmd += '-E %s ' % (str(self.evalue)) 
        return cmd

    def _hmm_search(self, output_path, database, hmmcutoff):
        '''
        Carry out a hmmsearch. 

        Parameters
        ----------   
        input_genome_path     - string. Path to file containing .faa file for 
                                an input genome
        output_path           - string. Path to file to output results into   
        databases             - string. Path to HMM to use for searching          
        '''
        
        input_genome_path = os.path.join(self.output_directory, self.GENOME_PROTEINS)
        cmd = "ls %s | sed 's/%s//g' | parallel -j %s hmmsearch --cpu %s -o /dev/null --noali --domtblout %s/{}%s " \
                          % (input_genome_path, self.PROTEINS_SUFFIX, self.parallel, 
                             self.threads, output_path, self.ANNOTATION_SUFFIX)
        if hmmcutoff:
            if(self.cut_ga or self.cut_nc or self.cut_tc):

                if self.cut_ga:
                    cmd += " --cut_ga "
                if self.cut_nc:
                    cmd += " --cut_nc "
                if self.cut_tc:
                    cmd += " --cut_tc "
            else:
                cmd += self._default_hmmsearch_options()
        else:
            cmd += self._default_hmmsearch_options()   

        cmd += "%s %s/{}.faa 2> /dev/null" % (database, input_genome_path)

        logging.debug(cmd)
        subprocess.call(cmd, shell = True)        
        logging.debug('Finished')
    
    def _generate_gff_files(self, genomes_list):
        '''
        Write GFF files for each of the genome objects in genomes_list

        Parameters
        ----------
        genomes_list - List. List of Genome objects
        '''
        output_directory_path = os.path.join(self.output_directory, 
                                             self.GENOME_GFF)
        os.mkdir(output_directory_path)      
        for genome in genomes_list:
            logging.info('    - Generating .gff file for %s' % genome.name)
            gff_output = os.path.join(output_directory_path, genome.name + self.GFF_SUFFIX)
            gg = GffGenerator()
            gg.write(genome, gff_output)

    def _rename_fasta(self, genomes_list):
        '''
        Rename the called proteins with annotation ids.
        
        Parameters
        ----------
        genomes_list - List. List of Genome objects
        '''
        seqio = SequenceIO()
        for genome in genomes_list:
            fd, fname = tempfile.mkstemp(suffix='.faa', text=True)
            with open(fname, 'w') as out_io:
                for description, sequence in seqio.each(open(genome.path)):
                    name = description.partition(' ')[0]
                    annotations = ' '.join(genome.sequences[name].all_annotations())
                    out_io.write( ">%s %s\n" % (name, annotations) )
                    out_io.write( str(sequence) + '\n' )
            os.close(fd)
            logging.debug('Moving %s to %s' % (fname, genome.path))
            shutil.move(fname, genome.path)

    def _pickle_objects(self, genomes_list):
        '''
        Store annotated genome objects as pickles.
        
        Parameters
        ----------
        genomes_list - List. List of Genome objects
        '''
        output_directory_path = os.path.join(self.output_directory,
                                             self.GENOME_OBJ)
        os.mkdir(output_directory_path)
        for genome in genomes_list:
            with open(os.path.join(output_directory_path, genome.name + self.PICKLE_SUFFIX), 'wb') as output:
                pickle.dump(genome, output)

    def list_splitter(self, input_list, chunk_number, chunk_max):
        list_size   = float( len(input_list) )
        chunk_size  = int(round( (list_size/chunk_number), 0 ))


        if chunk_size > chunk_max:
            chunk_size = chunk_max
        elif chunk_size < 1:
            chunk_size = list_size

        while list_size > 0:

            if len(input_list)<=chunk_size:
                yield input_list
                del input_list
            else:
                yield input_list[:chunk_size]
                del input_list[:chunk_size]

            try:
                list_size = len(input_list)
            except NameError:
                list_size = 0

    def parse_genome_inputs(self, genome_directory, protein_directory, genome_files, protein_files):
        '''
        Inputs
        ------
        
        Outputs
        -------
        
        '''

        prep_genomes_list   = list()  
        genomes_list        = list()

        if protein_directory:
            logging.info("Using provided proteins")
            protein_genome_list = list()
            
            for protein_file in os.listdir(protein_directory):
                protein_genome_list.append(os.path.join(protein_directory, protein_file))

            directory = self.prep_genome(protein_genome_list,
                                         os.path.join(self.output_directory,
                                                      self.GENOME_PROTEINS))
            
            for genome_proteins_file in os.listdir(directory):

                if genome_proteins_file.endswith(self.suffix):
                    genome = (self.light, os.path.join(directory, genome_proteins_file), None)
                    prep_genomes_list.append(genome)

        elif protein_files:
            logging.info("Using provided proteins")
            directory = self.prep_genome(protein_files,
                             os.path.join(self.output_directory,
                                          self.GENOME_PROTEINS))
            
            for protein_file in os.listdir(directory):
                prep_genomes_list.append((self.light, os.path.join(directory, os.path.basename(protein_file)), None))

        elif genome_directory:
            logging.info("Calling proteins for annotation")
            prep_genomes_list = self.call_proteins(genome_directory)
            directory = genome_directory
        
        elif genome_files:
            logging.info("Calling proteins for annotation")      
            directory = self.prep_genome(genome_files,
                                         os.path.join(self.output_directory,
                                                      self.GENOME_BIN))      
            prep_genomes_list = self.call_proteins(directory)
        
        
        for chunk in self.list_splitter(prep_genomes_list, self.chunk_number, self.chunk_max):
            genomes_list += self.pool.map(parse_genomes, chunk)

        return genomes_list

    def do(self, genome_directory, protein_directory, genome_files, protein_files):
        '''
        Run Annotate pipeline for enrichM

        Parameters
        ----------
        genome_directory    - String. Path to directory containing genomes
        protein_directory   - String. Path to directory containing proteins (.faa files) for genomes
        genome_files        - List. List of strings, each to a .fna genome file.
        protein_files       - List. List of strings, each to a .faa proteins file.
        '''

        logging.info("Running pipeline: annotate")
        logging.info("Setting up for genome annotation")
        genomes_list = self.parse_genome_inputs(genome_directory, protein_directory, genome_files, protein_files)
        if len(genomes_list)==0:
            logging.error('There were no genomes found with the suffix %s within the provided directory' \
                                        %  (self.suffix))
        else:
            logging.info("Starting annotation:")

            if (self.cluster or self.ortholog):
                logging.info('    - Annotating genomes with hypothetical clusters')
                cluster_ids, ortholog_ids = self.annotate_hypothetical(genomes_list, self.ortholog)

                logging.info('    - Generating hypotheticals frequency table') 
                mg = MatrixGenerator(MatrixGenerator.HYPOTHETICAL, cluster_ids)
                freq_table = os.path.join(self.output_directory, self.OUTPUT_HYPOTHETICAL_CLUSTER)
                mg.write_matrix(genomes_list, self.count_domains, freq_table)

                if self.ortholog:
                    mg = MatrixGenerator(MatrixGenerator.ORTHOLOG, ortholog_ids)
                    freq_table = os.path.join(self.output_directory, self.OUTPUT_HYPOTHETICAL_ORTHOLOG)
                    mg.write_matrix(genomes_list, self.count_domains, freq_table)
                

            if self.ko:
                annotation_type = AnnotationParser.BLASTPARSER
                logging.info('    - Annotating genomes with ko ids using DIAMOND')
                self.annotate_diamond(
                    genomes_list, self.databases.KO_DB, annotation_type, AnnotationParser.KO, self.GENOME_KO)

                logging.info('    - Generating ko frequency table')
                mg = MatrixGenerator(MatrixGenerator.KO)
                freq_table = os.path.join(self.output_directory, self.OUTPUT_KO)
                mg.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.ko_hmm:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with ko ids using HMMs')
                self.hmmsearch_annotation(genomes_list,
                                          os.path.join(
                                              self.output_directory, self.GENOME_KO_HMM),
                                          self.databases.KO_HMM_DB,
                                          AnnotationParser.KO,
                                          annotation_type)

                logging.info('    - Generating ko frequency table')
                mg = MatrixGenerator(MatrixGenerator.KO)
                freq_table = os.path.join(
                    self.output_directory, self.OUTPUT_KO_HMM)
                mg.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.ec:
                annotation_type = AnnotationParser.BLASTPARSER
                logging.info('    - Annotating genomes with ec ids')
                self.annotate_diamond(genomes_list, self.databases.EC_DB, annotation_type, AnnotationParser.EC, self.GENOME_EC)
                
                logging.info('    - Generating ec frequency table')
                mg = MatrixGenerator(MatrixGenerator.EC)
                freq_table = os.path.join(self.output_directory, self.OUTPUT_EC)
                mg.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.pfam:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with pfam ids')
                self.hmmsearch_annotation(genomes_list,
                                          os.path.join(self.output_directory, self.GENOME_PFAM),
                                          self.databases.PFAM_DB,
                                          AnnotationParser.PFAM,
                                          annotation_type)
                
                logging.info('    - Generating pfam frequency table')
                mg = MatrixGenerator(MatrixGenerator.PFAM)
                freq_table = os.path.join(self.output_directory, self.OUTPUT_PFAM)
                mg.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.tigrfam:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with tigrfam ids')
                self.hmmsearch_annotation(genomes_list,
                                          os.path.join(self.output_directory, self.GENOME_TIGRFAM),
                                          self.databases.TIGRFAM_DB,
                                          AnnotationParser.TIGRFAM,
                                          annotation_type)

                logging.info('    - Generating tigrfam frequency table')
                mg = MatrixGenerator(MatrixGenerator.TIGRFAM)
                freq_table = os.path.join(self.output_directory, self.OUTPUT_TIGRFAM)
                mg.write_matrix(genomes_list, self.count_domains, freq_table)
            
            if self.cazy:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with CAZY ids')
                self.hmmsearch_annotation(genomes_list,
                                          os.path.join(self.output_directory, self.GENOME_CAZY),
                                          self.databases.CAZY_DB,
                                          AnnotationParser.CAZY,
                                          annotation_type)

                logging.info('    - Generating CAZY frequency table')
                mg = MatrixGenerator(MatrixGenerator.CAZY)
                freq_table = os.path.join(self.output_directory, self.OUTPUT_CAZY)
                mg.write_matrix(genomes_list, self.count_domains, freq_table)

            if hasattr(list(genomes_list[0].sequences.values())[0], "prod_id"):
                logging.info('Generating .gff files:')
                self._generate_gff_files(genomes_list)

                logging.info('Renaming protein headers')
                self._rename_fasta(genomes_list)

            if not self.light:
                logging.info('Storing genome objects')
                self._pickle_objects(genomes_list)

            logging.info('Finished annotation')

