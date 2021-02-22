#!/usr/bin/env python3
# pylint: disable=line-too-long
"""
Various functions for annotating genomes.
"""
# System imports

import pickle
import shutil
import os
import tempfile
import logging
import subprocess
import multiprocessing as mp
from os import path, close, mkdir, listdir
from enrichm.genome import Genome, AnnotationParser
from enrichm.databases import Databases
from enrichm.sequence_io import SequenceIO
from enrichm.writer import Writer, MatrixGenerator
from enrichm.toolbox import list_splitter, run_command

def parse_genomes(params):
    '''
    Parses an input genome file into a Genome object. This is outside
    of a class to enable parallelisation
    '''

    genome = Genome(*params)
    return genome

class Annotate:
    '''
    Annotates proteins, and MAGs
    '''
    GENOME_BIN = 'genome_bin'
    GENOME_PROTEINS = 'genome_proteins'
    GENOME_GENES = 'genome_genes'
    GENOME_KO = 'annotations_ko'
    GENOME_KO_HMM = 'annotations_ko_hmm'
    GENOME_EC = 'annotations_ec'
    GENOME_PFAM = 'annotations_pfam'
    GENOME_TIGRFAM = 'annotations_tigrfam'
    GENOME_HYPOTHETICAL = 'annotations_hypothetical'
    GENOME_CAZY = 'annotations_cazy'
    GENOME_GFF = 'annotations_gff'
    GENOME_OBJ = 'annotations_genomes'
    OUTPUT_KO = 'ko_frequency_table.tsv'
    OUTPUT_KO_HMM = 'ko_hmm_frequency_table.tsv'
    OUTPUT_EC = 'ec_frequency_table.tsv'
    OUTPUT_PFAM = 'pfam_frequency_table.tsv'
    OUTPUT_TIGRFAM = 'tigrfam_frequency_table.tsv'
    OUTPUT_CAZY = 'cazy_frequency_table.tsv'
    OUTPUT_CLUSTER = 'cluster_frequency_table.tsv'
    OUTPUT_ORTHOLOG = 'ortholog_frequency_table.tsv'
    OUTPUT_HYPOTHETICAL_ANNOTATIONS = 'hypothetical_annotations.tsv'
    OUTPUT_DIAMOND = "DIAMOND_search"
    GFF_SUFFIX = '.gff'
    PROTEINS_SUFFIX = '.faa'
    ANNOTATION_SUFFIX = '.tsv'
    PICKLE_SUFFIX = '.pickle'

    def __init__(self, output_directory, annotate_ko, annotate_ko_hmm, annotate_pfam,
                 annotate_tigrfam, annoatate_cluster, annotate_ortholog, annotate_cazy, annotate_ec,
                 annotate_orthogroup, evalue, bit, percent_id_cutoff, aln_query, aln_reference, 
                 fraction_aligned, cut_ga_pfam, cut_nc_pfam, cut_tc_pfam, cut_ga_tigrfam, cut_nc_tigrfam,
                 cut_tc_tigrfam, cut_hmm, inflation, chunk_number, chunk_max,
                 count_domains, threads, parallel, suffix, light):


        # Define inputs and outputs
        self.output_directory = output_directory

        # Define type of annotation to be carried out
        self.annotate_ko = annotate_ko
        self.annotate_ko_hmm = annotate_ko_hmm
        self.annotate_pfam = annotate_pfam
        self.annotate_tigrfam = annotate_tigrfam
        self.annotate_cluster = annoatate_cluster
        self.annotate_ortholog = annotate_ortholog
        self.annotate_orthogroup = annotate_orthogroup
        self.annotate_cazy = annotate_cazy
        self.annotate_ec = annotate_ec

        # Cutoffs
        self.evalue = evalue
        self.bit = bit
        self.percent_id_cutoff = percent_id_cutoff
        self.aln_query = aln_query
        self.aln_reference = aln_reference
        self.fraction_aligned = fraction_aligned
        self.cut_ga_pfam = cut_ga_pfam
        self.cut_nc_pfam = cut_nc_pfam
        self.cut_tc_pfam = cut_tc_pfam
        self.cut_ga_tigrfam = cut_ga_tigrfam
        self.cut_nc_tigrfam = cut_nc_tigrfam
        self.cut_tc_tigrfam = cut_tc_tigrfam
        self.cut_hmm = cut_hmm
        self.inflation = inflation
        self.chunk_number = chunk_number
        self.chunk_max = chunk_max
        self.count_domains = count_domains

        # Parameters
        self.threads = threads
        self.parallel = parallel
        self.suffix = suffix
        self.light = light

        # Set up multiprocesses pool
        self.pool = mp.Pool(processes=int(self.parallel))

        # Load databases
        self.databases = Databases()

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
            mkdir(genome_directory)
            genome_paths = list()

            for genome_path in genome_file_list:

                if genome_path.endswith(self.suffix):
                    genome_paths.append(f"{genome_path}")

            cmd = f"xargs --arg-file=/dev/stdin cp --target-directory={genome_directory}"

            logging.debug(cmd)
            process = subprocess.Popen(["bash", "-c", cmd],
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       universal_newlines=True)
            process.communicate(input=str('\n'.join(genome_paths)))

        return genome_directory

    def call_proteins(self, genome_directory):
        '''
        Use prodigal to call proteins within the genomes

        Parameters
        ----------
        genome_directory  - string. Directory containing .fna files for each
                            input genome

        Outputs
        -------
        returns the directory containing an .faa file for each input genomes
        '''
        protein_directory_path = path.join(self.output_directory, self.GENOME_PROTEINS)
        gene_directory_path = path.join(self.output_directory, self.GENOME_GENES)
        mkdir(protein_directory_path)
        mkdir(gene_directory_path)
        genome_list = list()
        genome_paths = list()

        for genome in listdir(genome_directory):

            if genome.endswith(self.suffix):
                genome_paths.append(path.splitext(genome)[0])

        logging.info("    - Calling proteins for %i genomes", len(genome_paths))
        cmd = "ls %s/*%s | \
                    sed 's/%s//g' | \
                    grep -o '[^/]*$' | \
                    parallel -j %s \
                        prodigal \
                            -q \
                            -p meta \
                            -o /dev/null \
                            -d %s/{}%s \
                            -a %s/{}%s \
                            -i %s/{}%s \
                            > /dev/null 2>&1" \
                % (genome_directory, self.suffix, self.suffix, self.parallel, gene_directory_path,
                   self.suffix, protein_directory_path, self.PROTEINS_SUFFIX, genome_directory,
                   self.suffix)

        run_command(cmd)

        protein_directory_files = listdir(protein_directory_path)
        genome_directory_files = listdir(genome_directory)

        for genome_protein, genome_nucl in zip(protein_directory_files, genome_directory_files):
            genome_protein_base = genome_protein.replace(self.PROTEINS_SUFFIX, self.suffix)
            output_genome_protein_path = path.join(protein_directory_path, genome_protein)
            output_genome_nucl_path = path.join(genome_directory, genome_nucl)
            output_genome_gene_path = path.join(gene_directory_path, genome_protein_base)

            genome = (self.light, output_genome_protein_path, output_genome_nucl_path,
                      output_genome_gene_path)
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

        output_directory_path = path.join(self.output_directory,
                                          output_subdirectory)
        genome_dict = {genome.name:genome for genome in genomes_list}
        mkdir(output_directory_path)
        specific_cutoffs = None

        with tempfile.NamedTemporaryFile() as temp:

            to_write = str()

            for genome in genomes_list:
                to_write += f"sed \"s/>/>{genome.name}~/g\" {genome.path}\n"

            temp.write(str.encode(to_write))
            temp.flush()
            output_annotation_path = path.join(output_directory_path, self.OUTPUT_DIAMOND) + \
                                        self.ANNOTATION_SUFFIX
            logging.info('    - BLASTing genomes')
            self.diamond_search(temp.name, output_annotation_path, database)

            for genome_name, batch in self.get_batches(output_annotation_path):

                if batch:
                    genome = genome_dict[genome_name]
                    genome.add(batch, self.evalue, self.bit, self.aln_query, self.aln_reference,
                               specific_cutoffs, parser_type, ids_type)

    def get_batches(self, input_file):
        '''
        Separate DIAMOND blast results into batches, where a batch is all the hits for a genome.

        Parameters
        ----------
        input_file - string. Directory to search for blast results.
        '''

        last = None
        input_file_io = open(input_file)

        for line in input_file_io:
            split_line = line.strip().split('\t')
            genome_id = split_line[0].split('~')[0]

            if last is None:
                last = genome_id
                batch = [split_line]

            else:

                if last == genome_id:
                    batch.append(split_line)
                else:
                    yield last, batch
                    batch = [split_line]
                    last = genome_id

        if last is None:
            yield None, None
        else:
            yield last, batch

    def diamond_search(self, tmp_name, output_path, database):
        '''
        Carry out a diamond blastp search.

        Parameters
        ----------
        input_genome_path - string. Path to file containing .faa file for an input genome
        output_path - string. Path to file to output results into
        databases - string. Path to HMM to use for searching
        '''

        cmd = f'bash {tmp_name} | diamond blastp \
                                    --quiet \
                                    --outfmt 6 \
                                    --max-target-seqs 1 \
                                    --query /dev/stdin \
                                    --out {output_path} \
                                    --db {database} \
                                    --threads {self.threads} '
        if self.evalue:
            cmd += f'--evalue {self.evalue} '

        if self.bit:
            cmd += f'--min-score {self.bit} '

        if self.percent_id_cutoff:
            cmd += f'--id {self.percent_id_cutoff*100} '

        if self.aln_query:
            cmd += f"--query-cover {self.aln_query*100} "

        if self.aln_reference:
            cmd += f"--subject-cover {self.aln_reference*100} "

        run_command(cmd)

    def hmmsearch_annotation(self, genomes_list, output_directory_path, database, ids_type, parser):
        '''
        Annotate the proteins encoded by each genome with pfam ids using HMM searches.

        Parameters
        ----------
        genomes_list - list. list of Genome objects

        '''
        mkdir(output_directory_path)
        genome_dict = {genome.name: genome for genome in genomes_list}

        hmmcutoff = (ids_type in (AnnotationParser.TIGRFAM, AnnotationParser.PFAM))

        if ids_type == AnnotationParser.KO_HMM:
            specific_cutoffs = self.databases.parse_ko_cutoffs()
        else:
            specific_cutoffs = None

        self.hmm_search(output_directory_path, database, hmmcutoff)
        
        if ids_type == AnnotationParser.PFAM:
            pfam2clan = self.databases.pfam2clan()
        else:
            pfam2clan = None

        for genome_annotation in listdir(output_directory_path):
            genome_id = path.splitext(genome_annotation)[0]
            genome = genome_dict[genome_id]
            output_annotation_path = path.join(output_directory_path, genome_annotation)
            genome.add(output_annotation_path, self.evalue, self.bit, self.aln_query,
                       self.aln_reference, specific_cutoffs, parser, ids_type,
                       pfam2clan=pfam2clan)

    def annotate_hypothetical(self, genomes_list):
        '''
        Sort proteins coded by each genome into homologous clusters.

        Inputs
        ------
        genomes_list - list. list of Genome objects

        '''
        output_directory_path = path.join(self.output_directory, self.GENOME_HYPOTHETICAL)
        mkdir(output_directory_path)

        renamed_genomes = list()
        for genome in genomes_list:
            renamed_genome = next(tempfile._get_candidate_names())
            cmd = f"sed 's/>/>{genome.name}~/g' {genome.path} > {renamed_genome}"
            run_command(cmd)
            renamed_genomes.append(renamed_genome)

        
        tmp_dir = tempfile.mkdtemp()

        db_path = path.join(output_directory_path, "db")
        clu_path = path.join(output_directory_path, "clu")
        align_path = path.join(output_directory_path, "alignDb")
        blast_output_path = path.join(output_directory_path, "alignDb.m8")
        formatted_blast_output_path = path.join(output_directory_path, "alignDb.formatted.m8")

        clu_tsv_path = path.join(output_directory_path, "hypothetical_clusters.tsv")

        logging.info('    - Generating MMSeqs2 database')
        cmd = f"mmseqs createdb {' '.join(renamed_genomes)} {db_path}"
        run_command(cmd)
        for renamed_genome in renamed_genomes:
            os.remove(renamed_genome)

        logging.info('    - Clustering genome proteins')
        cmd = f"mmseqs cluster \
                    {db_path} \
                    {clu_path} \
                    {tmp_dir} \
                    --threads {self.threads} \
                    --min-seq-id {self.percent_id_cutoff} \
                    -c {self.fraction_aligned} \
                    -v 0"
        run_command(cmd)

        logging.info('    - Extracting clusters')
        cmd = f'mmseqs createtsv \
                    {db_path} \
                    {db_path} \
                    {clu_path} \
                    {clu_tsv_path} \
                    --threads {self.threads} \
                    -v 0'
        run_command(cmd)

        if self.annotate_ortholog:

            logging.info('    - Computing Smith-Waterman alignments for clustering results')
            cmd = f"mmseqs alignall \
                        {db_path} \
                        {clu_path} \
                        {align_path} \
                        --alignment-mode 3 \
                        --threads {self.threads} \
                        -v 0"
            run_command(cmd)

            logging.info('    - Converting to BLAST-like output')
            cmd = f"mmseqs createtsv \
                        {db_path} \
                        {db_path} \
                        {align_path} \
                        {blast_output_path} \
                        --threads {self.threads} \
                        -v 0"
            # --format-output query,target,bits
            run_command(cmd)

            logging.info('    - Reformatting BLAST output')
            cmd = "OFS=\"\t\" awk 'FNR==NR{a[$1]=$2;next}{$3=a[$3]; \
                                            $1=\"\"; for(i=2;i<NF;i++){printf(\"%s\t\",$i)} \
                                            printf(\"\\n\")}' %s %s | cut -f1,2,5 > %s" \
                % ("%s", db_path + '.lookup', blast_output_path, formatted_blast_output_path)
            run_command(cmd)

            ortholog_dict = self.run_mcl(formatted_blast_output_path,
                                            output_directory_path)
            ortholog_ids = ortholog_dict.keys()
        else:
            ortholog_dict = dict()
            ortholog_ids = list()
        cluster_ids = self.parse_cluster_results(clu_tsv_path,
                                                 genomes_list,
                                                 ortholog_dict,
                                                 output_directory_path)
        return cluster_ids, ortholog_ids

    def run_mcl(self, blast_abc, output_directory_path):
        '''
        Parse the protein clusters produced from Mmseqs2 using mcl

        Parameters
        ----------
        blast_abc - string. an abc file for mcl to run on. More information on the format of abc
                    files can be found at https://micans.org/mcl/man/clmprotocols.html
        output_directory_path - string. Path to write the results of mcl parsing to.
        '''

        dict_path = path.join(output_directory_path, "alignDb.dict")
        mci_path = path.join(output_directory_path, "alignDb.mci")
        cluster_path = path.join(output_directory_path, "mcl_clusters.tsv")
        output_path = path.join(output_directory_path, "mcl_clusters.convert.tsv")
        
        logging.info('    - Preparing network')
        ortholog_dict = dict()
        cmd = f"mcxload \
                    -abc {blast_abc} \
                    -write-tab {dict_path} \
                    -o {mci_path} \
                    --stream-mirror \
                    --stream-neg-log10"
        run_command(cmd)
        logging.info('    - Finding orthologs')
        ortholog_dict = dict()
        cmd = f'mcl \
                    {mci_path} \
                    -te {self.threads} \
                    -I {self.inflation} \
                    -o {cluster_path}'
        run_command(cmd)

        logging.info('    - Reformatting output')
        ortholog_dict = dict()
        cmd = f'mcxdump \
                    -icl {cluster_path} \
                    -o {output_path} \
                    -tabr {dict_path}'
        run_command(cmd)

        ortholog = 1
        for line in open(output_path):
            ortholog_idx = "ortholog_%i" % ortholog
            ortholog_dict[ortholog_idx] = set()

            for protein in line.strip().split('\t'):
                ortholog_dict[ortholog_idx].add(protein)

            ortholog += 1
        return ortholog_dict

    def parse_cluster_results(self, cluster_output_path, genomes_list, ortholog_dict,
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
        logging.info('    - Parsing input cluster file: %s', cluster_output_path)

        cluster_ids = set()
        previous_cluster_name = None
        counter = 0
        genome_dictionary = {genome.name:genome for genome in genomes_list}

        output_hypothetical_annotations = path.join(output_directory_path,
                                                    self.OUTPUT_HYPOTHETICAL_ANNOTATIONS)
        with open(output_hypothetical_annotations, 'w') as out_io:

            for line in open(cluster_output_path):

                cluster_id, member = line.strip().split('\t')
                genome_id, sequence_id = member.split('~')

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
        cmd = str()

        if self.bit:
            cmd += '-T %s ' % (str(self.bit))
        else:
            cmd += '-E %s ' % (str(self.evalue))

        return cmd

    def hmm_search(self, output_path, database, hmmcutoff):
        '''
        Carry out a hmmsearch.

        Parameters
        ----------
        input_genome_path     - string. Path to file containing .faa file for
                                an input genome
        output_path           - string. Path to file to output results into
        databases             - string. Path to HMM to use for searching
        '''

        input_genome_path = path.join(self.output_directory, self.GENOME_PROTEINS)
        cmd = "ls %s | sed 's/%s//g' | parallel -j %s\
                                                hmmsearch \
                                                    --cpu %s \
                                                    -o /dev/null \
                                                    --noali \
                                                    --domtblout %s/{}%s " \
                          % (input_genome_path, self.PROTEINS_SUFFIX, self.parallel,
                             self.threads, output_path, self.ANNOTATION_SUFFIX)
        if hmmcutoff:
            if (self.cut_ga_pfam or self.cut_nc_pfam or self.cut_tc_pfam) and 'pfam' in database:
                if self.cut_ga_pfam:
                    cmd += " --cut_ga "
                if self.cut_nc_pfam:
                    cmd += " --cut_nc "
                if self.cut_tc_pfam:
                    cmd += " --cut_tc "
            elif (self.cut_ga_tigrfam or self.cut_nc_tigrfam or self.cut_tc_tigrfam) and 'tigrfam' in database:
                if self.cut_ga_tigrfam:
                    cmd += " --cut_ga "
                if self.cut_nc_tigrfam:
                    cmd += " --cut_nc "
                if self.cut_tc_tigrfam:
                    cmd += " --cut_tc "
            else:
                cmd += self._default_hmmsearch_options()
        else:
            cmd += self._default_hmmsearch_options()

        cmd += "%s %s/{}.faa 2> /dev/null" % (database, input_genome_path)

        run_command(cmd)

    def generate_gff_files(self, genomes_list):
        '''
        Write GFF files for each of the genome objects in genomes_list

        Parameters
        ----------
        genomes_list - List. List of Genome objects
        '''
        output_directory_path = path.join(self.output_directory,
                                          self.GENOME_GFF)
        mkdir(output_directory_path)
        for genome in genomes_list:
            logging.info('    - Generating .gff file for %s', genome.name)
            gff_output = path.join(output_directory_path, genome.name + self.GFF_SUFFIX)
            Writer.write_gff(genome, gff_output)

    def rename_fasta(self, genomes_list):
        '''
        Rename the called proteins with annotation ids.

        Parameters
        ----------
        genomes_list - List. List of Genome objects
        '''
        seqio = SequenceIO()

        for genome in genomes_list:
            file_object, fname = tempfile.mkstemp(suffix='.faa', text=True)

            if genome.gene:
                fd_gene, fname_gene = tempfile.mkstemp(suffix='.fna', text=True)

                with open(fname_gene, 'w') as out_gene_io:

                    for description, sequence in seqio.each(open(genome.gene)):
                        name = description.partition(' ')[0]
                        annotations = ' '.join(genome.sequences[name].all_annotations())
                        out_gene_io.write(">%s %s\n" % (name, annotations))
                        out_gene_io.write(sequence + '\n')

                close(fd_gene)
                logging.debug('Moving %s to %s', fname_gene, genome.gene)
                shutil.move(fname_gene, genome.gene)

            with open(fname, 'w') as out_io:

                for description, sequence in seqio.each(open(genome.path)):
                    name = description.partition(' ')[0]
                    annotations = ' '.join(genome.sequences[name].all_annotations())
                    out_io.write(">%s %s\n" % (name, annotations))
                    out_io.write(str(sequence) + '\n')

            close(file_object)
            logging.debug('Moving %s to %s', fname, genome.path)
            shutil.move(fname, genome.path)

    def pickle_objects(self, genomes_list):
        '''
        Store annotated genome objects as pickles.

        Parameters
        ----------
        genomes_list - List. List of Genome objects
        '''
        output_directory_path = path.join(self.output_directory,
                                          self.GENOME_OBJ)
        mkdir(output_directory_path)

        for genome in genomes_list:
            genome_pickle_path = path.join(output_directory_path, genome.name + self.PICKLE_SUFFIX)

            with open(genome_pickle_path, 'wb') as output:
                pickle.dump(genome, output)

    def parse_genome_inputs(self, genome_directory, protein_directory, genome_files, protein_files):
        '''
        Inputs
        ------

        Outputs
        -------

        '''

        prep_genomes_list = list()
        genomes_list = list()

        if protein_directory:
            logging.info("Using provided proteins")
            protein_genome_list = list()

            for protein_file in listdir(protein_directory):
                protein_genome_list.append(path.join(protein_directory, protein_file))

            directory = self.prep_genome(protein_genome_list,
                                         path.join(self.output_directory,
                                                   self.GENOME_PROTEINS))

            for genome_proteins_file in listdir(directory):

                if genome_proteins_file.endswith(self.suffix):
                    genome = (self.light, path.join(directory, genome_proteins_file), None, None)
                    prep_genomes_list.append(genome)

        elif protein_files:
            logging.info("Using provided proteins")
            genome_proteins_path = path.join(self.output_directory, self.GENOME_PROTEINS)
            directory = self.prep_genome(protein_files, genome_proteins_path)

            for protein_file in listdir(directory):
                protein_file_path = path.join(directory, path.basename(protein_file))
                prep_genomes_list.append((self.light, protein_file_path, None, None))

        elif genome_directory:
            logging.info("Calling proteins for annotation")
            prep_genomes_list = self.call_proteins(genome_directory)
            directory = genome_directory

        elif genome_files:
            logging.info("Calling proteins for annotation")
            directory = self.prep_genome(genome_files,
                                         path.join(self.output_directory, self.GENOME_BIN))
            prep_genomes_list = self.call_proteins(directory)

        for chunk in list_splitter(prep_genomes_list, self.chunk_number, self.chunk_max):
            genomes_list += self.pool.map(parse_genomes, chunk)

        return genomes_list

    def annotate_pipeline(self, genome_directory, protein_directory, genome_files, protein_files):
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
        genomes_list = self.parse_genome_inputs(genome_directory, protein_directory,
                                                genome_files, protein_files)

        if genomes_list:
            logging.info("Starting annotation:")

            if (self.annotate_cluster or self.annotate_ortholog):
                logging.info('    - Annotating genomes with hypothetical clusters')
                cluster_ids, ortholog_ids = self.annotate_hypothetical(genomes_list)

                logging.info('    - Generating hypotheticals frequency table')
                matrix_generator = MatrixGenerator(MatrixGenerator.HYPOTHETICAL, cluster_ids)
                freq_table = path.join(self.output_directory, self.OUTPUT_CLUSTER)
                matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

                if self.annotate_ortholog:
                    matrix_generator = MatrixGenerator(MatrixGenerator.ORTHOLOG, ortholog_ids)
                    freq_table = path.join(self.output_directory, self.OUTPUT_ORTHOLOG)
                    matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.annotate_orthogroup:
                logging.warning(f"Not yet implemented")
                #self.annotate_orthogroup(genomes_list)

            if self.annotate_ko:
                annotation_type = AnnotationParser.BLASTPARSER
                logging.info('    - Annotating genomes with ko ids using DIAMOND')
                self.annotate_diamond(genomes_list, self.databases.KO_DB,
                                      annotation_type, AnnotationParser.KO,
                                      self.GENOME_KO)

                logging.info('    - Generating ko frequency table')
                matrix_generator = MatrixGenerator(MatrixGenerator.KO)
                freq_table = path.join(self.output_directory, self.OUTPUT_KO)
                matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.annotate_ko_hmm:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with ko ids using HMMs')
                self.hmmsearch_annotation(genomes_list,
                                          path.join(
                                              self.output_directory, self.GENOME_KO_HMM),
                                          self.databases.KO_HMM_DB,
                                          AnnotationParser.KO,
                                          annotation_type)

                logging.info('    - Generating ko frequency table')
                matrix_generator = MatrixGenerator(MatrixGenerator.KO)
                freq_table = path.join(
                    self.output_directory, self.OUTPUT_KO_HMM)
                matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.annotate_ec:
                annotation_type = AnnotationParser.BLASTPARSER
                logging.info('    - Annotating genomes with ec ids')
                self.annotate_diamond(genomes_list, self.databases.EC_DB, annotation_type,
                                      AnnotationParser.EC, self.GENOME_EC)

                logging.info('    - Generating ec frequency table')
                matrix_generator = MatrixGenerator(MatrixGenerator.EC)
                freq_table = path.join(self.output_directory, self.OUTPUT_EC)
                matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.annotate_pfam:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with pfam ids')
                self.hmmsearch_annotation(genomes_list,
                                          path.join(self.output_directory, self.GENOME_PFAM),
                                          self.databases.PFAM_DB,
                                          AnnotationParser.PFAM,
                                          annotation_type)

                logging.info('    - Generating pfam frequency table')
                matrix_generator = MatrixGenerator(MatrixGenerator.PFAM)
                freq_table = path.join(self.output_directory, self.OUTPUT_PFAM)
                matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.annotate_tigrfam:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with tigrfam ids')
                self.hmmsearch_annotation(genomes_list,
                                          path.join(self.output_directory, self.GENOME_TIGRFAM),
                                          self.databases.TIGRFAM_DB,
                                          AnnotationParser.TIGRFAM,
                                          annotation_type)

                logging.info('    - Generating tigrfam frequency table')
                matrix_generator = MatrixGenerator(MatrixGenerator.TIGRFAM)
                freq_table = path.join(self.output_directory, self.OUTPUT_TIGRFAM)
                matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

            if self.annotate_cazy:
                annotation_type = AnnotationParser.HMMPARSER
                logging.info('    - Annotating genomes with CAZY ids')
                self.hmmsearch_annotation(genomes_list,
                                          path.join(self.output_directory, self.GENOME_CAZY),
                                          self.databases.CAZY_DB,
                                          AnnotationParser.CAZY,
                                          annotation_type)

                logging.info('    - Generating CAZY frequency table')
                matrix_generator = MatrixGenerator(MatrixGenerator.CAZY)
                freq_table = path.join(self.output_directory, self.OUTPUT_CAZY)
                matrix_generator.write_matrix(genomes_list, self.count_domains, freq_table)

            if hasattr(list(genomes_list[0].sequences.values())[0], "prod_id"):
                logging.info('Generating .gff files:')
                self.generate_gff_files(genomes_list)

            logging.info('Renaming protein headers')
            self.rename_fasta(genomes_list)

            if not self.light:
                logging.info('Storing genome objects')
                self.pickle_objects(genomes_list)

            logging.info('Finished annotation')

        else:
            logging.error('No files found with %s suffix in input directory', self.suffix)
            raise Exception("No input files found")
