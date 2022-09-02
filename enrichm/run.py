#!/usr/bin/env python3

import logging
import sys
import os
import shutil
import time

from enrichm.data import Data
from enrichm.network_analyzer import NetworkAnalyser
from enrichm.enrichment import Enrichment
from enrichm.annotate import Annotate
from enrichm.classifier import Classify
from enrichm.generate import GenerateModel
from enrichm.predict import Predict
from enrichm.uses import Uses

####################################################################################################

debug = {1:logging.CRITICAL, 2:logging.ERROR, 3:logging.WARNING, 4:logging.INFO, 5:logging.DEBUG}

####################################################################################################

class Run:

    def __init__(self):

        self.DATA = 'data'
        self.ANNOTATE = 'annotate'
        self.CLASSIFY = 'classify'
        self.ENRICHMENT = 'enrichment'
        self.PREDICT = 'predict'
        self.GENERATE = 'generate'
        self.USES = 'uses'

    def _logging_setup(self, args):
        if args.verbosity not in range(1, 6):
            raise Exception("Logging verbosity must be a positive integer between 1 and 5.")

        logger = logging.getLogger('')
        logger.setLevel(debug[args.verbosity])
        log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                       datefmt="%Y-%m-%d %H:%M:%S %p")

        stream_logger = logging.StreamHandler(sys.stdout)
        stream_logger.setFormatter(log_format)
        stream_logger.setLevel(debug[args.verbosity])
        logger.addHandler(stream_logger)

        if args.subparser_name!=self.DATA:
            file_logger = logging.FileHandler(os.path.join(args.output, args.log), 'a')
            file_logger.setFormatter(log_format)
            logger.addHandler(file_logger)

    def _check_general(self, args):
        '''
        Check general input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        dependencies = {'hmmsearch':"http://hmmer.org/download.html",
                        'diamond':"https://github.com/bbuchfink/diamond",
                        'R':"https://www.r-project.org",
                        'parallel':"https://www.gnu.org/software/parallel",
                        'prodigal':"https://github.com/hyattpd/Prodigal/wiki/installation",
                        'mmseqs':"https://github.com/soedinglab/MMseqs2"}

        missing_dependencies = list()

        for dependency in dependencies.keys():

            if shutil.which(dependency) == None:
                missing_dependencies.append(dependency)

        if len(missing_dependencies)>0:
            dependency_string = '\n'.join(['\t%s\t%s' % (d, dependencies[d]) for d in missing_dependencies])
            raise Exception('The following dependencies need to be installed to run enrichm:\n%s' % (dependency_string))

        if args.subparser_name != self.DATA:
            # Set up working directory
            if not args.output:
                args.output = '%s-enrichm_%s_output' % (time.strftime("%Y-%m-%d_%H-%M"), args.subparser_name)

            if(os.path.isdir(args.output) or os.path.isfile(args.output)):

                if args.force:

                    if os.path.isdir(args.output):
                        shutil.rmtree(args.output)

                    else:
                        os.remove(args.output)

                else:
                    raise Exception("File '%s' exists." % args.output)

            os.mkdir(args.output)

    def _check_data(self, args):
        if not(args.create or args.uninstall):
            raise Exception("Only one of the following can be specified: --create, --uninstall")

        if not(os.access(Data.DATABASE_DIR, os.R_OK|os.W_OK)):
            raise Exception(f"EnrichM does not have read/write in database directory: {Data.DATABASE_DIR}")

    def _check_annotate(self, args):
        '''
        Check annotate input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        # ensure either a list of genomes or a directory of genomes have been specified
        if not(args.genome_files or args.genome_directory or args.protein_directory or args.protein_files):
            raise Exception("Input error: Either a list of genomes or a directory of genomes need to be specified.")

        if len([x for x in [args.genome_files, args.genome_directory, args.protein_directory, args.protein_files] if x]) != 1:
            raise Exception("Input error: Only one type of input can be specified (--genome_files, --genome_directory, --protein_directory, or --protein_files).")

        if not args.suffix:

            if(args.genome_directory or args.genome_files):
                args.suffix = '.fna'

            elif(args.protein_directory or args.protein_files):
                args.suffix = '.faa'

        if(args.id>1 or args.id<0):
            raise Exception("Identity (--id) must be between 0 and 1.")

        if(args.aln_query>1 or args.aln_query<0):
            raise Exception("Alignment to query cutoff (--aln_query) must be between 0 and 1")

        if(args.aln_reference>1 or args.aln_reference<0):
            raise Exception("Alignment to reference cutoff (--aln_reference) must be between 0 and 1")

        if any([args.cut_ga_pfam, args.cut_nc_pfam, args.cut_tc_pfam]):
            if len([x for x in [args.cut_ga_pfam, args.cut_nc_pfam, args.cut_tc_pfam] if x])>1:
                raise Exception("Only one of the following can be selected: --cut_ga_pfam, --cut_nc_pfam, --cut_tc_pfam")
            if args.evalue:
                logging.warning('selecting one of the following overrides evalue thresholds: --cut_ga_pfam, --cut_nc_pfam, --cut_tc_pfam, --cut_ga_tigrfam, --cut_nc_tigrfam, --cut_tc_tigrfam')

        if any([args.cut_ga_tigrfam, args.cut_nc_tigrfam, args.cut_tc_tigrfam]):
            if len([x for x in [args.cut_ga_tigrfam, args.cut_nc_tigrfam, args.cut_tc_tigrfam] if x])>1:
                raise Exception("Only one of the following can be selected: --cut_ga_tigrfam, --cut_nc_tigrfam, --cut_tc_tigrfam")
            if args.evalue:
                logging.warning('selecting one of the following overrides evalue thresholds: --cut_ga_pfam, --cut_nc_pfam, --cut_tc_pfam, --cut_ga_tigrfam, --cut_nc_tigrfam, --cut_tc_tigrfam')

    def _check_enrichment(self, args):
        '''
        Check enrichment input and output options are valid.

        Parameters

        args    - object. Argparse object

        Output
        ------
        '''
        ### ~ TODO: Check Multi test correction inputs...
        types = [args.ko, args.pfam, args.tigrfam, args.cluster, args.ortholog, args.cazy, args.ec, args.ko_hmm]

        if not args.abundance and args.abundance_metadata:
           raise Exception("Values for both --abundance and --abundance_metadata are required")
        if not args.transcriptome and args.transcriptome_metadata:
           raise Exception("Values for both --abundance and --abundance_metadata are required")
        if args.annotation_matrix and args.annotate_output:
            raise Exception("Use either --annotate_output or --annotation_matrix")

        if(not args.annotation_matrix and not args.annotate_output and not args.gff_files):
            raise Exception("Either --annotate_output, --annotation_matrix or --gff_files must be specified!")

        if args.annotation_matrix or args.annotate_output:
            if not args.abundance:
                if not args.metadata:
                    raise Exception("Genome groups need to be specified using the --metadata flag")

        if args.annotate_output:

            if not any(types):
                raise Exception(
                    "Input Error: One of the following flags must be specified: --ko --pfam --tigrfam --orthologs --orthogroup --clusters --ko_hmm --ec --cazy")

            if len([x for x in types if x]) > 1:
                raise Exception(
                    "Only one of the following flags may be specified: --ko --pfam --tigrfam --orthologs --orthogroup --clusters --ko_hmm --ec --cazy")

    def _check_classify(self, args):
        '''
        Check classify input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        if(args.cutoff>1 and args.cutoff<0):
            raise Exception("--cutoff needs to be between 0 and 1")

        if(args.gff_files and args.genome_and_annotation_matrix):
            raise Exception(f"Both --gff_files {args.gff_files} and --genome_and_annotation_matrix {args.genome_and_annotation_matrix} were specified. Please provide only one option.")

        if(args.module_rules_json and not args.gff_files):
            raise Exception(f"--gff_files must be provided to use --module_rules_json")

        if args.module_rules_json:
            if not os.path.isfile(args.module_rules_json):
                raise Exception(f"File does not exist: {args.module_rules_json}")

        if args.custom_modules:
            if not os.path.isfile(args.custom_modules):
                raise Exception(f"File does not exist: {args.custom_modules}")

    def _check_network(self, args):
        '''
        Check network (explore, pathway) input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        if not hasattr(args, 'enrichment_output'):
            args.enrichment_output = None

        if any([args.abundance, args.abundance_metadata]):
            if not (args.abundance and args.abundance_metadata):
                raise Exception("Both abundance and abundance metadata need to be specified")

        if any([args.tpm_values, args.tpm_metadata]):
            if not (args.tpm_values and args.tpm_metadata):
                raise Exception("Both --tpm_values and --tpm_metadata need to be specified")

        if args.subparser_name == NetworkAnalyser.PATHWAY:
            args.depth = None
            args.queries = None

        if args.subparser_name == NetworkAnalyser.EXPLORE:
            args.filter = None
            args.limit = None

            if not(args.queries):

                if args.depth:
                    logging.warning("--depth argument ignored without --queries flag")

    def _check_predict(self, args):
        '''
        Check general input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        pass

    def _check_uses(self, args):
        '''
        Check general input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        pass

    def _check_generate(self, args):
        '''
        Check general input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        pass

    def run_data(self, args):
        self._check_data(args)
        d = Data()
        d.do(args.uninstall, args.create)

    def run_annotate(self, args):
        self._check_annotate(args)

        annotate = Annotate(# Define inputs and outputs
                            args.output,
                            # Define type of annotation to be carried out
                            args.ko, args.ko_hmm, args.pfam, args.tigrfam,
                            args.clusters, args.orthologs, args.cazy,
                            args.ec, args.orthogroup,
                            # Cutoffs
                            args.evalue, args.bit, args.id, args.aln_query,
                            args.aln_reference, args.c, args.cut_ga_pfam,

                            args.cut_nc_pfam, args.cut_tc_pfam,
                            args.cut_ga_tigrfam, args.cut_nc_tigrfam,
                            args.cut_tc_tigrfam, args.cut_ko,
                            args.inflation, args.chunk_number, args.chunk_max,
                            args.count_domains,
                            # Parameters
                            args.threads, args.parallel, args.suffix, args.light)

        annotate.annotate_pipeline(args.genome_directory,
                                    args.protein_directory,
                                    args.genome_files,
                                    args.protein_files)

    def run_classify(self, args):
        self._check_classify(args)
        skip_database_check  = False
        if args.custom_modules:
            skip_database_check = True
        classify = Classify(skip_database_check)
        classify.classify_pipeline(args.custom_modules, args.cutoff, args.aggregate,
                                    args.genome_and_annotation_matrix, args.module_rules_json, 
                                    args.gff_files, args.output)

    def run_enrichment(self, args):
        self._check_enrichment(args)
        enrichment = Enrichment()
        enrichment.enrichment_pipeline(# Input options
                                        args.annotate_output, args.annotation_matrix, args.gff_files,
                                        args.metadata, args.abundance, args.abundance_metadata,
                                        args.transcriptome, args.transcriptome_metadata,
                                        # Runtime options
                                        args.pval_cutoff, args.proportions_cutoff, 
                                        args.threshold, args.multi_test_correction, 
                                        args.batchfile, args.processes, 
                                        args.allow_negative_values, args.ko, args.pfam, 
                                        args.tigrfam, args.cluster, args.ortholog, args.cazy,
                                        args.ec, args.ko_hmm, args.range, args.subblock_size,
                                        args.operon_mismatch_cutoff, args.operon_match_score_cutoff,
                                        # Outputs
                                        args.output)


    def run_network(self, args):
        self._check_network(args)
        network_analyser=NetworkAnalyser()
        network_analyser.network_pipeline(args.subparser_name, args.matrix, 
                                            args.genome_metadata, args.tpm_values,
                                            args.tpm_metadata, args.abundance, 
                                            args.abundance_metadata, args.metabolome,
                                            args.enrichment_output, args.depth, args.filter,
                                            args.limit, args.queries, args.output)


    def run_predict(self, args):
        self._check_predict(args)
        predict = Predict()
        predict.predict_pipeline(args.forester_model_directory,
                args.input_matrix,
                args.output)

    def run_generate(self, args):
        self._check_generate(args)
        generate_model = GenerateModel()
        generate_model.generate_pipeline(args.input_matrix,
                args.groups,
                args.model_type,
                args.testing_portion,
                args.grid_search,
                args.threads,
                args.output)


    def run_uses(self, args):
        self._check_uses(args)
        uses = Uses()
        uses.uses_pipeline(args.compounds_list,
                args.annotation_matrix,
                args.metadata,
                args.output,
                args.count)

    def get_pipeline(self, subparser_name):
        if subparser_name == self.DATA:
            pipeline = self.run_data
        elif subparser_name == self.ANNOTATE:
            pipeline = self.run_annotate
        elif subparser_name == self.CLASSIFY:
            pipeline = self.run_classify
        elif subparser_name == self.ENRICHMENT:
            pipeline = self.run_enrichment
        elif subparser_name == NetworkAnalyser.PATHWAY:
            pipeline = self.run_network
        elif subparser_name == NetworkAnalyser.EXPLORE:
            pipeline = self.run_network
        elif subparser_name == self.PREDICT:
            pipeline = self.run_predict
        elif subparser_name == self.GENERATE:
            pipeline = self.run_generate
        elif subparser_name == self.USES:
            pipeline = self.run_uses
        else:
            ValueError(subparser_name)


        return pipeline

    def run_enrichm(self, args, command):
        '''
        Parameters
        ----------

        Output
        ------
        '''
        self._check_general(args)
        self._logging_setup(args)

        logging.info("Command: %s" % ' '.join(command))
        logging.info("Running the %s pipeline" % args.subparser_name)
        pipeline = self.get_pipeline(args.subparser_name)
        pipeline(args)
        logging.info('Finished running EnrichM')
