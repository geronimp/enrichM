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

__author__      = "Joel Boyd"
__copyright__   = "Copyright 2018"
__credits__     = ["Joel Boyd"]
__license__     = "GPL3"
__version__     = "0.0.7"
__maintainer__  = "Joel Boyd"
__email__       = "joel.boyd near uq.net.au"
__status__      = "Development"

###############################################################################
# Imports
import logging
import sys
import os
import shutil
import time

# Local
from enrichm.data import Data
from enrichm.network_analyzer import NetworkAnalyser
from enrichm.enrichment import Enrichment
from enrichm.annotate import Annotate
from enrichm.classifier import Classify
from enrichm.generate import GenerateModel
from enrichm.predict import Predict
from enrichm.connect import Connect
from enrichm.aggregate import Aggregate
###############################################################################

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################

class Run:
    
    def __init__(self):
        
        self.ANNOTATE        = 'annotate'
        self.COMPARE         = 'compare'

        self.CLASSIFY        = 'classify'
        self.BUILD           = 'build'
        self.ENRICHMENT      = 'enrichment'
        self.MODULE_AB       = 'module_ab'
        self.DATA            = 'data'
        self.PREDICT         = 'predict'
        self.GENERATE        = 'generate'
        self.CONNECT         = 'connect'
        self.AGGREGATE       = 'aggregate'

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

        dependencies = {'hmmsearch': "http://hmmer.org/download.html",
                        'diamond': "https://github.com/bbuchfink/diamond",
                        'R': "https://www.r-project.org",
                        'parallel': "https://www.gnu.org/software/parallel",
                        'prodigal': "https://github.com/hyattpd/Prodigal/wiki/installation",
                        'mmseqs': "https://github.com/soedinglab/MMseqs2"}
        
        missing_dependencies = list()
        
        for dependency in dependencies.keys():

            if shutil.which(dependency) == None:
                missing_dependencies.append(dependency)
        
        if len(missing_dependencies)>0:
            dependency_string = '\n'.join(['\t%s\t%s' % (d, dependencies[d]) for d in missing_dependencies])
            raise Exception('The following dependencies need to be installed to run enrichm:\n%s' % (dependency_string))

        # We dont need an output directory for the DATA pipeline
        if args.subparser_name!=self.DATA:
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

    def _check_annotate(self, args):
        '''
        Check annotate input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        if args.cut_ko:
            if int(Data.CURRENT_VERSION.split('_')[-1].replace('v', '')) < 10:
                raise Exception("EnrichM database needs to be version 10 or higher to use KO HMM cutoffs. Please run enrichm data.")

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
        
        if any([args.cut_ga, args.cut_nc, args.cut_tc]):
            if len([x for x in [args.cut_ga, args.cut_nc, args.cut_tc] if x])>1:
                raise Exception("Only one of the following can be selected: --cut_ga, --cut_nc, --cut_tc")
            
            if args.evalue:
                logging.warning('selecting one of the following overrides evalue thresholds: --cut_ga, --cut_nc, --cut_tc')
            
    def _check_enrichment(self, args):
        '''
        Check enrichment input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object

        Output
        ------
        '''
        ### ~ TODO: Check Multi test correction inputs...
        types = [args.ko, args.pfam, args.tigrfam, args.hypothetical, args.cazy, args.ec, args.ko_hmm]
        
        if not any(types):
            
            raise Exception("Input Error: One of the following flags must be specified: --ko --pfam --tigrfam --hypothetical --cazy")
        
        if len([x for x in types if x])>1:
            
            raise Exception("Only one of the following flags may be specified: --ko --pfam --tigrfam --hypothetical --cazy")

        if not args.abundance and args.abundance_metadata:
           raise Exception("values for both --abundance and --abundance_metadata are required") 

    def _check_classify(self, args):  
        '''
        Check classify input and output options are valid.
        
        Parameters
        ----------
        args    - object. Argparse object
        '''
        # Ensure either an annotation matrix or list file has been specified:
        if not(args.genome_and_annotation_file or args.genome_and_annotation_matrix):
            raise Exception("Input error: An input file must be specified to either \
--genome_and_annotation_file or --genome_and_annotation_matrix")

        elif(args.aggregate and args.genome_and_annotation_file):
            raise Exception("--aggregate needs to be run with the genome and annotation matrix")

    def _check_build(self, args):
        '''
        Check build input and output options are valid.
        
        Parameters
        ----------
        args    - object. Argparse object
        '''
        if not(args.metadata or args.abundances):
            
            raise Exception("No metadata or abundance information provided to build.")

    def _check_network(self, args):
        '''
        Check network (explore, pathway) input and output options are valid.
        
        Parameters
        ----------
        args    - object. Argparse object
        '''

        if any([args.abundance, args.abundance_metadata]):
            if not (args.abundance and args.abundance_metadata):
                raise Exception("Both abundance and abundance metadata need to be specified")

        if args.subparser_name==NetworkAnalyser.PATHWAY:
            args.depth              = None
            args.queries            = None
            args.starting_compounds = None
            args.steps              = None
            args.number_of_queries  = None
        
        if args.subparser_name==NetworkAnalyser.TRAVERSE:
            args.depth              = None
            args.queries            = None
    
        if args.subparser_name==NetworkAnalyser.EXPLORE:
            args.filter             = None  
            args.limit              = None 
            args.starting_compounds = None
            args.steps              = None
            args.number_of_queries  = None

            if not(args.queries):

                if args.depth:
                    logging.warning("--depth argument ignored without --queries flag")

    def _check_predict(self, args):
        '''
        Inputs
        ------
        
        Outputs
        -------
        '''
        pass

    def _check_generate(self, args):
        '''
        Parameters
        ----------
        
        Output
        ------
        '''
        pass

    def _check_connect(self, args):
        '''
        Check connect
        
        Inputs
        ------
        
        Outputs
        -------
        
        '''
        if(args.cutoff > 1.0 and args.cutoff < 0.0):
           
           raise Exception("Cutoff needs to be between 0 - 1") 
    
    def main(self, args, command):
        '''
        Parameters
        ----------
        
        Output
        ------
        '''
        self._check_general(args)
        self._logging_setup(args)

        logging.info("Running command: %s" % ' '.join(command))

        if args.subparser_name == self.DATA:
            d = Data()
            d.do(args.uninstall)
        
        if args.subparser_name == self.ANNOTATE:
            self._check_annotate(args)
            a = Annotate(# Define inputs and outputs
                         args.output,
                         # Define type of annotation to be carried out
                         args.ko,
                         args.ko_hmm,
                         args.pfam,
                         args.tigrfam,
                         args.clusters,
                         args.orthologs,
                         args.cazy,
                         args.ec,
                         # Cutoffs
                         args.evalue,
                         args.bit,
                         args.id,
                         args.aln_query, 
                         args.aln_reference, 
                         args.c,
                         args.cut_ga,
                         args.cut_nc,
                         args.cut_tc,
                         args.cut_ko,
                         args.inflation,
                         args.chunk_number,
                         args.chunk_max,
                         args.count_domains,
                         # Parameters
                         args.threads,
                         args.parallel,
                         args.suffix,
                         args.light)
            a.do(args.genome_directory,
                 args.protein_directory, 
                 args.genome_files,
                 args.protein_files)

        
        elif args.subparser_name == self.CLASSIFY:
            self._check_classify(args)
            c = Classify()
            c.do(args.custom_modules, 
                 args.cutoff,
                 args.aggregate,
                 args.genome_and_annotation_file,
                 args.genome_and_annotation_matrix,
                 args.output)

        elif args.subparser_name == self.ENRICHMENT: 
            self._check_enrichment(args)
            e = Enrichment()
            e.do(# Input options
                 args.annotate_output,
                 args.metadata,
                 args.abundance,
                 args.abundance_metadata,
                 # Runtime options
                 args.genomes_to_compare_with_group,
                 args.pval_cutoff,
                 args.proportions_cutoff,
                 args.threshold,
                 args.multi_test_correction,
                 args.batchfile,
                 args.processes,
                 args.ko,
                 args.pfam,
                 args.tigrfam,
                 args.hypothetical,
                 args.cazy,
                 args.ec,
                 args.ko_hmm,
                 # Outputs
                 args.output)

        elif args.subparser_name == self.CONNECT:
            self._check_connect(args)
            c = Connect()
            c.do(args.annotate_output,
                 args.metadata,
                 args.custom_modules,
                 args.cutoff,
                 args.output)

        elif args.subparser_name == NetworkAnalyser.PATHWAY:
            self._check_network(args)
            na=NetworkAnalyser(args.metadata)
            na.do(args.matrix,
                  args.tpm_values,
                  args.abundance,
                  args.abundance_metadata,
                  args.metabolome,
                  args.enrichment_output,
                  args.depth,
                  args.filter,
                  args.limit,
                  args.queries,
                  args.starting_compounds, 
                  args.steps,
                  args.number_of_queries,
                  args.output)
        
        if args.subparser_name == self.PREDICT:
            self._check_predict(args)
            p = Predict()
            p.do(args.forester_model_directory,
                 args.input_matrix,
                 args.output)

        elif args.subparser_name == self.GENERATE:
            self._check_generate(args)
            gm = GenerateModel()
            gm.do(args.input_matrix,
                  args.groups,
                  args.model_type,
                  args.testing_portion,
                  args.grid_search,
                  args.threads,
                  args.output)
        
        logging.info('Done!')
