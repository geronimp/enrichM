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
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

###############################################################################
# Imports

import logging
import sys
import os
import shutil
import time

from data import Data
from network_analyzer import NetworkAnalyser
from enrichment import Enrichment
from annotate import Annotate
from classifier import Classify
from comparer import Compare

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

    def _logging_setup(self, args):
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

    def _check_data(self, args):
        '''
        Check data input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object
        '''
        pass

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
        if(args.genome_directory or args.genome_files):
            args.suffix = '.fna'
        elif(args.protein_directory or args.protein_files):
            args.suffix = '.faa'
        if(args.id>1 or args.id<0):
            raise Exception("Identity (--id) must be between 0 and 1.")

    def _check_enrichment(self, args):
        '''
        Check enrichment input and output options are valid.

        Parameters
        ----------
        args    - object. Argparse object

        Output
        ------
        '''
        if not any([args.do_all, args.do_gvg, args.do_ivg, args.do_ivi]):
            raise Exception("Input error: No comparisons were specified. You will need to tell enrichM \
which statistical tests to run using the --do_ivi --do_gvg --do_ivg, or --do_all flags")
        ### ~ TODO: Check Multi test correction inputs...
        if not(args.annotation_matrix or args.annotation_file):
            raise Exception("Input error: No input file was specified. Please specify annotations to either the --annotation_matrix --annotation_file flags")
    
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

    def _check_build(self, args):
        '''
        Check build input and output options are valid.
        
        Parameters
        ----------
        args    - object. Argparse object
        '''
        if not(args.metadata or args.abundances):
            raise Exception("No metadata or abundance information provided to build.")
    
    def _check_compare(self, args):
        pass ### ~ TODO: Dunno yet
    
    def _check_network(self, args):
        '''
        Check network (explore, pathway) input and output options are valid.
        
        Parameters
        ----------
        args    - object. Argparse object
        '''

        if args.subparser_name==NetworkAnalyser.EXPLORE:
            args.filter             = None  
            args.limit              = None 
            args.starting_compounds = None
            args.steps              = None
            args.number_of_queries  = None

            if not(args.queries):
                if args.depth:
                    logging.warning("--depth argument ignored without --queries flag")

        if args.subparser_name==NetworkAnalyser.PATHWAY:
            args.depth              = None
            args.queries            = None
            args.starting_compounds = None
            args.steps              = None
            args.number_of_queries  = None
        
        if args.subparser_name==NetworkAnalyser.TRAVERSE:
            args.depth              = None
            args.queries            = None


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
            self._check_data(args)
            d = Data()

            d.do()
        
        if args.subparser_name == self.ANNOTATE:
            self._check_annotate(args)
            a = Annotate(# Define inputs and outputs
                         args.output,
                         # Define type of annotation to be carried out
                         args.ko,
                         args.pfam,
                         args.tigrfam,
                         args.cog,
                         args.hypothetical,
                         # Cutoffs
                         args.evalue,
                         args.bit,
                         args.id,
                         args.aln_query, 
                         args.aln_reference, 
                         args.cascaded,
                         args.c,
                         # Parameters
                         args.threads,
                         args.parallel,
                         args.suffix)

            a.do(args.genome_directory, args.protein_directory, 
                 args.genome_files, args.protein_files)

        elif args.subparser_name == self.CLASSIFY:
            self._check_classify(args)
            c = Classify()
            c.do(args.custom_modules, 
                 args.cutoff,
                 args.genome_and_annotation_file,
                 args.genome_and_annotation_matrix,
                 args.output)

        elif args.subparser_name == self.ENRICHMENT: 
            self._check_enrichment(args)
            e = Enrichment()
            e.do(# Inputs
                 args.annotation_matrix,
                 args.annotation_file,
                 args.metadata,
                 args.modules,
                 args.abundances,
                 args.do_all,
                 args.do_ivi, 
                 args.do_gvg,
                 args.do_ivg,
                 args.pval_cutoff,
                 args.proportions_cutoff,
                 args.threshold,
                 args.multi_test_correction,
                 args.output)

        elif args.subparser_name == self.COMPARE:
            self._check_compare(args)
            c = Compare()
            c.do(args.enrichm_annotate_output)

        elif(args.subparser_name == NetworkAnalyser.PATHWAY or
             args.subparser_name == NetworkAnalyser.EXPLORE or
             args.subparser_name == NetworkAnalyser.TRAVERSE):

            self._check_network(args)
            na=NetworkAnalyser(args.metadata)
            na.do(args.matrix,
                  args.transcriptome,
                  args.metabolome,
                  args.depth,
                  args.filter,
                  args.limit,
                  args.queries,
                  args.subparser_name,
                  args.starting_compounds, 
                  args.steps,
                  args.number_of_queries,
                  args.output)
        logging.info('Done!')
