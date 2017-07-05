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

import logging
import sys
import os
import shutil
from kegg_module_grabber import KeggModuleGrabber
from network_analyzer import NetworkAnalyser
from metagenome_analyzer import MetagenomeAnalyzer
from enrichment import Enrichment
from annotate import Annotate
from classifier import Classify

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

        self.CLASSIFY        = 'classify'
        self.BUILD           = 'build'
        self.ENRICHMENT      = 'enrichment'
        self.MODULE_AB       = 'module_ab'

        self.PATHWAY         = 'pathway'
        self.EXPLORE         = 'explore'

    def _logging_setup(self, args):

        logger = logging.getLogger('')
        logger.setLevel(debug[args.verbosity])
        log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                       datefmt="%Y-%m-%d %H:%M:%S %p")

        stream_logger = logging.StreamHandler(sys.stdout)
        stream_logger.setFormatter(log_format)
        stream_logger.setLevel(debug[args.verbosity])
        logger.addHandler(stream_logger)

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
        # Set up working directory
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
        # ensure either a list of genomes or a directory of genomes have been specified
        if not(args.genome_files or args.genome_directory or args.proteins_directory):
            raise Exception("Input error: Either a list of genomes or a directory of genomes \
need to be specified.")
    
    def _check_enrichment(self, args):
        '''
        
        Parameters
        ----------
        
        Output
        ------
        '''
        if not(args.enrichm_annotations or args.annotation_matrix): 
            raise Exception("Input error: Either enrichm annotations (--enrichm_annotations) or an \
annotation matrix (--annotation_matrix) need to be specified.")
        if(args.enrichm_annotations and args.annotation_matrix): 
            raise Exception("Input error: Only one set of comparisons can be made at any time, so \
please specify either enrichM annotations (--enrichm_annotations) OR an annotation matrix \
(--annotation_matrix)")
        ### ~ TODO: Check Multi test correction inputs...
        
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
            
    def _check_network(self, input, output):
        '''
        Check network (explore, pathway) input and output options are valid.
        
        Parameters
        ----------
        args    - object. Argparse object
        '''
        if args.subparser_name==NetworkAnalyser.EXPLORE:
          if not(args.queries):
            if args.depth:
              logging.warning("--depth argument ignored without --queries flag")

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
        if args.subparser_name == self.ANNOTATE:
            self._check_annotate(args)
            a = Annotate(# Define inputs and outputs
                         args.genome_files,
                         args.output,
                         # Define type of annotation to be carried out
                         args.ko,
                         args.pfam,
                         args.tigrfam,
                         args.cog,
                         # Cutoffs
                         args.evalue,
                         args.bit,
                         args.id,
                         args.aln_query, 
                         args.aln_reference, 
                         # Parameters
                         args.threads
                         )
            a.do(args.genome_directory, args.proteins_directory)

        elif args.subparser_name == self.CLASSIFY:
            self._check_classify(args)
            c = Classify()
            c.do(args.custom_modules, 
                 args.cutoff,
                 args.genome_and_annotation_file,
                 args.genome_and_annotation_matrix,
                 args.output
                 )

        elif args.subparser_name == self.ENRICHMENT: 
            ### ~ TODO: Make sure an output directory is written, rather than to files.
            self._check_enrichment(args)
            e = Enrichment()
            e.do(# Inputs
                 args.enrichm_annotations,
                 args.annotation_matrix,
                 args.metadata,
                 args.modules,
                 args.abundances,
                 args.no_ivi, 
                 args.no_gvg,
                 args.no_ivg,
                 args.cutoff,
                 args.threshold,
                 args.multi_test_correction,
                 args.output
                 )

        elif(args.subparser_name == self.PATHWAY or args.subparser_name == self.EXPLORE):
            self._check_network(args)
            na=NetworkAnalyser(args.metadata)
            na.do(args.depth,
                  args.filter,
                  args.limit,
                  args.metabolome,
                  args.number_of_queries,
                  args.queries,
                  args.starting_compounds,
                  args.steps,
                  args.subparser_name,
                  args.transcriptome,
                  args.output_prefix
                  )

        logging.info('Done!')
