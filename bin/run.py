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
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"

###############################################################################

import logging
import os
import shutil
from kegg_module_grabber import KeggModuleGrabber
from network_analyzer import NetworkAnalyser
from metagenome_analyzer import MetagenomeAnalyzer
from build_enrichment_matrix import BuildEncrichmentMatrix
from annotate import Annotate
from classifier import Classifier

###############################################################################

class Run:
<<<<<<< HEAD
    
    def __init__(self):
        self.ANNOTATE        = 'annotate'

        self.CLASSIFY        = 'classify'
        self.BUILD           = 'build'
        self.ENRICHMENT      = 'enrichment'
        self.MODULE_AB       = 'module_ab'

        self.PATHWAY         = 'pathway'
        self.EXPLORE         = 'explore'
    
=======

    def __init__(self):

        self.ANNOTATE        = 'annotate'
        self.MATRIX          = 'matrix'
        self.BUILD           = 'build'
        self.NETWORK         = 'network'
        self.EXPLORE         = 'explore'
        self.PATHWAY         = 'pathway'
        self.ENRICHMENT      = 'enrichment'
        self.MODULE_AB       = 'module_ab'
        self.TRAVERSE        = 'traverse'
        self.CLASSIFY        = 'classify'

>>>>>>> 33bdee9900b7c3a8d294ca7c47b434e1310dfcae
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
                logging.warning("Removing existing directory or file with name: %s" \
                                % args.output )
                shutil.rmtree(args.output)
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
            raise Exception("Input error: Either a list of genomes or a directory of genomes need to be specified.")
        
<<<<<<< HEAD
    def _check_annotate(self, args):
        pass
    
    def _check_enrichment(self, args):
        pass
=======
    def _check_classify(self, args):
        '''
        Check classify input and output options are valid.
        
        Parameters
        ----------
        args    - object. Argparse object
        '''
        # Ensure either an annotation matrix or list file has been specified:
        if not(args.genome_and_annotation_file or args.genome_and_annotation_matrix):
            raise Exception("Input error: An input file must be specified to either --genome_and_annotation_file or --genome_and_annotation_matrix")

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
>>>>>>> 33bdee9900b7c3a8d294ca7c47b434e1310dfcae

    def main(self, args):

        self._check_general(args)
        
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
                         args.threads)
            a.do(args.genome_directory, args.proteins_directory)

<<<<<<< HEAD
        if args.subparser_name == self.CLASSIFY:
            self._check_annotate(args)
            kem = KeggModuleGrabber()
            kem.do(args.custom_modules, 
                   args.output_prefix, 
                   args.genome_and_annotation_file, 
                   args.genome_and_annotation_matrix,
                   args.cutoff)

        elif args.subparser_name == self.ENRICHMENT:
            self._check_enrichment(args)
            bem = BuildEncrichmentMatrix()
            bem.do(args.annotations,
                   args.abundances,
                   args.metadata,
                   args.modules,
                   args.output_prefix)


        #elif args.subparser_name in self.network_options:
        #    na=NetworkAnalyser(args.metadata)
        #    na.main(args)
        #elif args.subparser_name in self.annotation_options:
        #    kmg = KeggModuleGrabber()
        #    kmg.main(args)
        #elif args.subparser_name in self.metagenome_annotation_options:
        #    ma = MetagenomeAnalyzer()
        #    ma.main(args)
        logging.info('Done')
=======
        elif args.subparser_name == self.CLASSIFY:
            self._check_classify(args)
            c = Classify()
            c.do(args.custom_modules, 
                 args.cutoff,
                 args.genome_and_annotation_file,
                 args.genome_and_annotation_matrix,
                 args.output)

        elif args.subparser_name == self.BUILD:
            bem = BuildEncrichmentMatrix()
            bem.do(args.annotations,
                     args.abundances, 
                     args.metadata, 
                     args.modules, 
                     args.output_prefix)

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
                  args.output_prefix)

        logging.info('Done!')
>>>>>>> 33bdee9900b7c3a8d294ca7c47b434e1310dfcae
