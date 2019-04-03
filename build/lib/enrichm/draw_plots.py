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
__copyright__   = "Copyright 2017"
__credits__     = ["Joel Boyd"]
__license__     = "GPL3"
__version__     = "0.0.7"
__maintainer__  = "Joel Boyd"
__email__       = "joel.boyd near uq.net.au"
__status__      = "Development"

################################################################################
import logging
import os
import subprocess
from enrichm.databases import Databases
from enrichm.data import Data
################################################################################
class Plot:

	def __init__(self):
		path_to_scripts = os.path.split(os.path.realpath(__file__))[0]
		self.draw_pca_script_path 		= os.path.join(path_to_scripts, "PLOT_ko_pca.r") 
		self.draw_heatmap_script_path 	= os.path.join(path_to_scripts, "PLOT_ko_heatmap.r") 
		self.draw_barplots_script_path 	= os.path.join(path_to_scripts, "PLOT_ko_breakdown.r") 
		self.ko00000 					= os.path.join(Data.DATABASE_DIR, Databases.DB_VERSION, 'ko00000.tsv')
		self.output_pca_plot = 'presence_absence_pca_plot.svg'
		self.output_heatmap_plot = 'presence_absence_pca_plot.svg'

	def draw_pca_plot(self, annotation_matrix, metadata, output_directory):
		logging.info('	- Generating PCA plot')
		output_path = os.path.join(output_directory, self.output_pca_plot)
		
		cmd = "Rscript %s -i %s -m %s -o %s > /dev/null 2>&1" \
			% (self.draw_pca_script_path, annotation_matrix, metadata, output_path)
		subprocess.call(cmd, shell=True)

	def draw_barplots(self, annotation_matrix, pvalue, output_directory):
		logging.info('	- Generating KO breakdown plots')
		cmd = "Rscript %s -i %s -o %s -k %s -p %f > /dev/null 2>&1" \
			% (self.draw_barplots_script_path, annotation_matrix, output_directory, self.ko00000, pvalue)
		subprocess.call(cmd, shell=True)
