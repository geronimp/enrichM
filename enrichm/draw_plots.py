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

import logging
import os
from enrichm.databases import Databases
from enrichm.data import Data
from enrichm.toolbox import run_command
################################################################################
class Plot:

    def __init__(self):
        self.databases = Databases()
        path_to_scripts = os.path.split(os.path.realpath(__file__))[0]
        self.draw_pca_script_path = os.path.join(path_to_scripts, "PLOT_ko_pca.r")
        self.draw_heatmap_script_path = os.path.join(path_to_scripts, "PLOT_ko_heatmap.r")
        self.draw_barplots_script_path = os.path.join(path_to_scripts, "PLOT_ko_breakdown.r")
        self.ko00000 = os.path.join(Data.DATABASE_DIR, self.databases.DB_VERSION, 'ko00000.tsv')
        self.output_pca_plot = 'presence_absence_pca_plot.svg'
        self.output_heatmap_plot = 'presence_absence_pca_plot.svg'

    def draw_pca_plot(self, annotation_matrix, metadata, output_directory):
        logging.info('	- Generating PCA plot')
        output_path = os.path.join(output_directory, self.output_pca_plot)
        cmd = f"Rscript {self.draw_pca_script_path} \
                    -i {annotation_matrix} \
                    -m {metadata} \
                    -o {output_path} > /dev/null 2>&1"
        run_command(cmd)

    def draw_barplots(self, annotation_matrix, pvalue, output_directory):
        logging.info('	- Generating KO breakdown plots')
        cmd = f"Rscript {self.draw_barplots_script_path} \
                    -i {annotation_matrix} \
                    -o {output_directory} \
                    -k {self.ko00000} \
                    -p {pvalue} > /dev/null 2>&1"
        run_command(cmd)
