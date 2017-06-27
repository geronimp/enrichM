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
# Imports

import os

###############################################################################

DATABASE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'databases')

KO_DB = os.path.join(DATABASE_DIR, 'uniref100.dmnd')
PFAM_DB = os.path.join(DATABASE_DIR, 'pfam.hmm')
PFAM_CLAN_DB = os.path.join(DATABASE_DIR, 'pfam_clans.txt')
TIGRFAM_DB = os.path.join(DATABASE_DIR, 'tigrfam.hmm')
