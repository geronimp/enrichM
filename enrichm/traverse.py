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
 
###############################################################################

class NetworkTraverser:
    
    to_omit = set(["C00828",  # Menaquinone
                   "C00534",  # Pyridoxamine
                   "C00006",  # NADP+
                   "C00003",  # NAD+
                   "C00002",  # ATP
                   "C00314",  # Pyridoxine
                   "C00864",  # Pantothenate
                   "C00504",  # Folate
                   "C00032",  # Heme
                   "C05443",  # Vitamin D3
                   "C00253",  # Nicotinate
                   "C00250",  # Pyridoxal
                   "C11378",  # Ubiquinone-10
                   "C05777",  # Coenzyme F430
                   "C00072",  # Ascorbate
                   "C00378",  # Thiamine
                   "C00101",  # Tetrahydrofolate
                   "C00029",  # UDP-glucose
                   "C00068",  # Thiamin diphosphate
                   "C00061",  # FMN
                   "C00063",  # CTP
                   "C05776",  # Vitamin B12
                   "C00113",  # PQQ
                   "C18237",  # Molybdoenzyme molybdenum cofactor
                   "C00051",  # Glutathione
                   "C00010",  # CoA
                   "C00016",  # FAD
                   "C00018",  # Pyridoxal phosphate
                   "C00019",  # S-Adenosyl-L-methionine
                   "C00153",  # Nicotinamide
                   "C04628",  # Coenzyme B
                   "C00862",  # Methanofuran
                   "C15672",  # Heme O
                   "C15670",  # Heme A
                   "C02059",  # Phylloquinone
                   "C03576",  # Coenzyme M
                   "C05441",  # Vitamin D2
                   "C00272",  # Tetrahydrobiopterin
                   "C02477",  # alpha-Tocopherol
                   "C00473",  # Retinol
                   "C00120",  # Biotin
                   "C00725",  # Lipoate
                   "C00053",  # 3'-Phosphoadenylyl sulfate
                   "C00194",  # Cobamide coenzyme
                   "C00255",  # Riboflavin
                   'C00001',  # H2O
                   'C00008',  # ADP
                   'C00013',  # Diphosphate
                   'C00004',  # NADH
                   'C00005',  # NADPH
                   'C00080',  # H+
                   'C00009',  # Orthophosphate
                   'C00008',  # ADP
                   'C00004',  # NADH
                   'C00020',  # AMP
                   'C00007',  # Oxygen
                   'C00015']) # UDP

    @staticmethod
    def dijkstra_traverse(from_node, to_node):
        pass

    @staticmethod  
    def shortest_bfs_path(r2rpair_graph, c2r_graph, from_node, to_node):
    # Modified from:
    # http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
        def bfs_traverse(r2rpair_graph, c2r_graph, from_node, to_node):
            queue = [(from_node, [from_node])]
            while queue:
                (vertex, path) = queue.pop(0)     
                for reaction in c2r_graph[vertex]:
                    for rpair in r2rpair_graph[reaction]:
                        if vertex in rpair:
                            compound = [x for x in rpair.split('_') 
                                          if x != vertex].pop()
                            new_path = path + [reaction, compound]
                            if compound==to_node:
                                yield path + [reaction, compound]
                            else:
                                queue.append((compound, new_path))

        try:
            return next(bfs_traverse(r2rpair_graph, c2r_graph, from_node, 
                                     to_node))
        except StopIteration:
            return None