import subprocess
import os
import sys

class Convert:
	"""docstring for parse"""
	def do(self, genome_list):
		for genome in genome_list:
			g = os.path.basename(genome)
			cmd = """sed 's/>/>%s/g' %s """ % (g, genome)
			subprocess.call(cmd, shell=True)
		return cmd

Convert().do(sys.argv[1:])