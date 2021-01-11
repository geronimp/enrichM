#!/usr/bin/env python
# Imports
import unittest
import os.path
import sys
import subprocess

###############################################################################

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','enrichm')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

###############################################################################

class Tests(unittest.TestCase):

    def test_hello_world(self):
        cmd = '%s -h > /dev/null' % path_to_script
        subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    unittest.main()
