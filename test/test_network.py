import unittest
import filecmp
import tempfile
import os
import sys

sys.path = [os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '..', )]+sys.path

from enrichm.network_analyzer import NetworkAnalyser

class Tests(unittest.TestCase):

    def test_hello_pathway(self):
        # network_analyser = NetworkAnalyser()
        pass
    def test_hello_query(self):
        # network_analyser = NetworkAnalyser()
        pass

if __name__ == "__main__":
    unittest.main()
