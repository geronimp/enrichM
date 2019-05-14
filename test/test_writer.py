import unittest
import filecmp
import tempfile
import os
import sys

sys.path = [os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '..', )]+sys.path

from enrichm.writer import Writer

class Tests(unittest.TestCase):

    def test_write(self):
        output_file = tempfile.mktemp()
        output_lines = [["this", "is", "a", "header"], ["Content", 1, 2, 3]]
        expected_list = ["this\tis\ta\theader\n", "Content\t1\t2\t3\n"]

        Writer.write(output_lines, output_file)

        for idx, line in enumerate(open(output_file)):
            self.assertEqual(line, expected_list[idx])
            
        self.assertEqual(idx+1, len(expected_list))
        os.remove(output_file)

if __name__ == "__main__":
    unittest.main()
