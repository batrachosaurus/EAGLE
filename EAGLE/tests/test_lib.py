import os
import shutil
import unittest

from EAGLE.constants import TEST_DIR, EAGLE_logger
from EAGLE.lib import general, alignment, phylo, seqs


PACKAGE_DIR = 'lib'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestGeneral(unittest.TestCase):

    def test_compare_files(self):
        self.assertTrue(general.compare_files(os.path.join(INPUT_DIR, "f1.txt"), os.path.join(INPUT_DIR, "f1.txt")))
        self.assertFalse(general.compare_files(os.path.join(INPUT_DIR, "f1.txt"), os.path.join(INPUT_DIR, "f2.txt")))


class TestAlignment(unittest.TestCase):
    pass


class TestPhylo(unittest.TestCase):
    pass


class TestSeqs(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()