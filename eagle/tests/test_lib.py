import os
import unittest

import eaglib._utils.files
from eagle.constants import TEST_DIR
from eagle.lib import general

PACKAGE_DIR = 'eaglib'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestGeneral(unittest.TestCase):

    def test_compare_files(self):
        self.assertTrue(
            eaglib._utils.files.compare_files(os.path.join(INPUT_DIR, "f1.txt"), os.path.join(INPUT_DIR, "f1.txt")))
        self.assertFalse(
            eaglib._utils.files.compare_files(os.path.join(INPUT_DIR, "f1.txt"), os.path.join(INPUT_DIR, "f2.txt")))


class TestAlignment(unittest.TestCase):
    pass


class TestPhylo(unittest.TestCase):
    pass


class TestSeqs(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()