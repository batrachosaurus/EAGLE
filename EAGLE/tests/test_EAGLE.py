import os
import shutil
import unittest

from EAGLE.constants import TEST_DIR, EAGLE_logger
from EAGLE import eag_location_explorer


PACKAGE_DIR = 'EAGLE'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestEAGLocationExplorer(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
