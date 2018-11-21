import os
import shutil
import unittest

from EAGLE.constants import EAGLE_logger
from EAGLEdb.constants import TEST_DIR
from EAGLEdb.lib import db_creation


PACKAGE_DIR = 'lib'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestDBCreation(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
