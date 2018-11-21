import os
import shutil
import unittest

from EAGLE.constants import EAGLE_logger
from Research.constants import TEST_DIR
from Research import check_bacteria_response


PACKAGE_DIR = 'Research'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestCheckBacteriaResponse(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
