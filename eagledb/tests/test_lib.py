import os
import json
import shutil
import unittest

from eagle.constants import eagle_logger
from eagledb.constants import TEST_DIR, PROFILES_DB_NAME, BTAX_JSON_NAME
from eagledb.lib import db_creation


PACKAGE_DIR = 'lib'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestDBCreation(unittest.TestCase):

    db_dir = os.path.join(TEST_DIR, 'test_results', 'eagledb')
    with open(os.path.join(db_dir, BTAX_JSON_NAME)) as btax_json_f:
        btax_dict = json.load(btax_json_f)

    def test_create_profiles_db(self):
        exp_repr_profiles_path = os.path.join(self.db_dir, PROFILES_DB_NAME)
        repr_profiles_path = db_creation.create_profiles_db(btax_dict=self.btax_dict,
                                                            db_dir=self.db_dir,
                                                            profiles_db_name=PROFILES_DB_NAME,
                                                            method="hmmer",
                                                            logger=eagle_logger)
        self.assertEqual(exp_repr_profiles_path, repr_profiles_path)


if __name__ == "__main__":
    unittest.main()
