import os
import json
import shutil
import unittest

from eagle.constants import TEST_DIR, eagle_logger, CONSTANTS_PATH, conf_constants
from eagledb.scheme import DBInfo
from eagle import orfs_explorer, btax_scanner


PACKAGE_DIR = 'eagle'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

EAGLEDB_TEST_DIR = os.path.join(os.path.dirname(CONSTANTS_PATH), "eagledb", "tests")
EAGLEDB_TEST_DATA = os.path.join(EAGLEDB_TEST_DIR, "test_data")

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestORFsExplorer(unittest.TestCase):

    conf_constants.fastme_exec_path = os.path.join(os.path.dirname(CONSTANTS_PATH), "docker_build", "fastme")
    conf_constants.msaprobs_exec_path = os.path.join(os.path.dirname(CONSTANTS_PATH), "docker_build", "msaprobs")
    conf_constants.kaks_calculator_exec_path = \
        os.path.join(os.path.dirname(CONSTANTS_PATH), "docker_build", "KaKs_Calculator")

    db_info_dict = DBInfo().get_json()
    with open(os.path.join(EAGLEDB_TEST_DATA, "eagledb", "db_info.json")) as db_info_f:
        for k, v in json.load(db_info_f).items():
            if type(v) in (str, unicode):
                db_info_dict[k] = os.path.join(EAGLEDB_TEST_DIR, v)

    def test_explore_orfs(self):
        orfs_explorer.explore_orfs(
            in_fasta=os.path.join(INPUT_DIR, "NC_000913.fasta"),  # Escherichia coli K-12 MG1655
            db_json=self.db_info_dict,
            out_dir=OUTPUT_DIR,
            btax_name="Enterobacterales_related",
            num_threads=4,
            tblastn_result_path=os.path.join(OUTPUT_DIR, "NC_000913.fasta.bl"),
            save_alignments=True,
            save_trees=True,
        )
        self.assertTrue(True)


class TestBTaxScanner(unittest.TestCase):

    db_info_dict = TestORFsExplorer.db_info_dict
    with open(db_info_dict["btax_json"]) as btax_dict_f:
        btax_dict = json.load(btax_dict_f)

    def test_get_btax_name(self):
        exp_btax_name = "Enterobacterales_related"
        btax_name = btax_scanner.get_btax_name(
            in_fasta=os.path.join(INPUT_DIR, "NC_000913.fasta"),  # Escherichia coli K-12 MG1655
            profiles_db=self.db_info_dict["repr_profiles"],
            btax_names=self.btax_dict.keys(),
            work_dir=OUTPUT_DIR,
            num_threads=4
            )
        eagle_logger.info("got base taxon name '%s'", btax_name)
        self.assertEqual(exp_btax_name, btax_name)


if __name__ == "__main__":
    unittest.main()
