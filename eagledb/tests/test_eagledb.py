import os
import shutil
import unittest

from eagle.constants import eagle_logger
from eagledb.constants import TEST_DIR
from eagledb import bactdb_creation, files_utils


PACKAGE_DIR = 'eagledb'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestBactDBCreation(unittest.TestCase):

    test_taxonomy = ["Bacteria", "Firmicutes", "Bacilli", "Lactobacillales",
                     "Aerococcaceae", "Abiotrophia", "Abiotrophia_defectiva"]
    test_org_name = "Abiotrophia_defectiva_ATCC_49176"

    def test_get_bacteria_from_ncbi(self):
        analyzed_bacteria_f_path = os.path.join(INPUT_DIR, "prepared_bacteria.json")
        result = bactdb_creation.get_bacteria_from_ncbi(bactdb_dir=OUTPUT_DIR,
                                                        last_bact=30,
                                                        prepared_bacteria_f_path=analyzed_bacteria_f_path)
        self.assertIs(type(result), list)

    def test_get_taxonomy(self):
        tax_f_name = "GCF_000160075.2_ASM16007v2_wgsmaster.gbff.gz"

        result_exp = (self.test_taxonomy, self.test_org_name)

        result = bactdb_creation.get_taxonomy(tax_f_name, INPUT_DIR, remove_tax_f=False)
        self.assertEqual(result, result_exp)

    def test_get_16S_fasta(self):
        rna_f_name = "GCF_000160075.2_ASM16007v2_rna_from_genomic.fna.gz"
        test_16S_fasta_result = "Abiotrophia_defectiva_ATCC_49176_16S_rRNA.fasta"
        result = bactdb_creation.get_16S_fasta(rna_f_name, INPUT_DIR, strain=self.test_org_name,
                                               remove_rna_f=False)
        result_f_name = os.path.split(result)[1].encode()
        shutil.move(result, os.path.join(OUTPUT_DIR, result_f_name))
        self.assertEqual(result_f_name, test_16S_fasta_result)


class TestFilesUtils(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
