import os
import shutil
import unittest

from EAGLE.constants import EAGLE_logger
from EAGLEdb.constants import TEST_DIR
from EAGLEdb import bactdb_creating, files_utils


PACKAGE_DIR = 'EAGLEdb'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestBactDBCreating(unittest.TestCase):

    test_tax = {"family": u'Aerococcaceae',
                       "genus": u'Abiotrophia',
                       "species": u'Abiotrophia_defectiva',
                       "strain": u'Abiotrophia_defectiva_ATCC_49176'}

    def test_get_bacteria_from_ncbi(self):
        analyzed_bacteria_f_path = os.path.join(INPUT_DIR, "analyzed_bacteria.json")
        result = bactdb_creating.get_bacteria_from_ncbi(bactdb_dir=OUTPUT_DIR,
                                                        last_bact=20,
                                                        analyzed_bacteria_f_path=analyzed_bacteria_f_path)
        self.assertIs(type(result), list)

    def test_get_taxonomy(self):
        result = {
            "family": None,
            "genus": None,
            "species": None,
            "strain": None,
        }
        tax_f_name = "GCF_000160075.2_ASM16007v2_wgsmaster.gbff.gz"

        result["family"], result["genus"], result["species"], result["strain"] = \
            bactdb_creating.get_taxonomy(tax_f_name, INPUT_DIR, remove_tax_f=False)
        self.assertEqual(result, self.test_tax)

    def test_get_16S_fasta(self):
        rna_f_name = "GCF_000160075.2_ASM16007v2_rna_from_genomic.fna.gz"
        test_16S_fasta_result = "Abiotrophia_defectiva_ATCC_49176_16S_rRNA.fasta"
        result = bactdb_creating.get_16S_fasta(rna_f_name, INPUT_DIR, strain=self.test_tax["strain"],
                                               remove_rna_f=False)
        result_f_name = os.path.split(result)[1].encode()
        shutil.move(result, os.path.join(OUTPUT_DIR, result_f_name))
        self.assertEqual(result_f_name, test_16S_fasta_result)


class TestFilesUtils(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
