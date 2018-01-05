import os
import shutil
import unittest

from EAGLEdb.constants import constants_path
from EAGLEdb.bactdb import get_bacteria_from_ncbi
from EAGLEdb.bactdb.bactdb_creator import get_taxonomy, get_16S_fasta

class TestBactDBCreator(unittest.TestCase):

    test_dir = os.path.join(constants_path, "bactdb", "test")
    input_dir = os.path.join(test_dir, "test_data")
    output_dir = os.path.join(test_dir, "test_results")
    tax_f_name = "GCF_000160075.2_ASM16007v2_wgsmaster.gbff.gz"
    test_tax_result = {"family": u'Aerococcaceae',
                       "genus": u'Abiotrophia',
                       "species": u'Abiotrophia_defectiva',
                       "strain": u'Abiotrophia_defectiva_ATCC_49176'}
    rna_f_name = "GCF_000160075.2_ASM16007v2_rna_from_genomic.fna.gz"
    test_16S_fasta_result = "Abiotrophia_defectiva_ATCC_49176.fasta"
    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    def test_get_bacteria_from_ncbi(self):
        result = get_bacteria_from_ncbi(bactdb_dir=self.output_dir, n_first_bact=20)
        print len(result)
        print result
        self.assertIs(type(result), list)

    def test_get_taxonomy(self):
        result = {"family": None,
                  "genus": None,
                  "species": None,
                  "strain": None}
        result["family"], result["genus"], result["species"], result["strain"] = \
            get_taxonomy(self.tax_f_name, self.input_dir, remove_tax_f=False)
        print result
        self.assertEqual(result, self.test_tax_result)

    def test_get_16S_fasta(self):
        result = get_16S_fasta(self.rna_f_name, self.input_dir, strain=self.test_tax_result["strain"],
                               remove_rna_f=False)
        result_f_name = os.path.split(result)[1].encode()
        shutil.move(result, os.path.join(self.output_dir, result_f_name))
        print result_f_name
        self.assertEqual(result_f_name, self.test_16S_fasta_result)


if __name__ == "__main__":
    unittest.main()
