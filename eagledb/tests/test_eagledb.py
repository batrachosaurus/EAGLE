import os
import json
import shutil
import unittest

from eagle.constants import eagle_logger, conf_constants
from eagledb.constants import TEST_DIR, CONSTANTS_PATH, conf_constants_db, BTAX_JSON_NAME, BACTERIA_LIST_F_NAME
from eagledb.scheme import SeqProfileInfo
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
    test_genome_id = "GCF_000160075.2_ASM16007v2"
    conf_constants.fastme_exec_path = os.path.join(os.path.dirname(CONSTANTS_PATH), "docker_build", "fastme")
    conf_constants.msaprobs_exec_path = os.path.join(os.path.dirname(CONSTANTS_PATH), "docker_build", "msaprobs")
    conf_constants_db.btc_profile_aln_method = "MAFFT"
    conf_constants_db.k_max = 15
    conf_constants_db.k_min = 10

    def test_get_bacteria_from_ncbi(self, last_bact=30, use_prapared=True):
        if use_prapared:
            prepared_bacteria_f_path = os.path.join(INPUT_DIR, "prepared_bacteria.json")
        else:
            prepared_bacteria_f_path = ""
        genomes_list = bactdb_creation.get_bacteria_from_ncbi(bactdb_dir=OUTPUT_DIR,
                                                              last_bact=last_bact,
                                                              prepared_bacteria_f_path=prepared_bacteria_f_path)
        self.assertIs(type(genomes_list), list)
        return genomes_list

    def test_get_taxonomy(self):
        tax_f_name = self.test_genome_id + "_wgsmaster.gbff.gz"

        result_exp = (self.test_taxonomy, self.test_org_name)

        result = bactdb_creation.get_taxonomy(tax_f_name, INPUT_DIR, remove_tax_f=False)
        self.assertEqual(result, result_exp)

    def test_get_16S_fasta(self):
        rna_f_name = self.test_genome_id + "_rna_from_genomic.fna.gz"
        test_16S_fasta_result = self.test_genome_id + "_16S_rRNA.fasta"
        exp_seq_id_list = ["lcl|NZ_KI535341.1_rrna_33 [locus_tag=GCWU000182_RS08290] [product=16S ribosomal RNA] "
                           "[location=complement(179989..181551)] [gbkey=rRNA]",
                           "lcl|NZ_KI535342.1_rrna_59 [locus_tag=GCWU000182_RS09510] [product=16S ribosomal RNA] "
                           "[location=169454..>169911] [gbkey=rRNA]"]
        result = bactdb_creation.get_16S_fasta(rna_f_name, INPUT_DIR, genome_id=self.test_genome_id,
                                               remove_rna_f=False)
        result_f_name = os.path.split(result[0])[1].encode()
        shutil.move(result[0], os.path.join(OUTPUT_DIR, result_f_name))
        self.assertEqual(result_f_name, test_16S_fasta_result)
        self.assertEqual(frozenset(result[1]), frozenset(exp_seq_id_list))

    def test_get_btax_dict(self, use_test_results=True):
        if use_test_results:
            with open(os.path.join(OUTPUT_DIR, BACTERIA_LIST_F_NAME)) as genomes_list_f:
                genomes_list = json.load(genomes_list_f)
        else:
            genomes_list = self.test_get_bacteria_from_ncbi(last_bact=100, use_prapared=False)  # IMPORTANT!!!
        btc_profiles = [SeqProfileInfo(name="16S_rRNA", seq_type="nucl").get_json()]
        btax_dict = bactdb_creation.get_btax_dict(genomes_list,
                                                  btax_level=4,
                                                  btc_profiles=btc_profiles,
                                                  db_dir=OUTPUT_DIR,
                                                  build_tree=True,
                                                  num_threads=4,
                                                  save_alignments=True)
        with open(os.path.join(OUTPUT_DIR, BTAX_JSON_NAME), "w") as btax_dict_f:
            json.dump(btax_dict, btax_dict_f, indent=2)
        self.assertIsInstance(btax_dict, dict)
        return btax_dict

    def test_get_btax_blastdb(self, use_test_results=True):
        if use_test_results:
            with open(os.path.join(OUTPUT_DIR, BTAX_JSON_NAME)) as btax_dict_f:
                btax_dict = json.load(btax_dict_f)
        else:
            btax_dict = self.test_get_btax_dict()
        btax_dict = bactdb_creation.get_btax_blastdb(btax_dict=btax_dict, db_dir=OUTPUT_DIR)
        with open(os.path.join(OUTPUT_DIR, BTAX_JSON_NAME), "w") as btax_dict_f:
            json.dump(btax_dict, btax_dict_f, indent=2)
        self.assertIsInstance(btax_dict, dict)


class TestFilesUtils(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
