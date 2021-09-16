import os
import configparser

from eaglib.constants import conf_constants as conf_constants_lib


CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(CONSTANTS_PATH, "tests")

BACTERIA_LIST_F_NAME = "bacteria.json"
PREPARED_BACTERIA_F_NAME = "prepared_bacteria.json"
BACTERIA_GLOBAL_DIST_MATRIX = "global_dist_matr.phylip"
BACTERIA_SHORT_TO_FULL_ORG_NAMES = "short_to_full_org_names.json"
BACT_FAM_F_NAME = "bact_fam.json"

BTAX_JSON_NAME = "btax.json"

ORG_TABLES_DIR = os.path.join(CONSTANTS_PATH, "org_tables")
REFSEQ_BACTERIA_TABLE = os.path.join(ORG_TABLES_DIR, "refseq_bacteria_table.txt")###
BACTERIA_REFSEQ_SUMMARY_LINK = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"###
GENBANK_BACTERIA_TABLE = os.path.join(ORG_TABLES_DIR, "genbank_bacteria_table.txt")
BACTERIA_GENBANK_SUMMARY_LINK = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"###
DEFAULT_BACTDB_DIR = os.path.join("EAGLEdb", "bacteria")###

PROFILES_DB_NAME = "db_repr_profiles"  # it is inner name (not configurable)
DB_INFO_NAME = "db_info.json"


class ConfConstants(object):

    def __init__(self):
        # GENERAL
        self.num_threads = conf_constants_lib.num_threads
        # Bacteria db
        self.only_repr = False
        self.btax_level = 4
        self.k_max = 30
        self.k_min = 20
        self.btc_profile_aln_method = "MAFFT"

    def update_by_config(self, config_path):
        config = configparser.ConfigParser()
        config.read(config_path)

        self.num_threads = int(config.get(section="EAGLEdb", option="num_threads", fallback=self.num_threads))
        self.only_repr = bool(int(config.get(section="EAGLEdb", option="only_repr", fallback=self.only_repr)))
        self.btax_level = int(config.get(section="EAGLEdb", option="btax_level", fallback=self.btax_level))
        self.k_max = int(config.get(section="EAGLEdb", option="k_max", fallback=self.k_max))
        self.k_min = int(config.get(section="EAGLEdb", option="k_min", fallback=self.k_min))
        self.btc_profile_aln_method = config.get(section="EAGLEdb", option="btc_profile_aln_method",
                                                 fallback=self.btc_profile_aln_method)


conf_constants = ConfConstants()
