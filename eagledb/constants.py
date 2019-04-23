import os

from eagle.constants import DEFAULT_CONFIG
from eagle.lib.general import ConfBase

CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(CONSTANTS_PATH, "tests")

BACTERIA_LIST_F_NAME = "bacteria.json"
PREPARED_BACTERIA_F_NAME = "prepared_bacteria.json"
BACTERIA_GLOBAL_DIST_MATRIX = "global_dist_matr.phylip"
BACTERIA_SHORT_TO_FULL_ORG_NAMES = "short_to_full_org_names.json"
BACT_FAM_F_NAME = "bact_fam.json"

BTAX_JSON_NAME = "btax.json"
DB_INFO_NAME = "db_info.json"

ORG_TABLES_DIR = os.path.join(CONSTANTS_PATH, "org_tables")
DEFAULT_REFSEQ_BACTERIA_TABLE = os.path.join(ORG_TABLES_DIR, "refseq_bacteria_table.txt")
DEFAULT_GENBANK_BACTERIA_TABLE = os.path.join(ORG_TABLES_DIR, "genbank_bacteria_table.txt")
DEFAULT_BACTDB_DIR = os.path.join("EAGLEdb", "bacteria")

PROFILES_DB_NAME = "db_repr_profiles"  # it is inner name (not configurable)


class ConfConstants(ConfBase):

    def __init__(self, config_path=DEFAULT_CONFIG):
        # Bacteria db
        self.only_repr = False
        self.btax_level = 3
        self.btc_profile_aln_method = "MSAProbs"

        super(ConfConstants, self).__init__(config_path=config_path)


conf_constants_db = ConfConstants()
