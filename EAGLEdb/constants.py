import os

from EAGLE.constants import DEFAULT_CONFIG
from EAGLE.lib.general import ConfBase

CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(CONSTANTS_PATH, "tests")

BACTERIA_LIST_F_NAME = "bacteria.json"
ANALYZED_BACTERIA_F_NAME = "analyzed_bacteria.json"
BACT_FAM_F_NAME = "bact_fam.json"

ORG_TABLES_DIR = os.path.join(CONSTANTS_PATH, "org_tables")
DEFAULT_REFSEQ_BACTERIA_TABLE = os.path.join(ORG_TABLES_DIR, "refseq_bacteria_table.txt")
DEFAULT_GENBANK_BACTERIA_TABLE = os.path.join(ORG_TABLES_DIR, "genbank_bacteria_table.txt")
DEFAULT_BACTDB_DIR = os.path.join("EAGLEdb", "bacteria")

PROFILES_DB_NAME = "db_repr_profiles"  # it is inner name (not configurable)


class ConfConstants(ConfBase):

    def __init__(self, config_path=DEFAULT_CONFIG):
        # Bacteria db
        self.only_repr=False

        super(ConfConstants, self).__init__(config_path=config_path)


conf_constants_db = ConfConstants()
