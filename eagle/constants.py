import os
import configparser

from eaglib.constants import conf_constants as conf_constants_lib
from eagledb.constants import DB_INFO_NAME


CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(CONSTANTS_PATH, "tests")

# inner names (not configurable)
PROFILES_SCAN_OUT = "profiles_scan.hsr"
ORF_ALNS_DIR = "orf_alignments"
ORF_TREES_DIR = "orf_homolog_trees"


class ConfConstants(object):

    def __init__(self):
        self.num_threads = conf_constants_lib.num_threads
        self.mode = "genome"
        self.btax_det_method = "hmmer"
        self.min_orf_l = 300  # 100 amino acids

    def update_by_config(self, config_path):
        config = configparser.ConfigParser()
        config.read(config_path)

        self.num_threads = int(config.get(section="EAGLE", option="num_threads", fallback=self.num_threads))
        self.mode = config.get(section="EAGLE", option="mode", fallback=self.mode)
        self.btax_det_method = config.get(section="EAGLE", option="btax_det_method", fallback=self.btax_det_method)
        self.min_orf_l = int(config.get(section="EAGLE", option="min_orf_l", fallback=self.min_orf_l))


conf_constants = ConfConstants()
