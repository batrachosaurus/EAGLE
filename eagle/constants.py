import logging
import os

from eagle.lib.general import ConfBase, get_redis_server, setup_logging

CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(CONSTANTS_PATH, "tests")
CONF_DIR_NAME = 'configs'
CONF_DIR_PATH = os.path.join(CONSTANTS_PATH, CONF_DIR_NAME)
DEFAULT_CONFIG = os.path.join(CONF_DIR_PATH, "default_config.ini")
LOG_CONFIG_NAME = 'log_conf.yaml'
LOGGER_NAME = 'EAGLE_logger'

# inner names (not configurable)
PROFILES_SCAN_OUT = "profiles_scan.out"
ORF_ALNS_DIR = "orf_alignments"
ORF_TREES_DIR = "orf_homolog_trees"

setup_logging(os.path.join(CONF_DIR_PATH, LOG_CONFIG_NAME))
eagle_logger = logging.getLogger(LOGGER_NAME)


class ConfConstants(ConfBase):

    def __init__(self, config_path=DEFAULT_CONFIG):
        # GENERAL
        self.num_threads = 4
        self.redis_host = 'localhost'
        self.redis_port = 6379
        self.redis_queue_db = 0
        # ALIGNMENT
        self.muscle_exec_path = "muscle"
        self.mafft_exec_path = "mafft"
        self.msaprobs_exec_path = "msaprobs"
        self.emboss_inst_dir = ""
        self.hmmer_inst_dir = ""
        self.blast_inst_dir = ""
        self.cons_thr = 0.98
        self.unif_window_l= 10
        self.unif_windows_step = 5
        # Ka/Ks
        self.kaks_calculator_exec_path = "KaKs_Calculator"
        self.kaks_window_l = 150  # 50 codons
        self.kaks_top_fract = 0.5
        # PHYLO
        self.fastme_exec_path = "fastme"
        # eagle
        self.mode = "genome"
        self.min_orf_l = 300  # 100 amino acids

        super(ConfConstants, self).__init__(config_path=config_path)

    def get_redis_server(self, restart=True):
        eagle_logger.info(get_redis_server(self.redis_host, self.redis_port, restart=restart))


conf_constants = ConfConstants()
