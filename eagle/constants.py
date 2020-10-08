import os
import configparser


CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(CONSTANTS_PATH, "tests")

# inner names (not configurable)
DB_INFO_NAME = "db_info.json"
PROFILES_SCAN_OUT = "profiles_scan.hsr"
ORF_ALNS_DIR = "orf_alignments"
ORF_TREES_DIR = "orf_homolog_trees"


class ConfConstants(object):

    def __init__(self):
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
        self.infernal_inst_dir = ""
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
        self.btax_det_method = "hmmer"
        self.min_orf_l = 300  # 100 amino acids

    def update_by_config(self, config_path):
        pass


conf_constants = ConfConstants()
