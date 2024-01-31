import os
import configparser


CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(CONSTANTS_PATH, "tests")


class ConfConstants(object):

    def __init__(self):
        # GENERAL
        self.num_threads = 4
        # ALIGNMENT
        self.muscle_exec_path = "muscle"
        self.mafft_exec_path = "mafft"
        self.msaprobs_exec_path = "msaprobs"
        self.hmmalign_exec_path = "hmmalign"
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

    def update_by_config(self, config_path):
        config = configparser.ConfigParser()
        config.read(config_path)

        self.num_threads = int(config.get(section="EAGleLIB", option="num_threads", fallback=self.num_threads))
        self.muscle_exec_path = config.get(section="EAGleLIB", option="muscle_exec_path",
                                           fallback=self.muscle_exec_path)
        self.hmmalign_exec_path = config.get(section="EAGleLIB", option="hmmalign_exec_path",
                                           fallback=self.hmmalign_exec_path)
        self.mafft_exec_path = config.get(section="EAGleLIB", option="mafft_exec_path", fallback=self.mafft_exec_path)
        self.emboss_inst_dir = config.get(section="EAGleLIB", option="emboss_inst_dir", fallback=self.emboss_inst_dir)
        self.hmmer_inst_dir = config.get(section="EAGleLIB", option="hmmer_inst_dir", fallback=self.hmmer_inst_dir)
        self.infernal_inst_dir = config.get(section="EAGleLIB", option="infernal_inst_dir",
                                           fallback=self.infernal_inst_dir)
        self.blast_inst_dir = config.get(section="EAGleLIB", option="blast_inst_dir", fallback=self.blast_inst_dir)
        self.cons_thr = int(config.get(section="EAGleLIB", option="cons_thr", fallback=self.cons_thr))
        self.unif_window_l = int(config.get(section="EAGleLIB", option="unif_window_l", fallback=self.unif_window_l))
        self.unif_windows_step = int(config.get(section="EAGleLIB", option="unif_windows_step",
                                                fallback=self.unif_windows_step))
        self.kaks_calculator_exec_path = config.get(section="EAGleLIB", option="kaks_calculator_exec_path",
                                                    fallback=self.kaks_calculator_exec_path)
        self.kaks_window_l = int(config.get(section="EAGleLIB", option="kaks_window_l", fallback=self.kaks_window_l))
        self.kaks_top_fract = int(config.get(section="EAGleLIB", option="kaks_top_fract", fallback=self.kaks_top_fract))
        self.fastme_exec_path = config.get(section="EAGleLIB", option="fastme_exec_path",
                                           fallback=self.fastme_exec_path)


conf_constants = ConfConstants()
