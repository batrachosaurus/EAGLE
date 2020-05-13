import os
import shutil

from eagle.constants import conf_constants
from eagle.lib.general import generate_random_string


class SeqsProfileHandler(object):

    def __init__(self,
                 hmmer_inst_dir=None,
                 infernal_inst_dir=None,
                 tmp_dir=None,
                 config_path=None,
                 logger=None):

        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_hmm_tmp"

        self.tmp_dir = tmp_dir
        self.hmmer_inst_dir = hmmer_inst_dir
        self.infernal_inst_dir = infernal_inst_dir
        self.logger = logger
