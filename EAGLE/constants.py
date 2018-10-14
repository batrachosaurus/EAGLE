import logging
import os

from EAGLE.lib.general import ConfBase, get_redis_server, setup_logging

CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
CONF_DIR_NAME = 'configs'
CONF_DIR_PATH = os.path.join(CONSTANTS_PATH, CONF_DIR_NAME)
DEFAULT_CONFIG = os.path.join(CONF_DIR_PATH, "default_config.ini")
LOG_CONFIG_NAME = 'log_conf.yaml'
LOGGER_NAME = 'EAGLE_logger'

PROFILES_SCAN_OUT = "profiles_scan.out"

setup_logging(os.path.join(CONF_DIR_PATH, LOG_CONFIG_NAME))
EAGLE_logger = logging.getLogger(LOGGER_NAME)


class ConfConstants(ConfBase):

    def __init__(self, config_path=DEFAULT_CONFIG):
        # GENERAL
        self.num_threads = 4
        self.redis_host = 'localhost'
        self.redis_port = 6379
        self.redis_queue_db = 0
        # ALIGNMENT
        self.muscle_exec_path = "muscle"
        self.emboss_inst_dir = ""
        self.hmmer_inst_dir = ""
        self.blast_inst_dir = ""
        # PHYLO
        self.fastme_exec_path = "fastme"
        # EAGLE
        self.mode = "genome"

        super(ConfConstants, self).__init__(config_path=config_path)

    def get_redis_server(self, restart=True):
        EAGLE_logger.info(get_redis_server(self.redis_host, self.redis_port, restart=restart))


conf_constants = ConfConstants()
