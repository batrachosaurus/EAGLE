import ConfigParser
import logging
import logging.config
import os
import yaml
import redis

from EAGLE.lib.general import ConfBase, get_redis_server

constants_path = os.path.dirname(os.path.realpath(__file__))
conf_dir_name = 'configs'
conf_dir_path = os.path.join(constants_path, conf_dir_name)
DEFAULT_CONFIG = os.path.join(conf_dir_path, "default_config.ini")
log_config_name = 'log_conf.yaml'
logger_name = 'EAGLE_logger'


def _setup_logging(
        default_path=os.path.join(conf_dir_path, log_config_name),
        default_level=logging.INFO,
        env_key='LOG_CFG'):
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):

        with open(path, 'rt') as f:
            string = f.read()
            config = yaml.load(string)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)


_setup_logging()
EAGLE_logger = logging.getLogger(logger_name)

def _config_parser(config_path):
    """ Function parses config file and puts the result into an object of ConfigParser class
      :param config_path: path to config file
      :return: a ConfigParser object
      """
    config = ConfigParser.ConfigParser()
    config.read(config_path)
    return config


class ConfConstants(ConfBase):

    def __init__(self, config_path=DEFAULT_CONFIG):
        # GENERAL
        self.num_threads = 4
        self.redis_host = 'localhost'
        self.redis_port = 6379
        self.redis_queue_db = 0
        # ALIGNMENT
        self.muscle_exec_path = ""
        self.emboss_inst_dir = ""
        self.hmmer_inst_dir = ""
        self.blast_inst_dir = ""
        # PHYLO
        self.fastme_exec_path = ""

        super(ConfConstants, self).__init__(config_path=config_path)

    def get_redis_server(self, restart=True):
        EAGLE_logger.info(get_redis_server(self.redis_host, self.redis_port, restart=restart))


conf_constants = ConfConstants()
