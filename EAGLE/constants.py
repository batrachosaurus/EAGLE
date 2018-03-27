import os
import logging
import logging.config
import yaml
import ConfigParser
from ConfigParser import NoSectionError, NoOptionError


constants_path = os.path.dirname(os.path.realpath(__file__))
conf_dir_name = 'configs'
conf_dir_path = os.path.join(constants_path, conf_dir_name)
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


class ConfConstants:

    muscle_inst_dir = ""

    def __init__(self):
        pass

    def update_by_config(self, config_path):
        config = _config_parser(config_path=config_path)
        try:
            muscle_inst_dir = config.get("ALIGNMENT", "muscle_inst_dir")
            self.muscle_inst_dir = None
            self.muscle_inst_dir = muscle_inst_dir
        except (NoSectionError, NoOptionError):
            pass


conf_constants = ConfConstants()
