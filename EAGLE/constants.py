import os
import logging
import logging.config
import yaml


constants_path = os.path.dirname(os.path.realpath(__file__))
conf_dir_name = 'configs'
log_config_name = 'log_conf.yaml'
logger_name = 'EAGLE_logger'


def setup_logging(
        default_path=os.path.join(constants_path, conf_dir_name, log_config_name),
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


setup_logging()
EAGLE_logger = logging.getLogger(logger_name)
