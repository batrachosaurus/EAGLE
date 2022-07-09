import os
import logging
import logging.config

import yaml

from eaglib.constants import CONSTANTS_PATH


LOG_CONFIG_PATH = os.path.join(CONSTANTS_PATH, "_utils", "logging", "eagle_log_conf.yml")


def setup_logging(default_path,
                  default_level=logging.INFO,
                  env_key='LOG_CFG'):
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):

        with open(path, 'rt') as f:
            string = f.read()
            config = yaml.load(string, Loader=yaml.FullLoader)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)


def send_log_message(message, mes_type="info", logger=None):
    if isinstance(logger, str):
        logger_name = logger
        logger = None
        logger = logging.getLogger(logger_name)

    if isinstance(logger, logging.Logger):
        if mes_type.lower() in ("info", "i"):
            logger.info(message)
        if mes_type.lower() in ("warning", "warn", "w"):
            logger.warning(message)
        if mes_type.lower() in ("error", "err", "e"):
            logger.error(message)
    else:  # TODO: detect logger
        if mes_type.lower() in ("info", "i"):
            eagle_logger.info(message)
        if mes_type.lower() in ("warning", "warn", "w"):
            eagle_logger.warning(message)
        if mes_type.lower() in ("error", "err", "e"):
            eagle_logger.error(message)


setup_logging(LOG_CONFIG_PATH)
eagle_logger = logging.getLogger(name="eagle_logger")
