from EAGLE import conf_constants
from EAGLEdb import conf_constants_db


def update_by_config(config_path):
    conf_constants.update_by_config(config_path)
    conf_constants_db.update()
