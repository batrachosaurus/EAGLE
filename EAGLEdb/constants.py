import os

from EAGLE.constants import DEFAULT_CONFIG
from EAGLE.lib.general import ConfConstantsBase

CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))

BACTERIA_LIST_F_NAME = "bacteria.json"
ANALYZED_BACTERIA_F_NAME = "analyzed_bacteria.p"
BACT_FAM_F_NAME = "bact_fam.json"


class ConfConstants(ConfConstantsBase):

    def __init__(self, config_path=DEFAULT_CONFIG):
        # Bacteria db
        self.only_repr=False

        super(ConfConstants, self).__init__(config_path=config_path)


conf_constants_db = ConfConstants()
