# This code can have only standard Python imports
import ConfigParser
import sys


class ConfConstantsBase:

    def __init__(self, config_path):
        self.config = None
        self.update_by_config(config_path=config_path)

    def update_by_config(self, config_path):
        self.config = _config_parser(config_path=config_path)
        for param in self.__dict__.keys():
            new_value = None
            if self._get_from_config(self.config, 'GENERAL', param, fallback=None):
                new_value = self.config.get('GENERAL', param)
            elif self._get_from_config(self.config, 'crAtlasSNPdb', param, fallback=None):
                new_value = self.config.get('crAtlasSNPdb', param)
            elif self._get_from_config(self.config, 'dbSNPreview', param, fallback=None):
                new_value = self.config.get('dbSNPreview', param)
            if new_value:
                if type(self.__dict__[param]) is list:
                    self.__dict__[param] = [elm.strip() for elm in new_value.split(",")]
                elif type(self.__dict__[param]) is bool:
                    if new_value.lower() in ("0", "false", "f", "no", "n", "na", "none", "null"):
                        self.__dict__[param] = False
                    else:
                        self.__dict__[param] = True
                elif new_value.lower() in ("0", "false", "f", "no", "n", "na", "none", "null"):
                    self.__dict__[param] = None
                else:
                    self.__dict__[param] = type(self.__dict__[param])(new_value)

    @staticmethod
    def _get_from_config(config_obj, section, option, fallback):
        try:
            return config_obj.get(section, option)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            return fallback

def _config_parser(config_path):
    """ Function parses config file and puts the result into an object of ConfigParser class
      :param config_path: path to config file
      :return: a ConfigParser object
      """
    config = ConfigParser.ConfigParser()
    config.read(config_path)
    return config


def get_config_parameter(config, section, parameter, fallback=None):
    pass


def worker(kwargs, use_try=False):
    func = kwargs.pop('function', None)
    if 'try_err_message' in kwargs.keys():
        use_try = True
    logger = kwargs.get('logger', None)
    if func:
        if use_try:
            try:
                func(**kwargs)
            except:
                if logger:
                    logger.warning("%s %s" % (kwargs['try_err_message'], sys.exc_info()))
                else:
                    print "%s %s" % (kwargs['try_err_message'], sys.exc_info())
        else:
            func(**kwargs)
    else:
        if logger:
            logger.warning("No function to run")
        else:
            print "No function to run"


def read_fasta_to_dict(fasta_path):
    fasta_dict = dict()
    seq_list = list()
    title = None
    fasta_f = open(fasta_path)
    for line_ in fasta_f:
        line = None
        line = line_.strip()
        if not line:
            continue
        if line[0] == ">":
            if title:
                fasta_dict[title] = "".join(seq_list)
                seq_list = list()
                title = None
            title = line[1:]
        else:
            seq_list.append(line)
    if title:
        fasta_dict[title] = "".join(seq_list)
        seq_list = list()
        title = None
    return fasta_dict