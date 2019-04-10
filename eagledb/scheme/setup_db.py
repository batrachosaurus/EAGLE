from functools import reduce
from operator import getitem
from abc import ABCMeta

from eagle.constants import eagle_logger


class ObjInfo:
    __metaclass__ = ABCMeta  # probably not compatible with Python 3

    # json scheme (the keys must match attribute names defined in __init__)
    attr_scheme = dict()

    def __init__(self):
        pass

    def get_json(self):
        res_json = dict()
        for attr in self.attr_scheme:
            keys = self.attr_scheme[attr]
            keys_l = len(keys)
            for i in range(keys_l):
                if i > 0:
                    if i == keys_l - 1:
                        sub_res_json = reduce(getitem, keys[:i], res_json)
                        sub_res_json[keys[i]] = self.__dict__[attr]
                    else:
                        sub_res_json = reduce(getitem, keys[:i], res_json)
                        sub_res_json[keys[i]] = dict()
                else:
                    if i == keys_l - 1:
                        res_json[keys[i]] = self.__dict__[attr]
                    else:
                        res_json[keys[i]] = dict()
            keys = None
        return res_json

    @classmethod
    def load_from_dict(cls, in_dict):
        obj_info = cls()
        for attr in cls.attr_scheme:
            try:
                obj_info.__dict__[attr] = reduce(getitem, cls.attr_scheme[attr], in_dict)
            except KeyError:
                eagle_logger.warning("the path '%s' is absent in the input dict - can not find a value for '%s'" %
                                     (str(cls.attr_scheme[attr]), attr))
                continue
        return obj_info


class GenomeInfo(ObjInfo):

    # json keys
    org_name_key = "org_name"
    taxonomy_key = "taxonomy"
    ncbi_download_prefix_key = "ncbi_download_prefix"
    fna_path_key = "fna_path"
    btc_seqs_path_key = "btc_seqs_path"
    source_db_key = "source_db"
    is_repr_key = "is_repr"

    # json scheme (the keys must match attribute names defined in __init__)
    attr_scheme = {
        "org_name": (org_name_key,),
        "taxonomy": (taxonomy_key,),
        "ncbi_download_prefix": (ncbi_download_prefix_key,),
        "fna_path": (fna_path_key,),
        "btc_seqs_path": (btc_seqs_path_key,),
        "source_db": (source_db_key,),
        "is_repr": (is_repr_key,),
    }

    # default values
    org_name_0 = None
    taxonomy_0 = None
    ncbi_download_prefix_0 = None
    fna_path_0 = None
    btc_seqs_path_0 = None
    source_db_0 = None
    is_repr_0 = False

    def __init__(self,
                 org_name=org_name_0,
                 taxonomy=taxonomy_0,
                 ncbi_download_prefix=ncbi_download_prefix_0,
                 fna_path=fna_path_0,
                 btc_seqs_path=btc_seqs_path_0,
                 source_db=source_db_0,
                 is_repr=is_repr_0):

        # attribute names must match keys form GenomeInfo.attr_scheme
        self.org_name = org_name
        self.taxonomy = list()
        if taxonomy is not None:
            self.taxonomy = taxonomy
        self.ncbi_download_prefix = ncbi_download_prefix
        self.fna_path = fna_path
        self.btc_seqs_path = btc_seqs_path
        self.source_db = source_db
        self.is_repr = is_repr

        super(GenomeInfo, self).__init__()


class SeqProfileInfo(ObjInfo):

    # json keys
    name_key = "name"
    path_key = "path"
    type_key = "type"
    weight_key = "weight"

    # json scheme (the keys must match attribute names defined in __init__)
    attr_scheme = {
        "name": (name_key,),
        "path": (path_key,),
        "type": (type_key,),
        "weight": (weight_key,),
    }

    # default values
    name_0 = None
    path_0 = None
    type_0 = 'nucl'
    weight_0 = 1.0

    def __init__(self):

        super(SeqProfileInfo, self).__init__()
