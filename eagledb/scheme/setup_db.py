from functools import reduce
from operator import getitem
from abc import ABCMeta, abstractmethod

from eagle.constants import eagle_logger


class JsonEntry(object):
    __metaclass__ = ABCMeta  # probably not compatible with Python 3

    def __init__(self):
        pass

    @staticmethod
    @abstractmethod
    def attr_scheme():
        """
        json scheme (the keys must match attribute names defined in __init__)
        CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        don't include it to the result dict
        :return: {attr_name: (path, in, json, entry,), ...}
        """
        return dict()

    def get_json(self):
        res_json = dict()
        for attr in self.attr_scheme():
            keys = self.attr_scheme()[attr]
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
        json_entry = cls()
        for attr in cls.attr_scheme():
            try:
                json_entry.__dict__[attr] = reduce(getitem, cls.attr_scheme()[attr], in_dict)
            except KeyError:
                eagle_logger.warning("the path '%s' is absent in the input dict - can not find a value for '%s'" %
                                     (str(cls.attr_scheme()[attr]), attr))
                continue
        return json_entry


class GenomeInfo(JsonEntry):

    # json keys
    org_name_key = "org_name"
    taxonomy_key = "taxonomy"
    ncbi_download_prefix_key = "ncbi_download_prefix"
    fna_path_key = "fna_path"
    btc_seqs_key = "btc_seqs"
    btc_seqs_fasta_key = "fasta"
    btc_seqs_id_key = "ids"
    source_db_key = "source_db"
    is_repr_key = "is_repr"

    # default values
    org_name_0 = None
    taxonomy_0 = None
    ncbi_download_prefix_0 = None
    fna_path_0 = None
    btc_seqs_fasta_0 = None
    btc_seqs_id_0 = None  # must be a dict: {btc_seq_id: btc_seq_profile_name}
    source_db_0 = None
    is_repr_0 = False

    def __init__(self,
                 org_name=org_name_0,
                 taxonomy=taxonomy_0,
                 ncbi_download_prefix=ncbi_download_prefix_0,
                 fna_path=fna_path_0,
                 btc_seqs_fasta=btc_seqs_fasta_0,
                 btc_seqs_id=btc_seqs_id_0,
                 source_db=source_db_0,
                 is_repr=is_repr_0):

        # attribute names must match keys form GenomeInfo.attr_scheme()
        self.org_name = org_name
        self.taxonomy = list()
        if taxonomy is not None:
            self.taxonomy = taxonomy
        self.ncbi_download_prefix = ncbi_download_prefix
        self.fna_path = fna_path
        self.btc_seqs_fasta = btc_seqs_fasta
        self.btc_seqs_id = btc_seqs_id
        if self.btc_seqs_id is None:
            self.btc_seqs_id = dict()
        self.source_db = source_db
        self.is_repr = is_repr

        super(GenomeInfo, self).__init__()

    @staticmethod
    def attr_scheme():
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "org_name": (GenomeInfo.org_name_key,),
            "taxonomy": (GenomeInfo.taxonomy_key,),
            "ncbi_download_prefix": (GenomeInfo.ncbi_download_prefix_key,),
            "fna_path": (GenomeInfo.fna_path_key,),
            "btc_seqs_fasta": (GenomeInfo.btc_seqs_key, GenomeInfo.btc_seqs_fasta_key,),
            "btc_seqs_id": (GenomeInfo.btc_seqs_key, GenomeInfo.btc_seqs_id_key,),
            "source_db": (GenomeInfo.source_db_key,),
            "is_repr": (GenomeInfo.is_repr_key,),
        }


class SeqProfileInfo(JsonEntry):

    # json keys
    name_key = "name"
    path_key = "path"
    seq_type_key = "type"
    weight_key = "weight"

    # default values
    name_0 = None
    path_0 = None
    seq_type_0 = 'nucl'
    weight_0 = 1.0

    def __init__(self,
                 name=name_0,
                 path=path_0,
                 seq_type=seq_type_0,
                 weight=weight_0):

        # attribute names must match keys form SeqProfileInfo.attr_scheme()
        self.name = name
        self.path = path
        self.seq_type = seq_type
        self.weight = weight

        super(SeqProfileInfo, self).__init__()

    @staticmethod
    def attr_scheme():
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "name": (SeqProfileInfo.name_key,),
            "path": (SeqProfileInfo.path_key,),
            "seq_type": (SeqProfileInfo.seq_type_key,),
            "weight": (SeqProfileInfo.weight_key,),
        }
