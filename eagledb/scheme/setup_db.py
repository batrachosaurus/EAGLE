from functools import reduce
from operator import getitem
from abc import ABCMeta, abstractmethod
from collections import defaultdict

from eagle.constants import eagle_logger


class JsonEntry(object):
    __metaclass__ = ABCMeta  # probably not compatible with Python 3

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
                        if keys[i] not in sub_res_json:
                            sub_res_json[keys[i]] = dict()
                else:
                    if i == keys_l - 1:
                        res_json[keys[i]] = self.__dict__[attr]
                    elif keys[i] not in res_json:
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
    genome_id_key = "genome_id"
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
    genome_id_0 = None
    org_name_0 = None
    taxonomy_0 = None
    ncbi_download_prefix_0 = None
    fna_path_0 = None
    btc_seqs_fasta_0 = None
    btc_seqs_id_0 = None  # must be a dict: {btc_seq_id: btc_seq_profile_name}
    source_db_0 = None
    is_repr_0 = False

    def __init__(self,
                 genome_id=genome_id_0,
                 org_name=org_name_0,
                 taxonomy=taxonomy_0,
                 ncbi_download_prefix=ncbi_download_prefix_0,
                 fna_path=fna_path_0,
                 btc_seqs_fasta=btc_seqs_fasta_0,
                 btc_seqs_id=btc_seqs_id_0,
                 source_db=source_db_0,
                 is_repr=is_repr_0):

        # attribute names must match keys form GenomeInfo.attr_scheme()
        self.genome_id = genome_id
        self.org_name = org_name
        self.taxonomy = list()
        if taxonomy is not None:
            self.taxonomy = taxonomy
        self.ncbi_download_prefix = ncbi_download_prefix
        self.fna_path = fna_path
        self.btc_seqs_fasta = btc_seqs_fasta
        self.btc_seqs_id = btc_seqs_id
        if self.btc_seqs_id is None:
            self.btc_seqs_id = defaultdict(str)  # {btc_seq_id: btc_profile_name}
        self.source_db = source_db
        self.is_repr = is_repr

    @staticmethod
    def attr_scheme():
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "genome_id": (GenomeInfo.genome_id_key,),
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


class BtaxInfo(JsonEntry):

    # json keys
    name_key = "name"
    genomes_key = "genomes"
    btax_fna_key = "btax_fna"
    fna_id_key = "fna_id"
    blastdb_key = "blastdb"
    repr_profiles_key = "repr_profiles"
    ref_tree_key = "ref_tree"
    ref_tree_newick_key = "newick"
    ref_tree_full_names_key = "full_names"
    distance_key = "distance"
    mean_d_key = "mean"
    median_d_key = "median"

    # default values
    name_0 = None
    genomes_0 = None
    btax_fna_0 = None
    fna_id_0 = None
    blastdb_0 = None
    repr_profiles_0 = None
    ref_tree_newick_0 = None
    ref_tree_full_names_0 = None
    mean_d_0 = float()
    median_d_0 = float()

    def __init__(self,
                 name=name_0,
                 genomes=genomes_0,
                 btax_fna=btax_fna_0,
                 fna_id=fna_id_0,
                 blastdb=blastdb_0,
                 repr_profiles=repr_profiles_0,
                 ref_tree_newick=ref_tree_newick_0,
                 ref_tree_full_names=ref_tree_full_names_0,
                 mean_d=mean_d_0,
                 median_d=median_d_0):

        # attribute names must match keys form SeqProfileInfo.attr_scheme()
        self.name = name
        self.genomes=genomes
        if self.genomes is None:
            self.genomes = list()
        self.btax_fna = btax_fna
        self.fna_id = fna_id
        if self.fna_id is None:
            self.fna_id = defaultdict(str)  # {seq_id: fna_path}
        self.blastdb = blastdb
        self.repr_profiles = repr_profiles
        if self.repr_profiles is None:
            self.repr_profiles = list()
        self.ref_tree_newick = ref_tree_newick
        self.ref_tree_full_names = ref_tree_full_names
        if self.ref_tree_full_names is None:
            self.ref_tree_full_names = defaultdict(str)  # {short_name: full_name}###
        self.mean_d = mean_d
        self.median_d = median_d

    @staticmethod
    def attr_scheme():
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "name": (BtaxInfo.name_key,),
            "genomes": (BtaxInfo.genomes_key,),
            "btax_fna": (BtaxInfo.btax_fna_key,),
            "fna_id": (BtaxInfo.fna_id_key,),
            "blastdb": (BtaxInfo.blastdb_key,),
            "repr_profiles": (BtaxInfo.repr_profiles_key,),
            "ref_tree_newick": (BtaxInfo.ref_tree_key, BtaxInfo.ref_tree_newick_key,),
            "ref_tree_full_names": (BtaxInfo.ref_tree_key, BtaxInfo.ref_tree_full_names_key,),
            "mean_d": (BtaxInfo.distance_key, BtaxInfo.mean_d_key,),
            "median_d": (BtaxInfo.distance_key, BtaxInfo.median_d_key,),
        }


class DBInfo(JsonEntry):

    #json keys
    all_genomes_key = "all_genomes"
    btax_json_key = "btax_json"
    repr_profiles_key = "repr_profiles"
    global_dist_matrix_key = "global_dist_matrix"
    all_org_full_names_key = "all_org_full_names"

    # default values
    all_genomes_0 = None
    btax_json_0 = None
    repr_profiles_0 = None
    global_dist_matrix_0 = None
    all_org_full_names_0 = None

    def __init__(self,
                 all_genomes=all_genomes_0,
                 btax_json=btax_json_0,
                 repr_profiles=repr_profiles_0,
                 global_dist_matrix=global_dist_matrix_0,
                 all_org_full_names=all_org_full_names_0):

        self.all_genomes = all_genomes
        self.btax_json = btax_json
        self.repr_profiles = repr_profiles
        self.global_dist_matrix = global_dist_matrix
        self.all_org_full_names = all_org_full_names

    @staticmethod
    def attr_scheme():
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "all_genomes": (DBInfo.all_genomes_key,),
            "btax_json": (DBInfo.btax_json_key,),
            "repr_profiles": (DBInfo.repr_profiles_key,),
            "global_dist_matrix": (DBInfo.global_dist_matrix_key,),
            "all_org_full_names": (DBInfo.all_org_full_names_key,),
        }
