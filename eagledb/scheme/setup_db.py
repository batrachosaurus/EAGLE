from functools import reduce
from operator import getitem
from abc import ABCMeta, abstractmethod
from collections import defaultdict


class JsonEntry(object):
    __metaclass__ = ABCMeta  # not compatible with Python 3

    @classmethod
    @abstractmethod
    def attr_scheme(cls):
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
                print("WARNING: the path '%s' is absent in the input dict - can not find a value for '%s'" %
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
    taxonomy_0 = list()
    ncbi_download_prefix_0 = None
    fna_path_0 = None
    btc_seqs_fasta_0 = None
    btc_seqs_id_0 = defaultdict(str)  # {btc_seq_id: btc_seq_profile_name}
    source_db_0 = None
    is_repr_0 = False

    def __init__(self,
                 genome_id=genome_id_0,
                 org_name=org_name_0,
                 taxonomy=None,
                 ncbi_download_prefix=ncbi_download_prefix_0,
                 fna_path=fna_path_0,
                 btc_seqs_fasta=btc_seqs_fasta_0,
                 btc_seqs_id=None,
                 source_db=source_db_0,
                 is_repr=is_repr_0):

        # attribute names must match keys form GenomeInfo.attr_scheme()
        if taxonomy is None:
            taxonomy = self.taxonomy_0
        if btc_seqs_id is None:
            btc_seqs_id = self.btc_seqs_id_0
        
        self.genome_id = genome_id
        self.org_name = org_name
        self.taxonomy = taxonomy
        self.ncbi_download_prefix = ncbi_download_prefix
        self.fna_path = fna_path
        self.btc_seqs_fasta = btc_seqs_fasta
        self.btc_seqs_id = btc_seqs_id
        self.source_db = source_db
        self.is_repr = is_repr

    @classmethod
    def attr_scheme(cls):
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "genome_id": (cls.genome_id_key,),
            "org_name": (cls.org_name_key,),
            "taxonomy": (cls.taxonomy_key,),
            "ncbi_download_prefix": (cls.ncbi_download_prefix_key,),
            "fna_path": (cls.fna_path_key,),
            "btc_seqs_fasta": (cls.btc_seqs_key, cls.btc_seqs_fasta_key,),
            "btc_seqs_id": (cls.btc_seqs_key, cls.btc_seqs_id_key,),
            "source_db": (cls.source_db_key,),
            "is_repr": (cls.is_repr_key,),
        }

    @classmethod
    def org_name_from_dict(cls, in_dict):
        return in_dict["org_name"]


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

    @classmethod
    def attr_scheme(cls):
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "name": (cls.name_key,),
            "path": (cls.path_key,),
            "seq_type": (cls.seq_type_key,),
            "weight": (cls.weight_key,),
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
    genomes_0 = list()
    btax_fna_0 = None
    fna_id_0 = defaultdict(str)  # {seq_id: fna_path}
    blastdb_0 = None
    repr_profiles_0 = dict()
    ref_tree_newick_0 = None
    ref_tree_full_names_0 = defaultdict(str)  # {short_name: full_name}
    mean_d_0 = float()
    median_d_0 = float()

    def __init__(self,
                 name=name_0,
                 genomes=None,
                 btax_fna=btax_fna_0,
                 fna_id=None,
                 blastdb=blastdb_0,
                 repr_profiles=None,
                 ref_tree_newick=ref_tree_newick_0,
                 ref_tree_full_names=None,
                 mean_d=mean_d_0,
                 median_d=median_d_0):

        # attribute names must match keys form SeqProfileInfo.attr_scheme()
        if genomes is None:
            genomes = self.genomes_0
        if fna_id is None:
            fna_id = self.fna_id_0
        if repr_profiles is None:
            repr_profiles = self.repr_profiles_0
        if ref_tree_full_names is None:
            ref_tree_full_names = self.ref_tree_full_names_0
        
        self.name = name
        self.genomes=genomes
        self.btax_fna = btax_fna
        self.fna_id = fna_id
        self.blastdb = blastdb
        self.repr_profiles = repr_profiles
        self.ref_tree_newick = ref_tree_newick
        self.ref_tree_full_names = ref_tree_full_names
        self.mean_d = mean_d
        self.median_d = median_d

    @classmethod
    def attr_scheme(cls):
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "name": (cls.name_key,),
            "genomes": (cls.genomes_key,),
            "btax_fna": (cls.btax_fna_key,),
            "fna_id": (cls.fna_id_key,),
            "blastdb": (cls.blastdb_key,),
            "repr_profiles": (cls.repr_profiles_key,),
            "ref_tree_newick": (cls.ref_tree_key, cls.ref_tree_newick_key,),
            "ref_tree_full_names": (cls.ref_tree_key, cls.ref_tree_full_names_key,),
            "mean_d": (cls.distance_key, cls.mean_d_key,),
            "median_d": (cls.distance_key, cls.median_d_key,),
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

    @classmethod
    def attr_scheme(cls):
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "all_genomes": (cls.all_genomes_key,),
            "btax_json": (cls.btax_json_key,),
            "repr_profiles": (cls.repr_profiles_key,),
            "global_dist_matrix": (cls.global_dist_matrix_key,),
            "all_org_full_names": (cls.all_org_full_names_key,),
        }
