from functools import reduce
from operator import getitem
from abc import ABCMeta, abstractmethod
from collections import defaultdict

from jsondler import JsonEntry

from eaglib.alignment import SeqsProfileInfo
from eaglib.seqs import GenomeInfo


SeqProfileInfo = SeqsProfileInfo


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
    from_root_key = "from_root"

    # default values
    all_genomes_0 = None
    btax_json_0 = None
    repr_profiles_0 = None
    global_dist_matrix_0 = None
    all_org_full_names_0 = None
    from_root_0 = None

    def __init__(self,
                 all_genomes=all_genomes_0,
                 btax_json=btax_json_0,
                 repr_profiles=repr_profiles_0,
                 global_dist_matrix=global_dist_matrix_0,
                 all_org_full_names=all_org_full_names_0,
                 from_root=from_root_0):

        self.all_genomes = all_genomes
        self.btax_json = btax_json
        self.repr_profiles = repr_profiles
        self.global_dist_matrix = global_dist_matrix
        self.all_org_full_names = all_org_full_names
        self.from_root = from_root

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
            "from_root": (cls.from_root_key,),
        }
