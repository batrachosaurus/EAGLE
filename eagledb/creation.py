# each path for DB file should be stored as relative path from db_dir
# each time a DB file is used it should be accessed with os.path.join(db_dir, relative_f_path)

import os
import json
import argparse
from copy import deepcopy
import multiprocessing as mp
from collections import Iterable, defaultdict

import numpy as np
import pandas as pd
import networkx as nx

from eaglib._utils.logging import eagle_logger
from eaglib._utils.workers import process_worker
from eaglib.seqs import SeqsDict, load_fasta_to_dict, reduce_seq_names
from eaglib.alignment import SeqsProfileInfo, SeqsProfile
from eaglib.phylo import DistanceMatrix
from eagledb.constants import conf_constants, conf_constants_lib, GLOBAL_DIST_MATRIX, SHORT_TO_FULL_ORG_NAMES
from eagledb.scheme import BtaxInfo, GenomeInfo


def create_cmd():
    parser = argparse.ArgumentParser()
    pass
    # TODO: genomes_table_path -> genomes_table


def create(db_dir: str,
           genomes_table: pd.DataFrame,
           btax_class_profiles: Iterable,
           btax_repr_profiles: Iterable,
           btax_level=None,
           num_threads=None,
           config_path=None,
           **kwargs):
    """

    :param db_dir:
    :param genomes_table:
        id
        name
        taxonomy - fixed positions list of taxonomic units
        btc_seqs - list of paths to FASTA with sequences used for basic taxons classification
                   sequence names are profile names, NO paralogs
                   sequences can be joined into single file or distributed between several files
        fna_seq - list of paths or links to the genome files (FASTA or archived FASTA)
                  sequences can be joined into single file or distributed between several files
    :param btax_class_profiles:  required for btc profiles meta data
    :param btax_repr_profiles:
    :param btax_level:
    :param num_threads:
    :param config_path:
    :param kwargs:
    :return:
    """
    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_lib.update_by_config(config_path=config_path)
    if num_threads is not None:
        conf_constants.num_threads = num_threads
    if btax_level is not None:
        conf_constants.btax_level = btax_level

    btax_dict, short2full_seq_names = get_btax_dict(db_dir=db_dir,
                                                    genomes_list=genomes_table.to_records(),
                                                    btax_level=conf_constants.btax_level,
                                                    btc_profiles=btax_class_profiles,
                                                    num_threads=conf_constants.num_threads)

    btax_names = tuple(btax_dict.keys())
    to_complete_tasks = list()
    for btax_name in btax_names:
        to_complete_tasks = {"function": complete_btax,
                             "btax_info": btax_dict[btax_name],
                             "btr_profiles": btax_repr_profiles,
                             "db_dir": db_dir}###
    pool = mp.Pool(conf_constants.num_threads)
    complete_result = pool.map(process_worker, to_complete_tasks)
    pool.close()
    pool.join()
    for i, btax_info in enumerate(complete_result):
        btax_dict[btax_names[i]] = btax_info

    btax_dict = get_btax_blastdb(btax_dict,
                                 db_dir=db_dir,
                                 btr_profiles=btax_repr_profiles,
                                 num_threads=conf_constants.num_threads)

    repr_profiles_path = create_profiles_db(btax_dict,
                                            db_dir=db_dir,
                                            profiles_db_name=PROFILES_DB_NAME,
                                            method="hmmer",
                                            hmmer_inst_dir=conf_constants_lib.hmmer_inst_dir,
                                            config_path=config_path,
                                            logger=eagle_logger)
    with open(os.path.join(db_dir, BTAX_JSON_NAME), "w") as btax_json_f:
        json.dump(btax_dict, btax_json_f, indent=2)  # maybe btax_dict will be dumped in get_btax_dict
    db_info = DBInfo(
        all_genomes=os.path.join(db_dir, BACTERIA_LIST_F_NAME),
        btax_json=os.path.join(db_dir, BTAX_JSON_NAME),
        repr_profiles=repr_profiles_path,
        global_dist_matrix=os.path.join(db_dir, BACTERIA_GLOBAL_DIST_MATRIX),
        all_org_full_names=os.path.join(db_dir, BACTERIA_SHORT_TO_FULL_ORG_NAMES),
        from_root=from_root,
    ).get_json()
    with open(os.path.join(db_dir, DB_INFO_NAME), "w") as db_info_f:
        json.dump(db_info, db_info_f, indent=2)
    return db_info


def get_btax_dict(db_dir,
                  genomes_list,
                  btax_level,
                  btc_profiles,
                  **kwargs):

    btc_info_dict = defaultdict(SeqsProfileInfo)
    for btc_profile_dict in btc_profiles:
        btc_profile_info = SeqsProfileInfo.load_from_dict(btc_profile_dict)
        btc_info_dict[btc_profile_info.name] = btc_profile_info

    btax_dict, btc_fasta_dict, genome_keys = genomes2btax(genomes_list=genomes_list, btax_level=btax_level)
    short_to_full_seq_names = reduce_seq_names({gk: None for gk in genome_keys}, num_letters=10, num_words=3)[1]
    btc_aln_dict, btc_dist_dict = get_btc_alignments(btc_fasta_dict=btc_fasta_dict, btc_info_dict=btc_info_dict,
                                                     short_to_full_seq_names=short_to_full_seq_names, db_dir=db_dir,
                                                     **kwargs)

    global_dist_matr = get_global_dist(btc_dist_dict, btc_profiles, genome_keys)
    # global_dist_matr is needed only for standardize_btax
    global_dist_matr_path = os.path.join(db_dir, GLOBAL_DIST_MATRIX)
    short_to_full_seq_names_path = os.path.join(db_dir, SHORT_TO_FULL_ORG_NAMES)
    short_to_full_seq_names = global_dist_matr.dump(matrix_path=global_dist_matr_path, matr_format="phylip")
    with open(short_to_full_seq_names_path, "w") as short_to_full_org_names_f:
        json.dump(short_to_full_seq_names, short_to_full_org_names_f, indent=2)

    eagle_logger.info("base taxons standardisation started")
    btax_dict = standardize_btax(btax_dict=btax_dict, global_dist_matr=global_dist_matr)
    eagle_logger.info("base taxons standardisation finished")

    for btax_name in btax_dict:
        btax_keys = [GenomeInfo.key_from_dict(g_dict) for g_dict in btax_dict[btax_name].genomes]
        for btc_profile_name in btc_info_dict:
            btax_dict[btax_name].repr_alns[btc_profile_name] = btc_aln_dict[btc_profile_name][btax_keys]
            btax_dict[btax_name].repr_dists[btc_profile_name] = btc_dist_dict[btc_profile_name][btax_keys]

    return btax_dict, short_to_full_seq_names

"""
    full_to_short_seq_names = {v: k for k, v in short_to_full_seq_names.items()}
    for btax_name in btax_dict:
        btax_keys = [GenomeInfo.key_from_dict(g_dict) for g_dict in btax_dict[btax_name].genomes]
        ###
        if btr_profiles is not None:
            pass
        else:
            btax_dict[btax_name].mean_d = global_dist_matr[btax_orgs].mean_dist###
            btax_dict[btax_name].median_d = global_dist_matr[btax_orgs].median_dist###
            if len(btax_orgs) > 2:
                btax_dict[btax_name].ref_tree_newick = build_tree_by_dist(global_dist_matr[btax_orgs],###
                                                                          tree_name=btax_name+"_tree",
                                                                          options={"-T": num_threads}).newick
                btax_btc_aln_dict = dict()
                for btc_profile_name, btc_aln in btc_aln_dict.items():
                    btax_btc_aln = btc_aln[btax_orgs].improve_aln(inplace=False)
                    btax_btc_aln.aln_name = btax_name + "_" + btc_profile_name
                    btax_btc_aln_dict[btc_profile_name] = deepcopy(btax_btc_aln)
                btax_dict[btax_name].repr_profiles = generate_btax_profiles(btax_btc_aln_dict,
                                                                            db_dir=db_dir,
                                                                            btax_name=btax_name,
                                                                            method="hmmer")
        ###
        btax_dict[btax_name].ref_tree_full_names = \
            {full_to_short_seq_names[btax_key]: btax_key for btax_key in btax_keys}
        btax_dict[btax_name] = btax_dict[btax_name].get_json()
    return btax_dict
"""

def genomes2btax(genomes_list, btax_level):
    btax_dict = defaultdict(BtaxInfo)
    btc_fasta_dict = defaultdict(dict)
    genome_keys = set()
    for genome_dict in genomes_list:
        if not genome_dict:
            continue
        genome_info = GenomeInfo(
            genome_id=genome_dict["id"],
            org_name=genome_dict["name"],
            taxonomy=genome_dict["taxonomy"],
            fna_seq_fasta=genome_dict["fna_seq"],
            repr_seq_fasta=genome_dict["btc_seqs"]
        )
        genome_key = genome_info.key  # to guarantee unique seq names
        genome_keys.add(genome_key)
        try:
            btax_name = genome_info.taxonomy[-btax_level]
        except IndexError:
            btax_name = genome_info.taxonomy[0]
        btax_dict[btax_name].genomes.append(genome_info.get_json())
        if btax_dict[btax_name].name is None:
            btax_dict[btax_name].name = btax_name
        repr_fasta_dict = dict()
        for fasta_path in genome_info.repr_seq_fasta:
            repr_fasta_dict.update(load_fasta_to_dict(fasta_path))
        btc_seqs_dict = SeqsDict.load_from_dict(repr_fasta_dict)
        for btc_seq_id in btc_seqs_dict:
            btc_fasta_dict[btc_seq_id][genome_key] = btc_seqs_dict[btc_seq_id]
        del genome_key
    return btax_dict, btc_fasta_dict, genome_keys


def get_btc_alignments(btc_fasta_dict, btc_info_dict, short_to_full_seq_names, db_dir, **kwargs):
    btc_dist_dict = dict()
    btc_aln_dict = dict()
    for btc_profile_name in btc_fasta_dict:
        btc_mult_aln = SeqsDict.load_from_dict(btc_fasta_dict[btc_profile_name],
                                               seqs_type=btc_info_dict[btc_profile_name].seq_type,
                                               logger=eagle_logger,
                                               **kwargs).construct_mult_aln(
            aln_name=btc_profile_name+"_aln",
            tmp_dir=kwargs.get("aln_tmp_dir", "mult_aln_tmp"),
            method=conf_constants.btc_profile_aln_method,
            num_threads=conf_constants.num_threads,
            op=15.0,
            ep=0.1,
            **kwargs
        )
        btc_mult_aln.short_to_full_seq_names = short_to_full_seq_names.copy()
        if btc_mult_aln.aln_type.lower() in SeqsDict.nucl_types:
            btc_mult_aln.dist_matr_options.update({"--dna": "p", "-T": conf_constants.num_threads, "-f": 6})
        btc_mult_aln.improve_aln(inplace=True)
        btc_dist_dict[btc_profile_name] = btc_mult_aln.get_distance_matrix()
        if kwargs.get("save_alignments", False):
            btc_mult_aln.dump(fname=os.path.join(db_dir, btc_mult_aln.aln_name + ".fasta"), format="fasta")
        btc_aln_dict[btc_profile_name] = deepcopy(btc_mult_aln)
    return btc_aln_dict, btc_dist_dict


def get_global_dist(btc_dist_dict, btc_profiles, genome_keys: set):
    seqs_order = {gk: i for i, gk in enumerate(genome_keys)}
    nseqs = len(seqs_order)
    global_dist_matrix = DistanceMatrix(seqs_order=seqs_order,
                                        matr=np.zeros((nseqs, nseqs)),
                                        aln_type="btc_global")
    global_matrix_orgs = global_dist_matrix.seq_names
    sumw = 0.0
    for btc_profile in btc_profiles:
        btc_profile_info = SeqsProfileInfo.load_from_dict(btc_profile)
        btc_profile_matr = btc_dist_dict[btc_profile_info.name]
        btc_profile_orgs = btc_profile_matr.seq_names
        dist_0 = btc_profile_matr.mean_dist
        eagle_logger.info("%s distance mean %s" % (btc_profile_info.name, dist_0))
        for i, seq_name in enumerate(global_dist_matrix.seq_names):
            if seq_name in btc_profile_orgs:
                btc_seq_dist = btc_profile_matr[seq_name]
                seq_dists = pd.Series(([0.0] * i + [btc_seq_dist.get(seq_name_, dist_0)
                                       for seq_name_ in global_matrix_orgs[i:]]),
                                      index=global_matrix_orgs)
            else:
                seq_dists = pd.Series([0.0] * i + [dist_0] * (nseqs-i), index=global_dist_matrix.seq_names)
                seq_dists[seq_name] = 0.0
            global_dist_matrix[seq_name] += seq_dists * float(btc_profile_info.weight)
        sumw += float(btc_profile_info.weight)
    global_dist_matrix.matr = global_dist_matrix.matr / sumw
    return global_dist_matrix


def standardize_btax(btax_dict, global_dist_matr, k_max=None, k_min=None):
    if k_max is None:
        k_max = conf_constants.k_max
    if k_min is None:
        k_min = conf_constants.k_min

    assert isinstance(btax_dict, defaultdict) and btax_dict.default_factory is BtaxInfo, \
        "ERROR: the value for btax_dict should be defaultdict(BtaxInfo)"

    btax_to_merge = set()
    do_round = True
    while do_round:
        if btax_to_merge:
            btax_merge_edges = set()
            btax_dist_matr = DistanceMatrix(seqs_order={btax_name_: i for i, btax_name_ in enumerate(btax_dict.keys())},
                                            matr=np.zeros(shape=(len(btax_dict), len(btax_dict))))
            for btax_name in btax_to_merge:
                for btax_name_ in btax_dict:
                    if btax_name != btax_name_ and btax_dist_matr[btax_name][btax_name_] == 0.0:
                        btax_dists = btax_dist_matr[btax_name]
                        btax_dists[btax_name_] = get_btax_dist(
                            btax1_orgs=[GenomeInfo.key_from_dict(g_dict) for g_dict in btax_dict[btax_name].genomes],
                            btax2_orgs=[GenomeInfo.key_from_dict(g_dict) for g_dict in btax_dict[btax_name_].genomes],
                            global_dist_matr=global_dist_matr
                        )
                        btax_dist_matr[btax_name] = btax_dists

                btax_dists = btax_dist_matr[btax_name]
                btax_merge_edges.add(frozenset([btax_name, btax_dists[btax_dists > 0.0].idxmin()]))

            g = nx.Graph()
            g.add_nodes_from(btax_dist_matr.seq_names)
            g.add_edges_from(btax_merge_edges)

            for comp in nx.connected_components(g):
                if len(comp) >= 2:
                    root_btax_name = None
                    btax_names = list()
                    btax_genomes = list()
                    for btax_name in comp:
                        if root_btax_name is None and btax_name not in btax_to_merge:
                            root_btax_name = btax_name
                        btax_names.append(btax_dict[btax_name].name)
                        btax_genomes.extend(btax_dict.pop(btax_dict[btax_name].genomes))
                    if root_btax_name is None:
                        root_btax_name = sorted(g.degree(comp), key=lambda nd: nd[1], reverse=True)[0][0]
                    btax_dict[root_btax_name.replace("_related", "") + "_related"] = \
                        BtaxInfo(name=",".join(btax_names), genomes=btax_genomes)
                    del root_btax_name
            btax_to_merge = set()

        btax_names = tuple(btax_dict.keys())
        filt_tasks = list()
        for btax_name in btax_names:
            filt_tasks.append({"function": filter_btax,
                               "btax_info": btax_dict[btax_name],
                               "global_dist_matr": global_dist_matr,
                               "k_max": k_max})
        btax_filt_pool = mp.Pool(conf_constants.num_threads)
        filt_result = btax_filt_pool.map(process_worker, filt_tasks)
        btax_filt_pool.close()
        btax_filt_pool.join()
        for i, btax_info in enumerate(filt_result):
            btax_dict[btax_names[i]] = btax_info
            if len(btax_info.genomes) < k_min:
                btax_to_merge.add(btax_names[i])
        if btax_to_merge:
            do_round = True
    return btax_dict


def filter_btax(btax_info: BtaxInfo, global_dist_matr, k_max=None):
    if k_max is None:
        k_max = conf_constants.k_max

    btax_genome_keys = list()
    for i, genome_info_dict in enumerate(btax_info.genomes):
        btax_genome_keys.append(GenomeInfo.load_from_dict(genome_info_dict).key)
    btax_dist_matr = global_dist_matr[btax_genome_keys]
    genomes_to_remove = list()
    while len(btax_genome_keys) > k_max:
        min_dist_sum = None
        min_closest_dist = None
        i_to_remove = None
        for i, genome_key in enumerate(btax_genome_keys):
            genome_dist = btax_dist_matr[genome_key][btax_genome_keys]
            dist_sum = genome_dist.sum()
            closest_dist = genome_dist[~genome_dist.index.isin([genome_key])].min()
            if min_closest_dist is None or closest_dist < min_closest_dist or \
                    (closest_dist == min_closest_dist and dist_sum < min_dist_sum):
                min_closest_dist = closest_dist
                min_dist_sum = dist_sum
                i_to_remove = i
        if i_to_remove is not None:
            genomes_to_remove.append(i_to_remove)
            del btax_genome_keys[i_to_remove]
    for min_dist_i in genomes_to_remove:
        del btax_info.genomes[min_dist_i]
    return btax_info


def get_btax_dist(btax1_orgs, btax2_orgs, global_dist_matr):
    n = 0
    sum_dist = 0.0
    for org in btax1_orgs:
        for org_ in btax2_orgs:
            sum_dist += global_dist_matr[org][org_]
            n += 1
    return sum_dist / float(n)


def complete_btax(btax_info: BtaxInfo, btr_profiles: Iterable, db_dir):###
    btax_dir = os.path.join(db_dir, btax_info.name)
    os.makedirs(btax_dir)
    get_btax_fna(btax_info=btax_info, btax_dir=btax_dir)###
    add_repr_profiles(btax_info=btax_info, btr_profiles=btr_profiles, btax_dir=btax_dir)###
    # TODO: find, align and get distances for repr_profiles seqs
    # TODO: build btax repr_profiles
    # TODO: build btax reference tree
    # TODO: build blastdb
    return btax_info


def get_btax_fna(btax_info: BtaxInfo, btax_dir: str) -> None:###
    # Inplace modification of BtaxInfo object
    btax_fna_path = os.path.join(btax_dir, btax_info.name + ".fasta")
    btax_fna_dict = dict()
    for btax_genome in btax_info.genomes:
        genome_info = GenomeInfo.load_from_dict(btax_genome)
        for seq_fasta in genome_info.fna_seq_fasta:
            if seq_fasta.startswith("http://") or seq_fasta.startswith("https://") or seq_fasta.startswith("ftp://"):
                local_seq_fasta = ...###
                # TODO: download fna
            else:
                local_seq_fasta = seq_fasta
            if local_seq_fasta:# is archive?
                extracted_seq_fasta = ...###
                # TODO: extract from archive
            else:
                extracted_seq_fasta = local_seq_fasta
            for seq_id, seq in load_fasta_to_dict(extracted_seq_fasta).items():
                if seq_id in btax_fna_dict:
                    btax_seq_id = ...###
                else:
                    btax_seq_id = seq_id
                btax_fna_dict[btax_seq_id] = seq
                btax_info.fna_id2genome[btax_seq_id] = (genome_info.key, seq_id)
                del btax_seq_id
            if seq_fasta != local_seq_fasta:
                os.remove(local_seq_fasta)
            if seq_fasta != extracted_seq_fasta:
                os.remove(extracted_seq_fasta)
            del local_seq_fasta
            del extracted_seq_fasta

    SeqsDict.load_from_dict(btax_fna_dict).dump(btax_fna_path, format='fasta')
    btax_info.btax_fna = btax_fna_path


def add_repr_profiles(btax_info: BtaxInfo, btr_profiles: Iterable, btax_dir):
    for repr_profile_dict in btr_profiles:
        repr_profile = SeqsProfile(SeqsProfileInfo.load_from_dict(repr_profile_dict),
                                   tmp_dir=os.path.join(btax_dir, "tmp"))
    pass###


def get_btax_blastdb(btax_dict, db_dir, btr_profiles=None, num_threads=None, config_path=None):
    if config_path:
        conf_constants_lib.update_by_config(config_path=config_path)
        conf_constants.update_by_config(config_path=config_path)
    if not num_threads:
        num_threads = conf_constants.num_threads
    else:
        conf_constants.num_threads = num_threads

    # Can be parallel
    for btax_name in btax_dict:
        btax_info = BtaxInfo.load_from_dict(btax_dict[btax_name])
        btax_info.btax_fna, btax_info.fna_id, downloaded_fna = get_btax_fna(btax_genomes=btax_info.genomes,
                                                                            btax_name=btax_name,
                                                                            btax_dir=db_dir)
        for i, btax_genome in enumerate(btax_info.genomes):
            genome_info = GenomeInfo.load_from_dict(btax_genome)
            if genome_info.genome_id in downloaded_fna:
                genome_info.fna_path = downloaded_fna[genome_info.genome_id]
                btax_info.genomes[i] = genome_info.get_json()
        btax_info.blastdb = create_btax_blastdb(btax_fna_path=btax_info.btax_fna,
                                                btax_name=btax_name,
                                                db_dir=db_dir,
                                                blast_inst_dir=conf_constants_lib.blast_inst_dir,
                                                logger=eagle_logger)
        if btr_profiles is not None:
            # create repr profile
            pass
        btax_dict[btax_name] = btax_info.get_json()
    return btax_dict


def create_profiles_db(btax_dict,
                       db_dir,
                       profiles_db_name=PROFILES_DB_NAME,
                       method="hmmer",
                       hmmer_inst_dir="",
                       config_path=None,  # TODO: remove this argument
                       logger=None):
    # Maybe it will be two databases: prot and nucl

    profiles_list = list()
    for btax_name in btax_dict:
        btax_info = BtaxInfo.load_from_dict(btax_dict[btax_name])
        try:
            profiles_list.extend(btax_info.repr_profiles.values())
        except (KeyError, AttributeError, TypeError):
            continue
    profiles_db_path = os.path.join(db_dir, profiles_db_name)
    if method.lower() == "hmmer":
        seq_profiles_db = SeqProfilesDB.build(profiles=profiles_list, name=profiles_db_path,
                                              method="hmmer", hmmer_inst_dir=hmmer_inst_dir, logger=logger)
    return profiles_db_path
