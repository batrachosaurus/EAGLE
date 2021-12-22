# each path for DB file should be stored as relative path from db_dir
# each time a DB file is used it should be accessed with os.path.join(db_dir, relative_f_path)

import os
import argparse
from copy import deepcopy
from collections import Iterable, defaultdict

import pandas as pd

from eaglib._utils.logging import eagle_logger
from eaglib.seqs import SeqsDict, load_fasta_to_dict
from eaglib.alignment import SeqsProfileInfo
from eagledb.constants import conf_constants, conf_constants_lib
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

    btax_dict = get_btax_dict(db_dir=db_dir,
                              genomes_list=genomes_table.to_records(),
                              btax_level=conf_constants.btax_level,
                              btc_profiles=btax_class_profiles,
                              btr_profiles=btax_repr_profiles,
                              num_threads=conf_constants.num_threads)

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
                  btr_profiles=None,
                  **kwargs):

    btax_dict = defaultdict(BtaxInfo)
    btc_fasta_dict = defaultdict(dict)
    btc_info_dict = defaultdict(SeqsProfileInfo)
    for btc_profile_dict in btc_profiles:
        btc_profile_info = SeqsProfileInfo.load_from_dict(btc_profile_dict)
        btc_info_dict[btc_profile_info.name] = btc_profile_info

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
            btc_fasta_dict[btc_seq_id][genome_info.genome_id] = btc_seqs_dict[btc_seq_id]

    btc_dist_dict = dict()
    btc_aln_dict = dict()
    short_to_full_seq_names = dict()###
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

        # TODO: only the code from else block should be remained after moving 16S rRNA obtaining out from get_bacteria_from_ncbi
        if btc_profile_name == "16S_rRNA":
            btc_mult_aln.short_to_full_seq_names = \
                reduce_seq_names({re.sub("lcl\|(N(C|Z)_)?", "", seq_name): seq_name for seq_name in btc_mult_aln},
                                 num_letters=10, num_words=1)[0]
        else:
            btc_mult_aln.short_to_full_seq_names = short_to_full_seq_names.copy()###

        if btc_mult_aln.aln_type.lower() in SeqsDict.nucl_types:
            btc_mult_aln.dist_matr_options.update({"--dna": "p", "-T": conf_constants.num_threads, "-f": 6})
        btc_mult_aln.improve_aln(inplace=True)
        btc_dist_dict[btc_profile_name] = btc_mult_aln.get_distance_matrix()
        short_to_full_seq_names.update(btc_mult_aln.short_to_full_seq_names)###
        if kwargs.get("save_alignments", False):
            btc_mult_aln.dump(fname=os.path.join(db_dir, btc_mult_aln.aln_name + ".fasta"), format="fasta")
        btc_aln_dict[btc_profile_name] = deepcopy(btc_mult_aln)

    global_dist_matr = get_global_dist(btc_dist_dict, btc_profiles, seq_ids_to_orgs)
    global_dist_matr_path = os.path.join(db_dir, BACTERIA_GLOBAL_DIST_MATRIX)
    short_to_full_seq_names_path = os.path.join(db_dir, BACTERIA_SHORT_TO_FULL_ORG_NAMES)
    short_to_full_seq_names = global_dist_matr.dump(matrix_path=global_dist_matr_path, matr_format="phylip")
    with open(short_to_full_seq_names_path, "w") as short_to_full_org_names_f:
        json.dump(short_to_full_seq_names, short_to_full_org_names_f, indent=2)

    eagle_logger.info("base taxons standardisation started")
    btax_dict = standardize_btax(btax_dict=btax_dict, global_dist_matr=global_dist_matr)
    eagle_logger.info("base taxons standardisation finished")

    full_to_short_seq_names = {v: k for k, v in short_to_full_seq_names.items()}
    for btax_name in btax_dict:
        btax_orgs = set(GenomeInfo.load_from_dict(genome).org_name for genome in btax_dict[btax_name].genomes)
        if btr_profiles is not None:
            pass
        else:
            btax_dict[btax_name].mean_d = global_dist_matr[btax_orgs].mean_dist
            btax_dict[btax_name].median_d = global_dist_matr[btax_orgs].median_dist
            if len(btax_orgs) > 2:
                btax_dict[btax_name].ref_tree_newick = build_tree_by_dist(global_dist_matr[btax_orgs],
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
        btax_dict[btax_name].ref_tree_full_names = \
            {full_to_short_seq_names[btax_org]: btax_org for btax_org in btax_orgs}
        btax_dict[btax_name] = btax_dict[btax_name].get_json()
    return btax_dict


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
                                                                            db_dir=db_dir)
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
