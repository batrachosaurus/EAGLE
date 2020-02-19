import json
import os
import multiprocessing as mp

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from eagle.btax_scanner import get_btax_name
from eagle.constants import conf_constants, eagle_logger, ORF_ALNS_DIR, ORF_TREES_DIR, DB_INFO_NAME
from eagle.lib.alignment import BlastHandler, construct_mult_aln, get_kaks_gmean
from eagle.lib.general import worker
from eagle.lib.seqs import get_orfs, load_fasta_to_dict, read_blast_out, parse_orf_id
from eagle.lib.orthan import explore_ortho_group
from eagledb.scheme import BtaxInfo, DBInfo


def explore_orfs_cmd():
    # This function will parse cmd input with argparse and run explore_orfs
    pass


def explore_orfs(in_fasta,
                 eagle_db,
                 annotation=None,
                 cu_table='default',
                 out_dir="",
                 min_orf_l=None,
                 btax_name=None,
                 num_threads=None,
                 btax_det_method=None,
                 config_path=None,
                 **kwargs):

    if config_path:
        conf_constants.update_by_config(config_path)
    if num_threads is not None:
        conf_constants.num_threads = int(num_threads)
    if min_orf_l is not None:
        conf_constants.min_orf_l = min_orf_l
    if btax_det_method is not None:
        conf_constants.btax_det_method = btax_det_method

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if kwargs.get("save_alignments", False) and not os.path.exists(os.path.join(out_dir, ORF_ALNS_DIR)):
        os.makedirs(os.path.join(out_dir, ORF_ALNS_DIR))
    if kwargs.get("save_trees", False) and not os.path.exists(os.path.join(out_dir, ORF_TREES_DIR)):
        os.makedirs(os.path.join(out_dir, ORF_TREES_DIR))

    if type(eagle_db) is str:
        if os.path.isdir(eagle_db):
            db_json = os.path.join(eagle_db, DB_INFO_NAME)
        else:
            db_json = eagle_db
        with open(db_json) as db_json_f:
            db_info = DBInfo.load_from_dict(json.load(db_json_f))
        # TODO: db_info.adjust_paths(os.path.dirname(db_json)) using db_info.from_root
    elif isinstance(eagle_db, dict):
        db_info = DBInfo.load_from_dict(eagle_db)
    else:
        eagle_logger.error("Unsupported type of value for 'eagle_db' argument")
        return
    with open(db_info.btax_json) as btax_dict_f:
        btax_dict = json.load(btax_dict_f)
    if btax_name is None:
        btax_name = get_btax_name(in_fasta,
                                  db_info.repr_profiles,
                                  btax_names=btax_dict.keys(),
                                  work_dir=out_dir,
                                  num_threads=conf_constants.num_threads,
                                  method=conf_constants.btax_det_method,
                                  hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                                  config_path=config_path)

    orfs_fasta_path = os.path.join(out_dir, os.path.basename(in_fasta)+".orfs")
    res_gtf_json = get_orfs(in_fasta_path=in_fasta,
                            out_fasta_path=orfs_fasta_path,
                            minsize=conf_constants.min_orf_l)
    blast_handler = BlastHandler(inst_dir=conf_constants.blast_inst_dir,
                                 config_path=config_path,
                                 logger=eagle_logger)

    if btax_name == "Unclassified":
        eagle_logger.warning("Basic taxon was not detected - cannot run further analysis")
    else:
        btax_info = BtaxInfo.load_from_dict(btax_dict[btax_name])
        eagle_logger.info("Basic taxon '%s' will be used for the sequence from %s" % (btax_name, in_fasta))
        tblastn_out_path = kwargs.get("tblastn_result_path", None)  # for debug and testing
        if tblastn_out_path is None:
            tblastn_out_path = os.path.join(out_dir, os.path.basename(in_fasta) + ".bl")
            blast_handler.run_blast_search(blast_type="tblastn",
                                           query=orfs_fasta_path,
                                           db=btax_info.blastdb,
                                           out=tblastn_out_path,
                                           num_threads=conf_constants.num_threads)
        res_gtf_json = analyze_tblastn_out(tblastn_out_path=tblastn_out_path,
                                           orfs_fasta_path=orfs_fasta_path,
                                           in_fasta=in_fasta,
                                           btax_data=btax_dict[btax_name],
                                           res_gtf_json=res_gtf_json,
                                           num_threads=conf_constants.num_threads,
                                           work_dir=out_dir,
                                           save_alignments=kwargs.get("save_alignments", False),
                                           save_trees=kwargs.get("save_trees", False))
    res_gtf_df = pd.DataFrame(res_gtf_json.values())
    res_gtf_df.sort_values("start", inplace=True)
    res_gtf_df = res_gtf_df[["seqid", "source", "type", "start", "end", "score", "strand", "frame", "attribute"]]
    res_gtf_df.to_csv(os.path.join(out_dir, os.path.basename(in_fasta)+".orfs.gtf"),
                      sep="\t", index=False, quotechar="'")


def _get_queries_btax(queries_scores_dict):
    queries_btax = dict()
    for query in queries_scores_dict.keys():
        try:
            queries_btax[query] = sorted(queries_scores_dict[query].items(), key=lambda x: x[1], reverse=True)[0][0]
        except IndexError:
            queries_btax[query] = "Unclassified"
    return queries_btax


def analyze_tblastn_out(tblastn_out_path,
                        orfs_fasta_path,
                        in_fasta,
                        btax_data,
                        res_gtf_json,
                        num_threads=conf_constants.num_threads,
                        work_dir="",
                        save_alignments=False,
                        save_trees=False):

    in_fasta_dict = load_fasta_to_dict(fasta_path=in_fasta)
    btax_info = BtaxInfo.load_from_dict(btax_data)
    orfs_stats = mp.Manager().dict()
    seq_ids_to_orgs = btax_info.fna_id
    tblatn_out_dict = read_blast_out(blast_out_path=tblastn_out_path)
    orfs_fasta_dict = load_fasta_to_dict(fasta_path=orfs_fasta_path)
    params_list = list()
    for seq_id in orfs_fasta_dict:
        chr_id, c1, c2, ori = parse_orf_id(seq_id)
        if ori == "+":
            orf_nucl_seq = in_fasta_dict[chr_id][c1-1: c2]
        else:
            orf_nucl_seq = str(Seq(in_fasta_dict[chr_id][c1-1: c2]).reverse_complement())
        params_list.append({
            "function": get_orf_stats,
            "orf_id": seq_id,
            "orf_homologs_prot": {seq_id: orfs_fasta_dict[seq_id]},
            "orf_homologs_nucl": {seq_id: orf_nucl_seq},
            "homologs_list": tblatn_out_dict[seq_id],
            "btax_data": btax_data,
            "seq_ids_to_orgs": seq_ids_to_orgs,
            "orfs_stats": orfs_stats,
            "work_dir": work_dir,
            "save_alignment": save_alignments,
            "save_tree": save_trees
        })

    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
    # list(map(worker, params_list))  # for debug
    pool.close()
    pool.join()

    tblatn_out_dict = None
    orfs_fasta_dict = None
    eagle_logger.info("ORFs stats calculated")
    for orf_id in orfs_stats.keys():
        try:
            res_gtf_json[orf_id]["attribute"] = json.dumps(orfs_stats[orf_id])
        except KeyError:
            pass
    orfs_stats = None
    return res_gtf_json


def get_orf_stats(orf_id,
                  orf_homologs_prot,
                  orf_homologs_nucl,
                  homologs_list,
                  btax_data,
                  seq_ids_to_orgs,
                  orfs_stats,
                  work_dir,
                  **kwargs):

    orf_stats = {
        "uniformity_std90": -1.0,
        "uniformity_std95": -1.0,
        "uniformity_std98": -1.0,
        "uniformity_std100": -1.0,
        "phylo_diff": -1.0,
        "Ka/Ks": -1.0,
        "Ks": -1.0,
        "representation": 0.0,
        "relative_mean_btax_dist": -1.0,
        "relative_median_btax_dist": -1.0,
        "relative_mean_ORF_dist": -1.0,
        "relative_median_ORF_dist": -1.0,
        "stops_per_seq_median": -1.0,
        "seqs_with_stops_fract": -1.0,
    }

    if len(homologs_list) < 3:
        orfs_stats[orf_id] = orf_stats
        eagle_logger.warning("A few homologs number for ORF '%s'" % orf_id)
        return
    btax_info = BtaxInfo.load_from_dict(btax_data)
    seq_ids_to_orgs[orf_id] = "Input_Organism_X"
    btax_fna = load_fasta_to_dict(fasta_path=btax_info.btax_fna)
    for hom in homologs_list:
        if hom["subj_start"] <= hom["subj_end"]:
            nucl_seq = Seq(btax_fna[hom["subj_id"]][hom["subj_start"] - 1:hom["subj_end"]])
        else:
            nucl_seq = Seq(btax_fna[hom["subj_id"]][hom["subj_end"] - 1: hom["subj_start"]]).reverse_complement()
        try:
            orf_homologs_prot[hom["subj_id"]] = str(nucl_seq.translate())
        except TranslationError:
            continue
        orf_homologs_nucl[hom["subj_id"]] = str(nucl_seq)
    btax_fna = None
    eagle_logger.info("got homologs sequences for ORF '%s'" % orf_id)
    orf_mult_aln = construct_mult_aln(seq_dict=orf_homologs_prot,
                                      aln_name=orf_id.replace("|:", "_")+"_aln",
                                      aln_type="prot",
                                      method="MUSCLE",
                                      tmp_dir=os.path.join(work_dir, orf_id.replace("|:", "_")+"_aln_tmp"),
                                      logger=eagle_logger)
    eagle_logger.info("got multiple alignment for ORF '%s' homologs" % orf_id)

    orf_mult_aln.nucl_seqs_dict = orf_homologs_nucl
    orf_mult_aln.seq_ids_to_orgs = seq_ids_to_orgs
    orf_mult_aln, ortho_stats, orf_homs_tree = explore_ortho_group(homologs_mult_aln=orf_mult_aln,
                                                                   ref_tree_newick=btax_info.ref_tree_newick,
                                                                   phylo_params={"method": "FastME",
                                                                                 "options": {"-T": 1}},
                                                                   ref_tree_full_names=btax_info.ref_tree_full_names,
                                                                   kaks_for_only_first=True,
                                                                   orgs_tree=True,
                                                                   get_conflicting_taxons=False,
                                                                   work_dir=work_dir,
                                                                   logger=eagle_logger)
    orf_stats["representation"] = float(len(orf_mult_aln) - 1) / float(len(btax_info.genomes))
    dist_matrix = orf_mult_aln.get_distance_matrix().replace_negative(inplace=False)
    btax_dist_matrix = dist_matrix[list(filter(lambda seq_name: seq_name != orf_id, dist_matrix.seq_names))]
    if btax_info.median_d > 0.0:
        orf_stats["relative_mean_btax_dist"] = btax_dist_matrix.mean_dist / btax_info.mean_d
        orf_stats["relative_mean_ORF_dist"] = dist_matrix[orf_id].mean() / btax_info.mean_d
    if btax_info.median_d > 0.0:
        orf_stats["relative_median_btax_dist"] = btax_dist_matrix.median_dist / btax_info.median_d
        orf_stats["relative_median_ORF_dist"] = dist_matrix[orf_id].median() / btax_info.median_d
    stops_stats = orf_mult_aln.stop_codons_stats()
    orf_stats["stops_per_seq_median"] = stops_stats["stops_per_seq_median"]
    orf_stats["seqs_with_stops_fract"] = stops_stats["seqs_with_stops_fract"]
    orf_stats["Ka/Ks"], orf_stats["Ks"] = get_kaks_gmean(ortho_stats.pop("Ka/Ks"), ortho_stats.pop("Ks"), stype="neg")
    orf_stats.update(ortho_stats)

    orfs_stats[orf_id] = orf_stats
    eagle_logger.info("got ORF '%s' stats" % orf_id)
    if kwargs.get("save_alignment", False):
        orf_mult_aln.dump_alignment(os.path.join(work_dir, ORF_ALNS_DIR, orf_mult_aln.aln_name+".fasta"))
    if kwargs.get("save_tree", False):
        orf_homs_tree.dump_tree(tree_path=os.path.join(work_dir, ORF_TREES_DIR, orf_homs_tree.tree_name+".nwk"))
