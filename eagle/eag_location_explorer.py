import io
import json
import os
from collections import defaultdict
import multiprocessing as mp

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from eagle.constants import conf_constants, eagle_logger, PROFILES_SCAN_OUT
from eagle.lib.alignment import HmmerHandler, BlastHandler, MultAln, construct_mult_aln
from eagle.lib.phylo import PhyloTree, build_tree_by_dist, compare_trees, dump_tree_newick
from eagle.lib.general import filter_list, worker, reverse_dict
from eagle.lib.seqs import get_orfs, load_fasta_to_dict, read_blast_out
from eagledb.scheme import BtaxInfo


def explore_genes(in_fasta,
                  db_json,
                  out_dir="",
                  mode=None,
                  btax_name=None,
                  num_threads=None,
                  btax_det_method="hmmer",
                  config_path=None,
                  **kwargs):

    if config_path:
        conf_constants.update_by_config(config_path)
    if num_threads:
        conf_constants.num_threads = int(num_threads)
    if mode:
        conf_constants.mode = None
        conf_constants.mode = mode

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(db_json) as db_json_f:
        db_info = json.load(db_json_f)
    if btax_name is None or mode != "genome":
        btax_names = get_btax(in_fasta,
                              db_info["db_repr_profiles"],
                              btax_names=db_info.keys(),
                              work_dir=out_dir,
                              mode=conf_constants.mode,
                              num_threads=conf_constants.num_threads,
                              method=btax_det_method,
                              hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                              config_path=config_path)###
    else:
        btax_names = {"btax_name": btax_name}

    orfs_fasta_path = os.path.join(out_dir, os.path.basename(in_fasta)+".orfs")
    res_gtf_json = get_orfs(in_fasta_path=in_fasta,
                            out_fasta_path=orfs_fasta_path)
    blast_handler = BlastHandler(inst_dir=conf_constants.blast_inst_dir,
                                 config_path=config_path,
                                 logger=eagle_logger)
    if mode == "genome":
        if btax_name is None:
            btax_name = btax_names.items()[0][1]
        if btax_name == "Unclassified":
            eagle_logger.warning("The family was not detected - cannot run further analysis")
        else:
            eagle_logger.info("Family %s will be used for sequence from %s" % (btax_name, in_fasta))
            tblastn_out_path = os.path.join(out_dir, os.path.basename(in_fasta) + ".bl")
            blast_handler.run_blast_search(blast_type="tblastn",
                                           query=orfs_fasta_path,
                                           db=db_info[btax_name]["blastdb"],
                                           out=tblastn_out_path)
            res_gtf_json = analyze_tblastn_out(tblastn_out_path=tblastn_out_path,
                                               orfs_fasta_path=orfs_fasta_path,
                                               btax_data=db_info[btax_name],
                                               res_gtf_json=res_gtf_json,
                                               num_threads=conf_constants.num_threads,
                                               work_dir=out_dir)
    else:
        # TODO: write blast for contigs mode
        pass
    res_gtf_df = pd.DataFrame(res_gtf_json.values())
    res_gtf_df.sort_values("start", inplace=True)
    res_gtf_df = res_gtf_df[["seqid", "source", "type", "start", "end", "score", "strand", "frame", "attribute"]]
    res_gtf_df.to_csv(os.path.join(out_dir, os.path.basename(in_fasta)+".gtf"), sep="\t", index=False, quotechar="'")


def get_btax(in_fasta,
             profiles_db,
             btax_names,
             work_dir="",
             mode=conf_constants.mode,
             num_threads=conf_constants.num_threads,
             method="hmmer",
             hmmer_inst_dir=conf_constants.hmmer_inst_dir,
             config_path=None,
             remove_scan_out=True):

    if method.lower() == "hmmer":
        hmmer_handler = HmmerHandler(inst_dir=hmmer_inst_dir,
                                     tmp_dir=os.path.join(work_dir, "hmmer_tmp"),
                                     config_path=config_path,
                                     logger=eagle_logger)
        eagle_logger.info("hmmscan started")
        hmmer_handler.run_hmmscan(profiles_db,
                                  in_fasta,
                                  num_threads=num_threads,
                                  out_path=os.path.join(work_dir, PROFILES_SCAN_OUT))
        eagle_logger.info("hmmscan finished")
        queries_scores_dict = defaultdict(dict)
        query_scores_dict = defaultdict(float)
        lines_from_query = 0
        query = None
        profiles_scan_f = open(os.path.join(work_dir, PROFILES_SCAN_OUT))
        for line_ in profiles_scan_f:
            line = None
            line = line_.strip()
            if not line:
                continue
            if line[0: 6] == "Query:":
                line_list = filter_list(line.split())
                query = line_list[1]
            elif query and lines_from_query < 5:
                lines_from_query += 1
            elif line.startswith("Domain annotation for each model (and alignments):"):
                queries_scores_dict[query] = query_scores_dict
                query_scores_dict = defaultdict(float)
                query = None
                lines_from_query = 0
            elif query:
                try:
                    line_list = filter_list(line.split())
                    btax_name = _get_btax_name(line_list[8], btax_names)
                    if btax_name:
                        query_scores_dict[btax_name] += float(line_list[4])
                except:
                    pass
        if remove_scan_out:
            os.remove(os.path.join(work_dir, PROFILES_SCAN_OUT))

    if mode == "genome":
        queries_scores_dict = _aggregate_queries(in_fasta, queries_scores_dict)
    return _get_queries_btax(queries_scores_dict)


def _get_btax_name(profile_name, btax_names):
    for btax_name in btax_names:
        btax_name_list = btax_name.lower().split("_")
        if btax_name_list == profile_name.lower().split("_")[:len(btax_name_list)]:
            return btax_name


def _aggregate_queries(in_fasta, queries_scores_dict):
    if "." in in_fasta:
        aggr_key = ".".join(in_fasta.split(".")[: -1])
    else:
        aggr_key = in_fasta
    queries_scores_dict[aggr_key] = defaultdict(int)
    for query in queries_scores_dict.keys():
        if query == aggr_key:
            continue
        for btax_name in queries_scores_dict[query]:
            queries_scores_dict[aggr_key][btax_name] += queries_scores_dict[query][btax_name]
        queries_scores_dict.pop(query)
    return queries_scores_dict


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
                        btax_data,
                        res_gtf_json,
                        num_threads=conf_constants.num_threads,
                        work_dir=""):

    btax_info = BtaxInfo.load_from_dict(btax_data)
    orfs_stats = mp.Manager().dict()
    seq_ids_to_orgs = btax_info.fna_id
    tblatn_out_dict = read_blast_out(blast_out_path=tblastn_out_path)
    orfs_fasta_dict = load_fasta_to_dict(fasta_path=orfs_fasta_path)
    params_list = list()
    for seq_id in orfs_fasta_dict:
        params_list.append({
            "function": get_orf_stats,
            "orf_id": seq_id,
            "orf_homologs_seqs": {seq_id: orfs_fasta_dict[seq_id]},
            "homologs_list": tblatn_out_dict[seq_id],
            "btax_data": btax_data,
            "seq_ids_to_orgs": seq_ids_to_orgs,
            "orfs_stats": orfs_stats,
            "work_dir": work_dir,
        })

    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
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
                  orf_homologs_seqs,
                  homologs_list,
                  btax_data,
                  seq_ids_to_orgs,
                  orfs_stats,
                  work_dir,
                  **kwargs):

    orf_stats = dict()
    if len(homologs_list) < 3:
        orf_stats = {"uniformity_std": None, "phylo_diff": None}
        orfs_stats[orf_id] = orf_stats
        eagle_logger.warning("A few homologs number for ORF '%s'" % orf_id)
        return
    btax_info = BtaxInfo.load_from_dict(btax_data)
    seq_ids_to_orgs[orf_id] = "Input_Organism_X"
    btax_fna = load_fasta_to_dict(fasta_path=btax_info.btax_fna)
    for hom in homologs_list:
        try:
            if hom["subj_start"] <= hom["subj_end"]:
                orf_homologs_seqs[hom["subj_id"]] = \
                    str(Seq(btax_fna[hom["subj_id"]][hom["subj_start"]-1:hom["subj_end"]]).translate())
            else:
                orf_homologs_seqs[hom["subj_id"]] = \
                    str(Seq(btax_fna[hom["subj_id"]][hom["subj_end"]-1:
                                                     hom["subj_start"]]).reverse_complement().translate())
        except TranslationError:
            continue
    btax_fna = None
    eagle_logger.info("got homologs sequences for ORF '%s'" % orf_id)
    orf_mult_aln = construct_mult_aln(seq_dict=orf_homologs_seqs,
                                      aln_name=orf_id.replace("|:", "_")+"_aln",
                                      aln_type="prot",
                                      method="MUSCLE",
                                      tmp_dir=os.path.join(work_dir, orf_id.replace("|:", "_")+"_aln_tmp"),
                                      logger=eagle_logger)
    eagle_logger.info("got multiple alignment for ORF '%s' homologs" % orf_id)
    # Uniformity
    orf_mult_aln.remove_paralogs(seq_ids_to_orgs=seq_ids_to_orgs, method="min_dist", inplace=True)
    if len(orf_mult_aln.seqs) < 4:
        orf_stats = {"uniformity_std": None, "phylo_diff": None}
        orfs_stats[orf_id] = orf_stats
        eagle_logger.warning("A few homologs number for ORF '%s'" % orf_id)
        return
    orf_stats["uniformity_std"] = orf_mult_aln.improve_aln(inplace=False).estimate_uniformity(
        cons_thr=conf_constants.cons_thr,
        window_l=conf_constants.unif_window_l,
        windows_step=conf_constants.unif_windows_step
    )
    if np.isnan(orf_stats["uniformity_std"]):
        orf_stats["uniformity_std"] = None
    eagle_logger.info("got uniformity_std for ORF '%s'" % orf_id)
    # Ka/Ks
    # Phylo
    phylo_tmp_dir = os.path.join(work_dir, orf_id.replace("|:", "_")+"_phylo_tmp")
    try:
        del orf_mult_aln[orf_id]
    except KeyError:
        pass
    orf_homs_tree = build_tree_by_dist(
        dist_matrix=orf_mult_aln.get_distance_matrix(),
        method="FastME",
        full_seq_names=dict((short_id, seq_ids_to_orgs[full_id])
                            for short_id, full_id in orf_mult_aln.short_to_full_seq_names.items()),###
        tree_name=orf_id.replace("|:", "_")+"_tree",
        tmp_dir=phylo_tmp_dir,
        logger=eagle_logger
    )
    orf_homs_tree.set_full_names(inplace=True)
    btax_tree = PhyloTree.load_tree_from_str(
        tree_str=btax_info.ref_tree_newick,
        full_seq_names=btax_info.ref_tree_full_names,
        tree_name="btax_tree",
        tmp_dir=phylo_tmp_dir,
        logger=eagle_logger
    )
    btax_tree.set_full_names(inplace=True)
    orf_stats["phylo_diff"] = compare_trees(phylo_tree1=orf_homs_tree,
                                            phylo_tree2=btax_tree,
                                            method="Robinson-Foulds")
    eagle_logger.info("got phylo_diff for ORF '%s'" % orf_id)

    orfs_stats[orf_id] = orf_stats
    eagle_logger.info("got ORF '%s' stats" % orf_id)
