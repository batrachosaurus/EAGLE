import io
import json
import os
from collections import defaultdict

import pandas as pd

from EAGLE.constants import conf_constants, EAGLE_logger, PROFILES_SCAN_OUT
from EAGLE.lib.alignment import HmmerHandler, BlastHandler
from EAGLE.lib.general import filter_list
from EAGLE.lib.seqs import get_orfs


def explore_genes(in_fasta,
                  db_json,
                  out_dir="",
                  mode=None,
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

    with open(db_json) as db_json_f:
        db_info = json.load(db_json_f)
    btax_names = get_btax(in_fasta,
                          db_info["db_repr_profiles"],
                          btax_names=db_info.keys(),
                          working_dir=out_dir,
                          mode=conf_constants.mode,
                          num_threads=conf_constants.num_threads,
                          method=btax_det_method,
                          hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                          config_path=config_path)
    # for statistics
    resp_f = io.open("responsed_bacteria.json", 'a', newline="\n")
    resp_f.write(unicode("  "+json.dumps({in_fasta: btax_names.items()[0][1]})+",\n"))
    resp_f.close()
    # analysis continuation
    orfs_fasta_path = os.path.join(out_dir, os.path.basename(in_fasta)+".orfs")
    res_gtf_json = get_orfs(in_fasta_path=in_fasta,
                            out_fasta_path=orfs_fasta_path)
    blast_handler = BlastHandler(inst_dir=conf_constants.blast_inst_dir,
                                 config_path=config_path,
                                 logger=EAGLE_logger)
    if mode == "genome":
        fam_name = btax_names.items()[0][1]
        if fam_name == "Unclassified":
            EAGLE_logger.warning("The family was not detected - cannot run further analysis")
        else:
            EAGLE_logger.info("Family %s detected for sequence from %s" % (fam_name, in_fasta))
            tblastn_out_path = os.path.join(out_dir, os.path.basename(in_fasta) + ".bl")
            blast_handler.run_blast_search(blast_type="tblastn",
                                           query=orfs_fasta_path,
                                           db=db_info[fam_name]["blastdb"],
                                           out=tblastn_out_path)
            res_gtf_json = analyze_tblastn_out(tblastn_out_path, orfs_fasta_path, res_gtf_json)
    else:
        # TODO: write blast for contigs mode
        pass
    res_gtf_df = pd.DataFrame(res_gtf_json.values())
    res_gtf_df.sort_values("start", inplace=True)
    res_gtf_df.to_csv(os.path.join(out_dir, os.path.basename(in_fasta)+".gtf"), sep="\t", index=False)


def get_btax(in_fasta,
             profiles_db,
             btax_names,
             working_dir="",
             mode=conf_constants.mode,
             num_threads=conf_constants.num_threads,
             method="hmmer",
             hmmer_inst_dir=conf_constants.hmmer_inst_dir,
             config_path=None,
             remove_scan_out=True):

    if method.lower() == "hmmer":
        hmmer_handler = HmmerHandler(inst_dir=hmmer_inst_dir,
                                     tmp_dir=os.path.join(working_dir, "hmmer_tmp"),
                                     config_path=config_path,
                                     logger=EAGLE_logger)
        EAGLE_logger.info("hmmscan started")
        hmmer_handler.run_hmmscan(profiles_db,
                                  in_fasta,
                                  num_threads=num_threads,
                                  out_path=os.path.join(working_dir, PROFILES_SCAN_OUT))
        EAGLE_logger.info("hmmscan finished")
        queries_scores_dict = defaultdict(dict)
        query_scores_dict = defaultdict(float)
        lines_from_query = 0
        query = None
        profiles_scan_f = open(os.path.join(working_dir, PROFILES_SCAN_OUT))
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
            os.remove(os.path.join(working_dir, PROFILES_SCAN_OUT))

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


def analyze_tblastn_out(tblastn_out_path, orfs_fasta_path, res_gtf_json):
    
    return res_gtf_json
