import io
import json
import shutil
from collections import defaultdict

from EAGLE.constants import conf_constants, EAGLE_logger, PROFILES_SCAN_OUT
from EAGLE.lib.alignment import HmmerHandler
from EAGLE.lib.general import filter_list


def explore_genes(in_fasta,
                  db_json,
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

    db_info = json.load(open(db_json))
    btax_names = get_btax(in_fasta,
                          db_info["db_repr_profiles"],
                          btax_names=db_info.keys(),
                          mode=conf_constants.mode,
                          num_threads=conf_constants.num_threads,
                          method=btax_det_method,
                          hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                          config_path=config_path)
    # for statistics
    resp_f = io.open("responsed_bacteria.json", 'a', newline="\n")
    resp_f.write(unicode("  "+json.dumps({in_fasta: btax_names.items()[0][1]})+",\n"))
    resp_f.close()


def get_btax(in_fasta,
             profiles_db,
             btax_names,
             mode=conf_constants.mode,
             num_threads=conf_constants.num_threads,
             method="hmmer",
             hmmer_inst_dir=conf_constants.hmmer_inst_dir,
             config_path=None,
             remove_scan_out=True):

    if method.lower() == "hmmer":
        hmmer_handler = HmmerHandler(inst_dir=hmmer_inst_dir, config_path=config_path, logger=EAGLE_logger)
        EAGLE_logger.info("hmmscan started")
        hmmer_handler.run_hmmscan(profiles_db, in_fasta, num_threads=num_threads, out_path=PROFILES_SCAN_OUT)
        EAGLE_logger.info("hmmscan finished")
        queries_scores_dict = defaultdict(dict)
        query_scores_dict = defaultdict(float)
        lines_from_query = 0
        query = None
        profiles_scan_f = open(PROFILES_SCAN_OUT)
        for line_ in profiles_scan_f:
            line = None
            line = line_.strip()
            if not line:
                continue
            if line[0: 6] == "Query:":
                line_list = filter_list(line.split())
                query = line_list[1]
            elif query and lines_from_query < 6:
                lines_from_query += 1
            elif line.startswith("Domain annotation for each model (and alignments):"):
                queries_scores_dict[query] = query_scores_dict
                query_scores_dict = dict()
                query = None
                lines_from_query = 0
            elif query:
                line_list = filter_list(line.split())
                btax_name = _get_btax_name(line_list[8], btax_names)
                if btax_name:
                    query_scores_dict[btax_name] += float(line_list[4])
        if remove_scan_out:
            shutil.rmtree(PROFILES_SCAN_OUT)

    if mode == "genome":
        queries_scores_dict = _aggregate_queries(in_fasta, queries_scores_dict)
    return _get_queries_btax(queries_scores_dict)


def _get_btax_name(profile_name, btax_names):
    for btax_name in btax_names:
        if btax_name.lower() in profile_name.lower():
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
            queries_scores_dict[aggr_key][btax_name] += queries_scores_dict.pop(query)[btax_name]
    return queries_scores_dict


def _get_queries_btax(queries_scores_dict):
    queries_btax = dict()
    for query in queries_scores_dict.keys():
        queries_btax[query] = sorted(queries_scores_dict.items(), key=lambda x: x[1], reverse=True)[0][0]
    return queries_btax
