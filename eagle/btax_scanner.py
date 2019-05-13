import os
from collections import defaultdict

from eagle.constants import conf_constants, eagle_logger, PROFILES_SCAN_OUT
from eagle.lib.alignment import HmmerHandler
from eagle.lib.general import filter_list


def get_btax_name_cmd():
    # This function will parse cmd input with argparse and run get_btax_name
    pass


def get_btax_name(in_fasta,
                  profiles_db,
                  btax_names,
                  work_dir="",
                  num_threads=conf_constants.num_threads,
                  min_score=100.0,
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
            elif query and lines_from_query < 4:
                lines_from_query += 1
            elif line.startswith("Domain annotation for each model (and alignments):"):
                queries_scores_dict[query] = query_scores_dict.copy()
                query_scores_dict = defaultdict(float)
                query = None
                lines_from_query = 0
            elif query:
                try:
                    line_list = filter_list(line.split())
                    btax_name = _parse_btax_name(line_list[8], btax_names)
                    if btax_name and float(line_list[4]) >= min_score:
                        query_scores_dict[btax_name] += float(line_list[4])
                except:
                    pass
        if remove_scan_out:
            os.remove(os.path.join(work_dir, PROFILES_SCAN_OUT))

    btax_scores_dict = _aggregate_queries(queries_scores_dict)
    return sorted(btax_scores_dict.items(), key=lambda x: x[1], reverse=True)[0][0]


def _parse_btax_name(profile_name, btax_names):
    for btax_name in btax_names:
        btax_name_list = btax_name.lower().split("_")
        if btax_name_list == profile_name.lower().split("_")[:len(btax_name_list)]:
            return btax_name


def _aggregate_queries(queries_scores_dict):
    btax_scores_dict = defaultdict(float)
    for query in queries_scores_dict:
        for btax_name in queries_scores_dict[query]:
            btax_scores_dict[btax_name] += queries_scores_dict[query][btax_name]
    if not btax_scores_dict:
        btax_scores_dict["Unclassified"] = 1.0
    return btax_scores_dict
