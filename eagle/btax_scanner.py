import os
from collections import defaultdict

from eagle.constants import conf_constants, PROFILES_SCAN_OUT
from eaglib.logging import eagle_logger


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

    query_scores_dict = defaultdict(lambda: defaultdict(float))
    if method.lower() == "hmmer":
        hmmer_handler = HmmerHandler(inst_dir=hmmer_inst_dir,
                                     tmp_dir=os.path.join(work_dir, "hmmer_tmp"),
                                     config_path=config_path,
                                     logger=eagle_logger)
        eagle_logger.info("hmmscan started")
        scan_out_path = hmmer_handler.run_hmmscan(profiles_db,
                                                  in_fasta,
                                                  num_threads=num_threads,
                                                  out_path=os.path.join(work_dir, PROFILES_SCAN_OUT),
                                                  shred_in_fasta=True)
        eagle_logger.info("hmmscan finished")
        scan_result_dict = hmmer_handler.read_hsr(scan_out_path, is_shred=True)

        for query in scan_result_dict:
            for profile_name in scan_result_dict[query]:
                btax_name = _parse_btax_name(profile_name, btax_names)
                score = float(scan_result_dict[query][profile_name]["score"])
                if btax_name and score >= min_score:
                    query_scores_dict[query][btax_name] += score

    btax_scores_dict = _aggregate_queries(query_scores_dict)
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
