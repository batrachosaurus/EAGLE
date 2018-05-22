import json
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
    btax_name = get_btax(in_fasta,
                         db_info["db_repr_profiles"],
                         btax_names=db_info.keys(),
                         mode=conf_constants.mode,
                         num_threads=conf_constants.num_threads,
                         method=btax_det_method,
                         hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                         config_path=config_path)


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
        hmmer_handler.run_hmmscan(profiles_db, in_fasta, num_threads=num_threads, out_path=PROFILES_SCAN_OUT)
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
            else:
                line_list = filter_list(line.split())
                btax_name = _get_btax_names(line_list[8], btax_names)
                if btax_name:
                    query_scores_dict[btax_name] += float(line_list[4])

    if mode == "genome":
        pass
    else:
        pass


def _get_btax_names(profile_name, btax_names):
    pass
