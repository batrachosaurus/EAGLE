import json

from EAGLE.constants import conf_constants, EAGLE_logger, PROFILES_SCAN_OUT
from EAGLE.lib.alignment import HmmerHandler


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
                         mode=conf_constants.mode,
                         num_threads=conf_constants.num_threads,
                         method=btax_det_method,
                         hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                         config_path=config_path)


def get_btax(in_fasta,
             profiles_db,
             mode=conf_constants.mode,
             num_threads=conf_constants.num_threads,
             method="hmmer",
             hmmer_inst_dir=conf_constants.hmmer_inst_dir,
             config_path=None,
             remove_scan_out=True):

    if method.lower() == "hmmer":
        hmmer_handler = HmmerHandler(inst_dir=hmmer_inst_dir, config_path=config_path, logger=EAGLE_logger)
        hmmer_handler.run_hmmscan(profiles_db, in_fasta, num_threads=num_threads, out_path=PROFILES_SCAN_OUT)

