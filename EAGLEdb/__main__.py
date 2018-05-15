import sys
import argparse
import json

from EAGLEdb.bactdb_creator import get_bacteria_from_ncbi, get_families_dict
from EAGLE.constants import conf_constants
from EAGLEdb.constants import DEFAULT_BACTDB_DIR, conf_constants_db, ANALYZED_BACTERIA_F_NAME
from EAGLEdb import join_bacteria_lists


def _parse_cmd_args(*args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-dbt",
                        "--dbtype",
                        help="The type of database to create (bacteria or archea or eukaryota)",
                        required=True)
    parser.add_argument("-irefseq",
                        "--input-table-refseq",
                        help="Path to a table with organisms to download from NCBI refseq",
                        required=False,
                        default=None)
    parser.add_argument("-igenbank",
                        "--input-table-genbank",
                        help="Path to a table with organisms to download from NCBI genbank",
                        required=False,
                        default=None)
    parser.add_argument("-d",
                        "--db-dir",
                        help="Path to a directory to collect database files",
                        required=False,
                        default=DEFAULT_BACTDB_DIR)
    parser.add_argument("-nt",
                        "--num-threads",
                        help="Threads number (can be set in config file)",
                        required=False,
                        default=conf_constants.num_threads)
    parser.add_argument("-anorgs",
                        "--analyzed-organisms",
                        help='Path to a json with organisms not to analyze listed. '
                             'Format as follows: {"org_name": true}',
                        required=False,
                        default=ANALYZED_BACTERIA_F_NAME)
    parser.add_argument("-c",
                        "--config-path",
                        help="Path to a config file",
                        required=False,
                        default=None)

    cmd_args = parser.parse_args(args)
    if cmd_args.config_path:
        conf_constants.update_by_config(config_path=cmd_args.config_path)
        conf_constants_db.update_by_config(config_path=cmd_args.config_path)
        cmd_args.num_threads = conf_constants.num_threads
    return cmd_args.__dict_


def create_bactdb(input_table_refseq=None,
                  input_table_genbank=None,
                  bactdb_dir=DEFAULT_BACTDB_DIR,
                  num_threads=None,
                  analyzed_bacteria_f_path=ANALYZED_BACTERIA_F_NAME,
                  analyzed_bacteria_info_list_path=None,
                  config_path=None,
                  **kwargs):

    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)

    bacteria_list = get_bacteria_from_ncbi(refseq_bacteria_table=input_table_refseq,
                                           genbank_bacteria_table=input_table_genbank,
                                           bactdb_dir=bactdb_dir,
                                           num_threads=num_threads,
                                           analyzed_bacteria_f_path=analyzed_bacteria_f_path)
    if analyzed_bacteria_info_list_path:
        bacteria_list = join_bacteria_lists(bacteria_list_1=bacteria_list,
                                            bacteria_list_2=json.load(open(analyzed_bacteria_info_list_path)))
    families_dict = get_families_dict(bacteria_list=bacteria_list,
                                      num_threads=num_threads,
                                      db_dir=bactdb_dir,
                                      only_repr=True)
    pass


def main():
    args_dict = _parse_cmd_args(*sys.argv[1:])
    if args_dict["dbtype"].lower() in ("bacteria", "bact", "b"):
        create_bactdb(**args_dict)


if __name__ == "__main__":
    main()
