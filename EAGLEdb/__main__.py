import argparse
import sys

from EAGLE.constants import conf_constants
from EAGLEdb import create_bactdb
from EAGLEdb.constants import conf_constants_db


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
                        default=None)
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
                        default=None)
    parser.add_argument("-anorgsinf",
                        "--analyzed-organisms-info",
                        help='Path to a json with info for organisms just analyzed',
                        required=False,
                        default=None)
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


def main():
    args_dict = _parse_cmd_args(*sys.argv[1:])
    if args_dict["dbtype"].lower() in ("bacteria", "bact", "b"):
        create_bactdb(**args_dict)


if __name__ == "__main__":
    main()
