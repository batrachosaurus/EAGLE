import argparse
import sys

from EAGLE.constants import conf_constants, EAGLE_logger
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
    parser.add_argument("-icustom",
                        "--input-table-custom",
                        help="Path to a table with custom genomes and their taxonomy (not implemented yet)",
                        required=False,
                        default=None)
    parser.add_argument("-btcp",
                        "--btax-class-profile",
                        help="The path to HMM profile of sequences that should be used for base taxons classification "
                             "while the db construction (not implemented yet)",
                        required=False,
                        default=None)
    parser.add_argument("-btrp",
                        "--btax-rep-profile",
                        help="The path to HMM profile of sequences that should be used for a base taxon response "
                             "while essential and advantageous genes exploration (not implemented yet)",
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
    return cmd_args.__dict__


def main():
    args_dict = _parse_cmd_args(*sys.argv[1:])
    if args_dict["dbtype"].lower() in ("bacteria", "bact", "b"):
        EAGLE_logger.info("'Bacteria' mode selected")
        create_bactdb(**args_dict)
    # metagenomes mode
    if args_dict["dbtype"].lower() in ("eukaryota", "eukaryot", "eukar", "euk", "eu", "e"):
        print("'Eukaryota' mode is not implemented yet")
    # mito mode
    if args_dict["dbtype"].lower() in ("other", "o"):
        print("'Other' mode is not implemented yet")
        if not args_dict["input_table_custom"] or not args_dict["btax_class_profile"]:
            print("'input-table-custom' and 'btax-class-profile' parameters required for 'Other' mode")
            return
    else:
        print("Unsupported mode '%s'" % args_dict["dbtype"])


if __name__ == "__main__":
    main()
