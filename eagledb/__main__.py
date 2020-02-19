import argparse
import sys

from eagle.constants import conf_constants, eagle_logger
from eagledb import create_bactdb
from eagledb.constants import conf_constants_db


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
                        help="Path to a table with custom genomes and their taxonomy (not implemented yet). "
                             "The table consists of two necessary columns (1 - genome path; 2 - genome taxonomy) "
                             "and one optional column (3 - organism name)",
                        required=False,
                        default=None)
    parser.add_argument("-btl",
                        "--btax-level",
                        help="The taxonomic level to split input genomes into basic taxons "
                             "(1 - species, 2 - genus, 3 - family, etc)",
                        type=int,
                        required=False,
                        default=None)
    parser.add_argument("-btcp",
                        "--btax-class-profile",
                        help="The path to HMM profile of sequences that should be used for basic taxons classification "
                             "while the db construction (not implemented yet)",
                        required=False,
                        default=None)
    parser.add_argument("-btrp",
                        "--btax-rep-profile",
                        help="The path to HMM profile of sequences that should be used for a basic taxon response "
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
                        type=int,
                        required=False,
                        default=None)
    parser.add_argument("-po",
                        "--prepared-organisms",
                        help='Path to a json with organisms not to prepare listed. '
                             'Format as follows: {"org_name": true}',
                        required=False,
                        default=None)
    parser.add_argument("-poinf",
                        "--prepared-organisms-info",
                        help='Path to a json with info for organisms just prepared',
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
        eagle_logger.info("'Bacteria' mode selected")
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
