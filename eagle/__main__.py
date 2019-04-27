import sys
import argparse

from eagle.constants import conf_constants
from eagle import explore_orfs, classify_orfs


def _parse_cmd_args(*args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-i",
                        "--in-fasta",
                        help="Path to input fasta file",
                        required=True)
    parser.add_argument("-db",
                        "--db-json",
                        help="Path to json with eagledb to use description",
                        required=True)
    parser.add_argument("-o",
                        "--out-dir",
                        help="Path to the directory for output",
                        required=False,
                        default="")
    parser.add_argument("-m",
                        "--mode",
                        help="Mode to run eagle: 'genome' - parses input fasta file as single genome "
                             "even if there not only one aequence; 'contigs' - parses input fasta files as "
                             "independent contigs (default: 'genome'). NOTE: sequences in one file can not be " 
                             "analyzed with different modes. If you need this split your input fasta into " 
                             "several files each of them could be analyzed using proper mode",
                        required=False,
                        default='genome')
    parser.add_argument("-btn",
                        "--btax-name",
                        help="The name of base taxon. If specified eagle will not scan the eagledb and "
                             "will work straight with this base taxon. Applicable only with 'genome' mode",
                        required=False,
                        default=None)
    parser.add_argument("-nt",
                        "--num-threads",
                        help="Number of threads",
                        required=False,
                        default=conf_constants.num_threads)
    parser.add_argument("-btd",
                        "--btax-det-method",
                        help="Method name to detect base taxon for input sequence (default: 'hmmer')",
                        required=False,
                        default="hmmer")
    parser.add_argument("-c",
                        "--config-path",
                        help="Path to a config file",
                        required=False,
                        default=None)

    cmd_args = parser.parse_args(args)
    if cmd_args.config_path:
        conf_constants.update_by_config(config_path=cmd_args.config_path)
        cmd_args.num_threads = conf_constants.num_threads
    return cmd_args.__dict__


def main():
    args_dict = _parse_cmd_args(*sys.argv[1:])
    explore_orfs(**args_dict)
    classify_orfs()  # currently not written


if __name__ == "__main__":
    main()
