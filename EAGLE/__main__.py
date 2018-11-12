import sys
import argparse

from EAGLE.constants import conf_constants
from EAGLE import explore_genes


def _parse_cmd_args(*args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-i",
                        "--in-fasta",
                        help="Path to input fasta file",
                        required=True)
    parser.add_argument("-db",
                        "--db-json",
                        help="Path to json with EAGLEdb to use description",
                        required=True)
    parser.add_argument("-o",
                        "--out-dir",
                        help="Path to the directory for output",
                        required=False,
                        default="")
    parser.add_argument("-m",
                        "--mode",
                        help="Mode to run EAGLE: 'genome' - parses input fasta file as single genome "
                             "even if there not only one aequence; 'contigs' - parses input fasta files as "
                             "independent contigs (default: 'genome')",
                        required=False,
                        default='genome')
    parser.add_argument("-btn",
                        "--btax-name",
                        help="The name of base taxon. If specified EAGLE will not scan the EAGLEdb and "
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
                        "--config_path",
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
    explore_genes(**args_dict)


if __name__ == "__main__":
    main()
