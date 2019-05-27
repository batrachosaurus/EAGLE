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
    parser.add_argument("-l",
                        "--min-orf-l",
                        help="Minimal length for ORF to analyze",
                        required=False)
    parser.add_argument("-btn",
                        "--btax-name",
                        help="The name of basic taxon. If specified, EAGLE will not scan the EAGLEdb profiles and "
                             "will work straight with this basic taxon",
                        required=False,
                        default=None)
    parser.add_argument("-nt",
                        "--num-threads",
                        help="Number of threads",
                        required=False,
                        default=conf_constants.num_threads)
    parser.add_argument("-btd",
                        "--btax-det-method",
                        help="Method name to detect basic taxon for input sequence (default: 'hmmer')",
                        required=False,
                        default="hmmer")
    parser.add_argument("-c",
                        "--config-path",
                        help="Path to a config file",
                        required=False,
                        default=None)
    parser.add_argument("-tbnr",
                        "--tblastn-result-path",
                        help="Path to tblastn result (outfmt 7) if it exists",
                        required=False,
                        default=None)
    parser.add_argument("-sa",
                        "--save-alignments",
                        help="Set it '1' if ORFs multiple alignments are needed to be saved",
                        required=False,
                        default=False)
    parser.add_argument("-st",
                        "--save-trees",
                        help="Set it '1' if ORFs phylogenetic trees are needed to be saved",
                        required=False,
                        default=False)

    cmd_args = parser.parse_args(args)
    if cmd_args.config_path:
        conf_constants.update_by_config(config_path=cmd_args.config_path)
        cmd_args.num_threads = conf_constants.num_threads
    if cmd_args.save_alignments:
        cmd_args.save_alignments = bool(int(cmd_args.save_alignments))
    if cmd_args.save_trees:
        cmd_args.save_trees = bool(int(cmd_args.save_trees))
    return cmd_args.__dict__


def main():
    args_dict = _parse_cmd_args(*sys.argv[1:])
    explore_orfs(**args_dict)
    classify_orfs()  # currently not written


if __name__ == "__main__":
    main()
