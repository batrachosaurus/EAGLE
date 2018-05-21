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
    parser.add_argument("-db",
                        "--db-json",
                        help="Path to json with EAGLEdb to use description",
                        required=True)




def main():
    args_dict = _parse_cmd_args(*sys.argv[1:])
    explore_genes(**args_dict)


if __name__ == "__main__":
    main()
