import sys
import argparse

from EAGLEdb.bactdb_creator import get_bacteria_from_ncbi, get_families_dict
from EAGLE.constants import conf_constants
from EAGLEdb.constants import DEFAULT_BACTDB_DIR


def create_bactdb(input_table_refseq=None,
                   input_table_genbank=None,
                   bactdb_dir=DEFAULT_BACTDB_DIR,
                   num_threads=None,
                   config_path=None):

    bacteria_list = get_bacteria_from_ncbi(refseq_bacteria_table=input_table_refseq,
                                           genbank_bacteria_table=input_table_genbank,
                                           num_threads=num_threads,
                                           bactdb_dir=bactdb_dir)
    families_dict = get_families_dict(bacteria_list=bacteria_list,
                                      num_threads=num_threads,
                                      db_dir=bactdb_dir,
                                      only_repr=True)
    pass


def main():

    pass


if __name__ == "__main__":
    main()
