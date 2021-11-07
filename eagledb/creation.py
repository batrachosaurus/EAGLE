# each path for DB file should be stored as relative path from db_dir
# each time a DB file is used it should be accessed with os.path.join(db_dir, relative_f_path)

import argparse
from collections import Iterable

import pandas as pd


def create_cmd():
    parser = argparse.ArgumentParser()
    pass
    # TODO: genomes_table_path -> genomes_table


def create(db_dir: str,
           genomes_table: pd.DataFrame,
           btax_class_profiles: Iterable,
           btax_repr_profiles: Iterable,
           btax_level=None,
           num_threads=None,
           config_path=None,
           **kwargs):
    """

    :param db_dir:
    :param genomes_table:
        id
        name
        taxonomy - fixed positions list of taxonomic units
        btc_seqs - list of paths to FASTA with sequences used for basic taxons classification
                   sequences can be joined into single file or distributed between several files
        genome - list of paths or links to the genome files (FASTA or archived FASTA)
                 sequences can be joined into single file or distributed between several files
    :param btax_class_profiles:
    :param btax_repr_profiles:
    :param btax_level:
    :param num_threads:
    :param config_path:
    :param kwargs:
    :return:
    """
    pass
