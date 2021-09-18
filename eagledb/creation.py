# each path for DB file should be stored as relative path from db_dir
# each time a DB file is used it should be accessed with os.path.join(db_dir, relative_f_path)
import argparse
from collections import Iterable


def create_cmd():
    parser = argparse.ArgumentParser()
    pass


def create(db_dir: str,
           genomes_table_path: str,
           btax_class_profiles: Iterable,
           btax_repr_profiles: Iterable,
           btax_level=None,
           num_threads=None,
           config_path=None,
           **kwargs):
    pass
