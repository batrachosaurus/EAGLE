import os
import argparse
from urllib.request import urlretrieve, urlopen

import pandas as pd

from eagledb.creation import create


BACTERIA_SUMMARY_LINKS = {
    "refseq": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
    "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
}
DB_DIR = os.path.join("EAGLEdb", "bacteria")
# TODO: add reference profiles


def create_ncbi_cmd():
    parser = argparse.ArgumentParser()
    pass


def create_ncbi(db_name: str, db_dir=None, num_threads=None, config_path=None):
    # conf_constants update

    pd.read_table(urlopen(BACTERIA_SUMMARY_LINKS[db_name.lower()]), sep="\t", header=1, dtype=str)


    # download BACTERIA_SUMMARY_LINK
    # prepare BACTERIA_SUMMARY
    # create(...)
    pass


def collect_orgs_ncbi(orgs_table, db_dir=None):
    pass
