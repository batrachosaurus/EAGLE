import os
import argparse


BACTERIA_REFSEQ_SUMMARY_LINK = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
BACTERIA_GENBANK_SUMMARY_LINK = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
DB_DIR = os.path.join("EAGLEdb", "bacteria")
# TODO: add reference profiles


def create_refseq_cmd():
    parser = argparse.ArgumentParser()
    pass


def create_refseq(db_dir=None, num_threads=None, config_path=None):
    # conf_constants update

    # download BACTERIA_REFSEQ_SUMMARY_LINK
    # prepare BACTERIA_REFSEQ_SUMMARY
    # create_ncbi(...)
    pass
