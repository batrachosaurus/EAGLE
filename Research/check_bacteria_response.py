import os
import io
import sys
import json

import pandas


from EAGLE import explore_genes
from EAGLE.constants import EAGLE_logger
from EAGLE.lib.general import gunzip
from EAGLEdb.lib.db_creation import download_organism_files
from EAGLEdb.bactdb_creation import get_taxonomy


SUMMARY_TABLES = sys.argv[1].split(",")
DB_JSON = sys.argv[2]
WORKING_DIR = "check_bactdb"
NUM_THREADS = 6
PROCESSED_BACT_JSON = "processed_bact.json"

# USABLE_FAMILIES = json.load(open(DB_JSON)).keys()
USABLE_FAMILIES = [  # list of families that consists of at list 4 genera and 8 species
    "Erwiniaceae",
    "Rhodospirillaceae",
    "Helicobacteraceae",
    "Halomonadaceae",
    "Alteromonadaceae",
    "Sphingomonadaceae",
    "Rhodobacteraceae",
    "Pseudonocardiaceae",
    "Bradyrhizobiaceae",
    "Burkholderiaceae",
    "Vibrionaceae",
    "Oxalobacteraceae",
    "Comamonadaceae",
    "Enterobacteriaceae",
    "Alcaligenaceae",
    "Flavobacteriaceae",
    "Erythrobacteraceae",
    "Micrococcaceae",
    "Microbacteriaceae",
    "Planococcaceae",
    "Ectothiorhodospiraceae",
    "Spirochaetaceae",
    "Acetobacteraceae",
    "Peptococcaceae",
    "Morganellaceae",
    "Lachnospiraceae",
    "Xanthomonadaceae",
    "Cytophagaceae",
    "Bifidobacteriaceae",
    "Thermaceae",
]


def get_ncbi_links_list(summary_table):
    summary_df = pandas.read_csv(summary_table, header=1, sep="\t", dtype=str)
    ncbi_links_list = summary_df.apply(lambda df_row: df_row['ftp_path'], axis=1)
    return filter(None, ncbi_links_list)


ncbi_db_links = list()
for summary_table in SUMMARY_TABLES:
    ncbi_db_links.__iadd__(get_ncbi_links_list(summary_table))
ncbi_db_links = list(set(ncbi_db_links))

if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR)
processed_bact_f = io.open(PROCESSED_BACT_JSON, 'w', newline="\n")
processed_bact_f.write(u"[\n")
for ncbi_db_link in ncbi_db_links:
    assembly_id = None
    tax_f_name = None
    fna_f_path = None
    in_fasta = None
    bacterium_info = {"family": None,
                      "genus": None,
                      "species": None,
                      "strain": None}
    assembly_id = ncbi_db_link.split("/")[-1]
    tax_f_name = assembly_id + "_genomic.gbff.gz"
    download_prefix = (ncbi_db_link + "/" + assembly_id).replace("https", "ftp")
    download_organism_files(download_prefix, "_genomic.gbff.gz", download_dir=WORKING_DIR, logger=EAGLE_logger)
    if os.path.exists(os.path.join(WORKING_DIR, tax_f_name)):
        bacterium_info["family"], bacterium_info["genus"], bacterium_info["species"], bacterium_info["strain"] = \
            get_taxonomy(tax_f_name, f_dir=WORKING_DIR)
        if bacterium_info["family"] in USABLE_FAMILIES:
            download_organism_files(download_prefix, "_genomic.fna.gz", download_dir=WORKING_DIR, logger=EAGLE_logger)
            fna_f_path = os.path.join(WORKING_DIR, assembly_id+"_genomic.fna.gz")
            if os.path.exists(fna_f_path):
                in_fasta = fna_f_path[:-3]
                gunzip(in_path=fna_f_path, out_path=in_fasta)
                explore_genes(in_fasta=in_fasta, db_json=DB_JSON, out_dir=WORKING_DIR, num_threads=NUM_THREADS)
                processed_bact_f.write(unicode("  "+json.dumps({in_fasta: bacterium_info["family"]})+",\n"))
                os.remove(in_fasta)

processed_bact_f.write(u"  {}\n]")
processed_bact_f.close()
