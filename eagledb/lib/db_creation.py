# btax - base taxon

import os
import io
import subprocess
import urllib.request
import gzip
import shutil
from functools import reduce
import operator

from eaglib.constants import conf_constants as conf_constants_lib
from eaglib.logging import eagle_logger
from eaglib.general import join_files, gunzip, download_file
from eaglib.alignment import BlastDB, SeqProfilesDB
from eagledb.constants import PROFILES_DB_NAME
from eagledb.scheme import BtaxInfo, GenomeInfo


def get_links_from_html(html_link, n_tries=3, debug=False):
    n_t = 0
    links = dict()
    while n_t < n_tries:
        html_file = urllib.request.urlopen(html_link)
        _read_html_file_links(html_file=html_file.read().split("\n"), links=links, debug=debug)
        n_t += 1
    links_list = list(links.keys())
    links_list.sort()
    return links_list


def _read_html_file_links(html_file, links, **kwargs):
    for lines in html_file:
        line = None
        line = lines.strip()
        if not line: continue
        if "<a href" not in line: continue
        if "parent directory" in line.lower() or ".." in line: continue
        line_list = line.split("<a href")
        links[(line_list[1].split('">')[0].strip(' ="/'))] = True


def download_organism_files(org_prefix, suffixes, download_dir="./", logger=None):
    if type(suffixes) is str:
        suffixes_list = [suffixes]
    else:
        suffixes_list = list(suffixes)
    for suffix in suffixes_list:
        file_link = None
        file_link = org_prefix + suffix
        download_file(file_link=file_link, download_dir=download_dir, logger=logger)


def get_tree_from_dict(input_dict, stop_level=2, special_keys=tuple()):
    tree = dict()
    if stop_level == 1:
        return filter(lambda key: key not in special_keys, input_dict.keys())
    for key in input_dict.keys():
        if key in special_keys:
            continue
        tree[key] = get_tree_from_dict(input_dict[key], stop_level=stop_level-1)
    return tree


def clean_btax_data(btax_data, orgs_to_remain, stop_level=2, special_keys=tuple()):
    cleaned_btax_data = dict()
    if stop_level == 1:
        for btax_k in btax_data.keys():
            if btax_k in special_keys or btax_k in orgs_to_remain:
                cleaned_btax_data[btax_k] = btax_data[btax_k]
    else:
        for btax_k in btax_data.keys():
            if btax_k in special_keys:
                cleaned_btax_data[btax_k] = btax_data[btax_k]
                continue
            try:
                cleaned_btax_data[btax_k] = clean_btax_data(btax_data[btax_k], orgs_to_remain, stop_level=stop_level-1,
                                                            special_keys=special_keys)
            except AttributeError:
                continue
    return cleaned_btax_data


def download_btax_files(download_pref_dict, suffix, download_dir="./", logger=None):
    downloaded_fna = dict()
    for download_pref in download_pref_dict:
        downloaded_f_path = None
        downloaded_f_path = os.path.join(download_dir,
                                         download_pref.split("/")[-1] + suffix)
        download_organism_files(download_pref, suffix, download_dir, logger)
        if os.path.exists(downloaded_f_path):
            if downloaded_f_path[-3:] == ".gz":
                gunzip(in_path=downloaded_f_path,
                       out_path=downloaded_f_path[:-3])
                downloaded_fna[downloaded_f_path[:-3]] = download_pref_dict[download_pref]
            else:
                downloaded_fna[downloaded_f_path] = download_pref_dict[download_pref]
    return downloaded_fna


def get_from_btax_data(key, btax_data, key_path=None):
    if key_path is None:
        key_path = list()

    try:
        btax_data_keys = btax_data.keys()
        if key in btax_data_keys:
            return [(btax_data[key], key_path)]
        else:
            key_data_list = list()
            for btax_key in btax_data_keys:
                key_data_list.__iadd__(get_from_btax_data(key, btax_data[btax_key],
                                                          key_path=key_path.__add__([btax_key])))
            return filter(None, key_data_list)
    except AttributeError:
        return list()


def get_btax_fna(btax_genomes, btax_name, db_dir):
    chr_id_dict = dict()
    btax_fna_path = os.path.join(db_dir, btax_name + ".fasta")
    fna_to_orgs = dict()
    fna_to_download = dict()
    genome_id_to_org_name = dict()
    for btax_genome in btax_genomes:
        genome_info = GenomeInfo.load_from_dict(btax_genome)
        if genome_info.fna_path is None and genome_info.ncbi_download_prefix is not None:
            fna_to_download[genome_info.ncbi_download_prefix] = genome_info.genome_id
            genome_id_to_org_name[genome_info.genome_id] = genome_info.org_name
        else:
            fna_to_orgs[genome_info.fna_path] = genome_info.org_name
    if fna_to_download:
        downloaded_fna = download_btax_files(fna_to_download, suffix="_genomic.fna.gz", download_dir=db_dir)
        for fna_path, genome_id in downloaded_fna.items():
            fna_to_orgs[fna_path] = genome_id_to_org_name[genome_id]

    join_files(in_files_list=list(filter(None, fna_to_orgs.keys())),
               out_file_path=btax_fna_path,
               files_transform=transform_seq_id,
               **{"seq_id_dict": chr_id_dict, "fna_to_orgs": fna_to_orgs})
    return btax_fna_path, chr_id_dict, {v: k for k, v in downloaded_fna.items()}


def transform_seq_id(fna_f, seq_id_dict, fna_to_orgs, **kwargs):
    transf_fna_lines = list()
    for line_ in fna_f:
        line = None
        line = line_.strip()
        if not line:
            continue
        if line[0] == b">":
            seq_id = None
            seq_id = line[1:].split()[0]
            if seq_id not in seq_id_dict:
                transf_fna_lines.append(b">"+seq_id)
                seq_id_dict[seq_id] = fna_to_orgs[fna_f.name]
            else:
                i = 1
                while seq_id + "_" + str(i) in seq_id_dict:
                    i += 1
                transf_fna_lines.append(b">"+seq_id+b"_"+str(i).encode("utf-8"))
                seq_id_dict[seq_id+"_"+str(i)] = fna_to_orgs[fna_f.name]
        else:
            transf_fna_lines.append(line)
    return io.BytesIO(b"\n".join(transf_fna_lines)+b"\n")


def create_btax_blastdb(btax_fna_path, btax_name, db_dir, blast_inst_dir=None, logger=None):
    if blast_inst_dir is None:
        blast_inst_dir = conf_constants_lib.blast_inst_dir

    blastdb_dir = os.path.join(db_dir, btax_name+"_blastdb")
    if not os.path.exists(blastdb_dir):
        os.makedirs(blastdb_dir)
    blast_db_path = os.path.join(blastdb_dir, btax_name)
    blast_db = BlastDB.make_blastdb(in_seqs=btax_fna_path, dbtype="nucl", db_name=blast_db_path,
                                    blast_inst_dir=blast_inst_dir, logger=logger)
    return blast_db_path


def generate_btax_profiles(source, db_dir, btax_name, method="hmmer"):###
    # TODO: what is source?
    btax_profiles_dir = os.path.join(db_dir, btax_name+"_profiles")
    if not os.path.exists(btax_profiles_dir):
        os.makedirs(btax_profiles_dir)
    btax_profiles_dict = dict()
    for profile_name in source.keys():
        profile_path = None
        profile_path = os.path.join(btax_profiles_dir, btax_name+"_"+profile_name+".hmm")
        source[profile_name].get_hmm_profile(profile_path=profile_path, method=method)
        btax_profiles_dict[profile_name] = profile_path
    return btax_profiles_dict


def create_profiles_db(btax_dict,
                       db_dir,
                       profiles_db_name=PROFILES_DB_NAME,
                       method="hmmer",
                       hmmer_inst_dir="",
                       config_path=None,  # TODO: remove this argument
                       logger=None):
    # Maybe it will be two databases: prot and nucl

    profiles_list = list()
    for btax_name in btax_dict:
        btax_info = BtaxInfo.load_from_dict(btax_dict[btax_name])
        try:
            profiles_list.extend(btax_info.repr_profiles.values())
        except (KeyError, AttributeError, TypeError):
            continue
    profiles_db_path = os.path.join(db_dir, profiles_db_name)
    if method.lower() == "hmmer":
        seq_profiles_db = SeqProfilesDB.build(profiles=profiles_list, name=profiles_db_path,
                                              method="hmmer", hmmer_inst_dir=hmmer_inst_dir, logger=logger)
    return profiles_db_path
