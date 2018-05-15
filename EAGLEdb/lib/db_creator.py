import os
import io
import multiprocessing as mp
import platform
import subprocess
import urllib2
import gzip
import shutil

import wget

from EAGLE.lib.general import worker, join_files
from EAGLE.lib.alignment import BlastHandler


def get_links_from_html(html_link, num_threads=1, n_tries=3, debug=False):
    n_t = 0
    links = mp.Manager().dict()
    params_list = []
    while n_t < n_tries:
        html_file = urllib2.urlopen(html_link)
        params_list.append({'function': _read_html_file_links,
                            'html_file': html_file.read().split("\n"),
                            'links': links,
                            'debug': debug})
        n_t += 1
    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
    pool.close()
    pool.join()
    links_list = links.keys()
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
    # TODO: rename and move to EAGLEdb.lib
    if type(suffixes) is str:
        suffixes_list = [suffixes]
    else:
        suffixes_list = list(suffixes)
    for suffix in suffixes_list:
        file_link = None
        file_link = org_prefix + suffix
        if platform.system() == 'Windows':
            try:
                wget.download(file_link, out=download_dir)
            except IOError:
                if logger:
                    logger.warning("'%s' file has not been found" % file_link)
                else:
                    print("'%s' file has not been found" % file_link)
        else:
            subprocess.call("wget " + file_link + " -P " + download_dir + "/", shell=True)


def get_tree_from_dict(input_dict, stop_level=2):
    tree = dict()
    if stop_level == 1:
        return input_dict.keys()
    for key in input_dict.keys():
        tree[key] = get_tree_from_dict(input_dict[key], stop_level=stop_level-1)
    return tree


def clean_btax_data(btax_data, orgs_to_remain, stop_level=2):
    cleaned_btax_data = dict()
    if stop_level == 1:
        return dict(filter(None, map(lambda btax_k: (btax_k, btax_data[btax_k]) if btax_k in orgs_to_remain else None,
                                     btax_data.keys())))
    for btax_k in  btax_data.keys():
        cleaned_btax_data[btax_k] = clean_btax_data(btax_data[btax_k], orgs_to_remain, stop_level=stop_level-1)
    return cleaned_btax_data


def download_btax_files(key_prefix_pairs, btax_data, download_dir="./", logger=None):
    for key in key_prefix_pairs.keys():
        downloaded_f_path = None
        downloaded_f_path = os.path.join(download_dir,
                                         btax_data["download_prefix"].split("/")[-1] + key_prefix_pairs[key])
        download_organism_files(btax_data["download_prefix"], key_prefix_pairs[key], download_dir, logger)
        if os.path.exists(downloaded_f_path):
            if downloaded_f_path[-3:] == ".gz":
                btax_data[key] = downloaded_f_path[:-3]
                with gzip.open(downloaded_f_path, 'rb') as downloaded_f_gz, \
                        io.open(btax_data[key], 'wb', newline=b'\n') as downloaded_f:
                    shutil.copyfileobj(downloaded_f_gz, downloaded_f)
                    downloaded_f.close()
                os.remove(downloaded_f_path)
            else:
                btax_data[key] = downloaded_f_path
        else:
            btax_data[key] = None
    return btax_data


def create_btax_blastdb(btax_data, btax_name, db_dir, blast_inst_dir="", logger=None):
    blastdb_dir = os.path.join(db_dir, btax_name+"_blastdb")
    if not os.path.exists(blastdb_dir):
        os.makedirs(blastdb_dir)
    fna_list = get_from_btax_data("fna_file", btax_data)
    btax_fna_path = os.path.join(db_dir, btax_name+".fasta")
    blast_db_path = os.path.join(blastdb_dir, btax_name)
    join_files(fna_list, btax_fna_path)
    blast_handler = BlastHandler(inst_dir=blast_inst_dir, logger=logger)
    blast_handler.make_blastdb(btax_fna_path, dbtype="nucl", db_name=blast_db_path)
    return blast_db_path


def get_from_btax_data(key, btax_data):
    try:
        return [btax_data[key]]
    except KeyError:
        key_data_list = list()
        for btax_key in btax_data:
            key_data_list.__iadd__(get_from_btax_data(key, btax_data[btax_key]))
        return key_data_list


def generate_btax_profile(source, db_dir, btax_name, method="hmmer"):
    btax_profiles_tmp_dir = os.path.join(db_dir, btax_name+"_profiles_tmp")
    if not os.path.exists(btax_profiles_tmp_dir):
        os.makedirs(btax_profiles_tmp_dir)
    aln_profiles_list = list()
    for aln_key in source.keys:
        aln_profile_path = None
        aln_profile_path = os.path.join(btax_profiles_tmp_dir, aln_key+".hmm")
        source[aln_key].get_hmm_profile(profile_path=aln_profile_path, method=method)
        aln_profiles_list.append(aln_profile_path)
    btax_profile_path = os.path.join(db_dir, btax_name+".hmm")
    join_files(aln_profiles_list, btax_profile_path)
    shutil.rmtree(btax_profiles_tmp_dir)
    return btax_profile_path
