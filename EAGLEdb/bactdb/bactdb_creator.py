import gzip
import io
import json
import multiprocessing as mp
import os
import pickle
import platform
import subprocess

import wget

from EAGLE.constants import EAGLE_logger, conf_constants
from EAGLE.lib.alignment import construct_mult_aln
from EAGLE.lib.general import worker, load_fasta_to_dict
from EAGLE.lib.phylo import build_tree_by_dist
from EAGLEdb.constants import BACTERIA_LIST_F_NAME, ANALYZED_BACTERIA_F_NAME, BACT_FAM_F_NAME, conf_constants_db
from EAGLEdb.lib import get_links_from_html


def get_bacteria_from_ncbi(refseq_bacteria_link="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria",
                           genbank_bacteria_link="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria",
                           bactdb_dir="EAGLEdb/bacteria",
                           num_threads=None,
                           first_bact=None,
                           last_bact=None,
                           analyzed_bacteria=ANALYZED_BACTERIA_F_NAME,
                           remove_bact_list_f=False,
                           config_path=None):

    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)
    if not num_threads:
        num_threads = conf_constants.num_threads
    else:
        conf_constants.num_threads = num_threads

    try:
        os.makedirs(bactdb_dir)
    except OSError:
        EAGLE_logger.info("bactdb directory exists")
    refseq_list = get_links_from_html(refseq_bacteria_link, num_threads=num_threads)
    genbank_list = get_links_from_html(genbank_bacteria_link, num_threads=num_threads)
    EAGLE_logger.info("20 first redseq bacteria: %s" % "; ".join(refseq_list[:20]))
    EAGLE_logger.info("20 first genbank bacteria: %s" % "; ".join(genbank_list[:20]))
    try:
        analyzed_bacteria = pickle.load(open(os.path.join(bactdb_dir, analyzed_bacteria), 'rb'))
    except IOError:
        analyzed_bacteria = mp.Manager().dict()
    bacteria_list_f_path = os.path.join(bactdb_dir, BACTERIA_LIST_F_NAME)
    bacteria_list_f = io.open(bacteria_list_f_path, 'w', newline="\n")
    bacteria_list_f.write(u"[\n")
    bacteria_list_f.close()
    n = 1
    i = 0
    j = 0
    proc_list = list()
    while i < len(refseq_list) or j < len(genbank_list):
        if first_bact and n < first_bact: continue
        if last_bact and n > last_bact: break
        if genbank_list[j] < refseq_list[i]:
            p = mp.Process(target=worker,
                           args=({'function': get_bacterium,
                                  'ncbi_db_link': genbank_bacteria_link,
                                  'bacterium_name': genbank_list[j],
                                  'analyzed_bacteria': analyzed_bacteria,
                                  'db_dir': bactdb_dir,
                                  'source_db': "genbank",
                                  'try_err_message': "%s is not prepared: " % genbank_list[j],
                                  'logger': EAGLE_logger},
                                 ))
            j += 1
        else:
            p = mp.Process(target=worker,
                           args=({'function': get_bacterium,
                                  'ncbi_db_link': refseq_bacteria_link,
                                  'bacterium_name': refseq_list[i],
                                  'analyzed_bacteria': analyzed_bacteria,
                                  'db_dir': bactdb_dir,
                                  'source_db': "refseq",
                                  'try_err_message': "%s is not prepared: " % refseq_list[i],
                                  'logger': EAGLE_logger},
                                 ))
            i += 1
        if genbank_list[j] == refseq_list[i-1]:
            j += 1
        p.start()
        proc_list.append(p)
        n += 1
        if n % num_threads == 0:
            for proc in proc_list:
                proc.join()
            proc_list = list()
    for proc in proc_list:
        proc.join()
    proc_list = None
    analyzed_bacteria_f = open(os.path.join(bactdb_dir, ANALYZED_BACTERIA_F_NAME), 'wb')
    pickle.dump(analyzed_bacteria, analyzed_bacteria_f)
    analyzed_bacteria_f.close()
    bacteria_list_f = io.open(bacteria_list_f_path, 'a', newline="\n")
    bacteria_list_f.write(u"  {}\n]")
    bacteria_list_f.close()
    return json.load(open(bacteria_list_f_path))


def get_bacterium(ncbi_db_link, bacterium_name, analyzed_bacteria, db_dir, source_db=None, **kwargs):
    bacterium_info = {"family": None,
                      "genus": None,
                      "species": None,
                      "strain": None,
                      "download_prefix": None,
                      "16S_rRNA_file": None,
                      "source_db": source_db,
                      "repr": False}
    bacterium_link = ncbi_db_link + "/" + bacterium_name
    EAGLE_logger.info('bacterium link: %s' % bacterium_link)
    if analyzed_bacteria.get(bacterium_name, None):
        return 0
    bacterium_list = get_links_from_html(bacterium_link)
    if "representative" in bacterium_list:
        next_page = bacterium_link + "/" + "representative"
        bacterium_info["repr"] = True
    else:
        next_page = bacterium_link + "/" + "latest_assembly_versions"
    assemblies_list = get_links_from_html(next_page)
    if not assemblies_list:
        EAGLE_logger.warning("Assemblies not loaded")
        return 0
    bacterium_prefix = (next_page + "/" + assemblies_list[-1] + "/" + assemblies_list[-1]).replace("https", "ftp")
    bacterium_info["download_prefix"] = bacterium_prefix
    download_bacterium_files(bacterium_prefix, ["_wgsmaster.gbff.gz", "_rna_from_genomic.fna.gz"], db_dir)
    tax_f_name = assemblies_list[-1] + "_wgsmaster.gbff.gz"
    if not os.path.exists(os.path.join(db_dir, tax_f_name)):
        tax_f_name = None
        download_bacterium_files(bacterium_prefix, "_genomic.gbff.gz", db_dir)
        tax_f_name = assemblies_list[-1] + "_genomic.gbff.gz"
    bacterium_info["family"], bacterium_info["genus"], bacterium_info["species"], bacterium_info["strain"] = \
        get_taxonomy(tax_f_name, db_dir)
    EAGLE_logger.info("got %s taxonomy" % bacterium_info["strain"])
    #if not os.path.exists():
    #
    bacterium_info["16S_rRNA_file"] = get_16S_fasta(assemblies_list[-1] + "_rna_from_genomic.fna.gz",
                                                    db_dir,
                                                    bacterium_info["strain"])
    EAGLE_logger.info("got %s 16S rRNA" % bacterium_info["strain"])
    f = io.open(os.path.join(db_dir, BACTERIA_LIST_F_NAME), 'a', newline="\n")
    f.write(unicode("  "+json.dumps(bacterium_info)+",\n"))
    f.close()
    analyzed_bacteria[bacterium_name] = True


def download_bacterium_files(bact_prefix, suffixes, download_dir="./"):
    # TODO: rename and move to EAGLEdb.lib
    if type(suffixes) is str:
        suffixes_list = [suffixes]
    else:
        suffixes_list = list(suffixes)
    for suffix in suffixes_list:
        file_link = None
        file_link = bact_prefix + suffix
        if platform.system() == 'Windows':
            try:
                wget.download(file_link, out=download_dir)
            except IOError:
                EAGLE_logger.warning("'%s' file has not been found" % file_link)
        else:
            subprocess.call("wget " + file_link + " -P " + download_dir + "/", shell=True)


def get_taxonomy(f_name, f_dir, remove_tax_f=True):
    family = None
    genus = None
    species = None
    strain = None
    org = False
    tax_list = []
    f_path = os.path.join(f_dir, f_name)
    f = gzip.open(f_path, 'rb')
    for line_ in f:
        line = None
        line = line_.decode("utf-8").strip()
        if not line: continue
        if line[:9] == "REFERENCE" or line[:7] == "COMMENT" or line[:8] == "FEATURES":
            family = get_family(tax_list, genus, species, strain)
            break
        if line[:8] == "ORGANISM":
            org = True
            line_list = line.split()
            genus = line_list[1]
            species = genus + "_" + line_list[2]
            strain = "_".join(line_list[1:])
        elif org:
            tax_list += list(prepare_tax_line(line))
    f.close()
    if remove_tax_f:
        os.remove(f_path)
    return family, genus, species, strain


def get_family(tax_list, g, sp, st):
    fam = None
    n = -1
    while -n <= len(tax_list):
        tax_u = None
        tax_u = tax_list[n].replace(" ", "_")
        n = n - 1
        if tax_u == st or tax_u == sp or tax_u == g:
            continue
        else:
            fam = tax_u
            break
    return fam


def prepare_tax_line(tax_line):
    tax_line_list = tax_line.split(";")
    for elm_ in tax_line_list:
        elm = None
        elm = elm_.strip(" .\t")
        if elm: yield elm


def get_16S_fasta(f_name, f_dir, strain, remove_rna_f=True):
    fasta_path = os.path.join(f_dir, strain + "_16S_rRNA.fasta")
    fasta_f = open(fasta_path, 'w')
    f_path = os.path.join(f_dir, f_name)
    rRNA = False
    seq_list = []
    f = gzip.open(f_path, 'rb')
    for line_ in f:
        line = None
        line = line_.decode("utf-8").strip()
        if not line: continue
        if line[0] == ">":
            if rRNA:
                fasta_f.write("".join(seq_list) + "\n")
                rRNA = False
            if "[product=16S ribosomal RNA]" in line:
                rRNA = True
                fasta_f.write(line + "\n")
        elif rRNA:
            seq_list.append(line)
    if rRNA:
        fasta_f.write("".join(seq_list) + "\n")
        rRNA = False
    fasta_f.close()
    f.close()
    if remove_rna_f:
        os.remove(f_path)
    return fasta_path


def get_families_dict(bacteria_list, db_dir, num_threads=None, only_repr=False, config_path=None):
    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)
    if not num_threads:
        num_threads = conf_constants.num_threads
    else:
        conf_constants.num_threads = num_threads
    if not only_repr:
        only_repr = conf_constants_db.only_repr
    else:
        conf_constants_db.only_repr = only_repr

    families_dict = dict()
    for bacterium in bacteria_list:
        bacterium_data = {"download_prefix": bacterium["download_prefix"],
                          "16S_rRNA_file": bacterium["16S_rRNA_file"],
                          "fna_file": None,
                          "source_db": bacterium["source_db"],
                          "repr": bacterium['repr']}
        if only_repr and not bacterium['repr']: continue
        if families_dict[bacterium['family']]:
            if families_dict[bacterium['family']][bacterium['genus']]:
                if families_dict[bacterium['family']][bacterium['genus']][bacterium['species']]:
                    families_dict[bacterium['family']][bacterium['genus']][bacterium['species']][bacterium['strain']] =\
                        bacterium_data
                else:
                    families_dict[bacterium['family']][bacterium['genus']][bacterium['species']] = \
                        {bacterium['strain']: bacterium_data}
            else:
                families_dict[bacterium['family']][bacterium['genus']] = \
                    {bacterium['species']:
                        {bacterium['strain']: bacterium_data}
                     }
        else:
            families_dict = \
                {bacterium['genus']:
                    {bacterium['species']:
                         {bacterium['strain']: bacterium_data}
                     }
                # "16S_rRNA_tree": None,
                # "WGS_tree": None,
                # "16S_rRNA_gtf": os.path.join(db_dir, bacterium['family']+"_16S_rRNA.gtf"),
                # "WGS_gtf": os.path.join(db_dir, bacterium['family']+"_WGS.gtf"),
                # "16S_rRNA_profile": None,
                # "WGS_profile": None,
                 }

    bact_fam_f_path = os.path.join(db_dir, BACT_FAM_F_NAME)
    prepare_families(families_dict, db_dir, num_threads=num_threads)

    return json.load(open(bact_fam_f_path))


def prepare_families(families_dict, db_dir, num_threads=4):
    bact_fam_f_path = os.path.join(db_dir, BACT_FAM_F_NAME)
    bact_fam_f = io.open(bact_fam_f_path, 'w', newline="\n")
    bact_fam_f.write(u"[\n")
    bact_fam_f.close()

    n = 0
    proc_list = list()
    for family in families_dict.keys():
        p = mp.Process(target=worker,
                       args=({'function': prepare_family,
                              'family_name': family,
                              'family_data': families_dict[family],
                              'db_dir': db_dir},
                             ))
        p.start()
        proc_list.append(p)
        n += 1
        if n % num_threads == 0:
            for proc in proc_list:
                proc.join()
            proc_list = list()
    for proc in proc_list:
        proc.join()
    proc_list = None

    bact_fam_f = io.open(bact_fam_f_path, 'a', newline="\n")
    bact_fam_f.write(u"  {}\n]")
    bact_fam_f.close()


def prepare_family(family_name, family_data, db_dir):
    rRNA_seqs_dict = dict()  # {seq_id: seq}
    ids_to_org_dict = dict()  # {seq_id: bacterium_name}
    for genus in family_data.keys():
        for species in family_data[genus].keys():
            for strain in family_data[genus][species].keys():
                bacterium_rRNA_dict = load_fasta_to_dict(family_data[genus][species][strain]["16S_rRNA_file"])
                for rRNA_id in bacterium_rRNA_dict.keys():
                    new_rRNA_id = None
                    new_rRNA_id = rRNA_id.split(" ")[0].split("|")[1]
                    ids_to_org_dict[new_rRNA_id] = strain
                    rRNA_seqs_dict[new_rRNA_id] = bacterium_rRNA_dict[rRNA_id]

    # reduce and modify ids
    # TODO: follows
    ### This section will be upgraded with my own alignment method but now MUSCLE and hmmer 3 are used
    tmp_fam_dir = os.path.join(db_dir, family_name+"_tmp")
    rRNA_aln = construct_mult_aln(seq_dict=rRNA_seqs_dict,
                                  method="MUSCLE",
                                  aln_type="nucl",
                                  aligner_inst_dir=conf_constants.muscle_inst_dir,
                                  hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                                  tmp_dir=tmp_fam_dir,
                                  logger=EAGLE_logger)
    rRNA_aln.remove_paralogs(ids_to_org_dict, method="min_dist", inplace=True)  # If I use my own alignment method: method="spec_pos"
    family_data["16S_rRNA_gtf"] = os.path.join(db_dir, family_name+"_16S_rRNA.gtf")
    family_data["16S_rRNA_fasta"] = os.path.join(db_dir, family_name+"_16S_rRNA.fasta")
    rRNA_aln.get_blocks_tsv(gtf_path=family_data["16S_rRNA_gtf"],
                            fasta_path=family_data["16S_rRNA_fasta"],
                            meta_dict=ids_to_org_dict)
    family_data["16S_rRNA_tree"] = build_tree_by_dist(rRNA_aln.get_distance_matrix(), tmp_dir=tmp_fam_dir).newick()
    family_data["16S_rRNA_profile"] = os.path.join(db_dir, family_name+".hmm")
    # profiles should not be here
    rRNA_aln.get_hmm_profile(method='hmmer', profile_path=family_data["16S_rRNA_profile"])  # hmmer will be replaced with my own method

    family_json_f = open(os.path.join(db_dir, family_name+".json"), 'w')
    json.dump(family_data, family_json_f)
    family_json_f.close()
