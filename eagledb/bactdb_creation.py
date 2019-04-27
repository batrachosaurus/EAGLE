import gzip
import io
import json
import os
from copy import deepcopy
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pandas

from eagle.constants import eagle_logger, conf_constants
from eagle.lib.alignment import construct_mult_aln, DistanceMatrix
from eagle.lib.general import worker, get_un_fix, bool_from_str
from eagle.lib.phylo import build_tree_by_dist
from eagle.lib.seqs import load_fasta_to_dict, reduce_seq_names
from eagledb import join_genomes_lists
from eagledb.constants import BACTERIA_LIST_F_NAME, PREPARED_BACTERIA_F_NAME, BACT_FAM_F_NAME, conf_constants_db, \
    DEFAULT_REFSEQ_BACTERIA_TABLE, DEFAULT_GENBANK_BACTERIA_TABLE, DEFAULT_BACTDB_DIR, PROFILES_DB_NAME, \
    BACTERIA_GLOBAL_DIST_MATRIX, BACTERIA_SHORT_TO_FULL_ORG_NAMES, BTAX_JSON_NAME, DB_INFO_NAME
from eagledb.lib.db_creation import download_organism_files, clean_btax_data, download_btax_files, \
    create_btax_blastdb, generate_btax_profile, create_profiles_db, get_btax_fna
from eagledb.scheme import GenomeInfo, SeqProfileInfo, BtaxInfo, DBInfo


def get_bacteria_from_ncbi(refseq_bacteria_table=None,
                           genbank_bacteria_table=None,
                           bactdb_dir=DEFAULT_BACTDB_DIR,
                           num_threads=None,
                           first_bact=None,
                           last_bact=None,
                           prepared_bacteria_f_path=PREPARED_BACTERIA_F_NAME,
                           remove_bact_list_f=False,
                           config_path=None):

    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)
    if not num_threads:
        num_threads = conf_constants.num_threads
    else:
        conf_constants.num_threads = num_threads
    if refseq_bacteria_table is None and genbank_bacteria_table is None:
        refseq_bacteria_table = DEFAULT_REFSEQ_BACTERIA_TABLE
        genbank_bacteria_table = DEFAULT_GENBANK_BACTERIA_TABLE
    try:
        os.makedirs(bactdb_dir)
    except OSError:
        eagle_logger.warning("bactdb directory exists")
    prepared_bacteria = mp.Manager().dict()
    if os.path.exists(prepared_bacteria_f_path):
        eagle_logger.info("loading prepared bacteria from '%s'" % prepared_bacteria_f_path)
        prepared_bacteria_f = open(prepared_bacteria_f_path)
        prepared_bacteria.update(json.load(prepared_bacteria_f))
        prepared_bacteria_f.close()
        eagle_logger.info("prepared bacteria loaded")
    bacteria_list_f_path = os.path.join(bactdb_dir, BACTERIA_LIST_F_NAME)
    bacteria_list_f = io.open(bacteria_list_f_path, 'w', newline="\n")
    bacteria_list_f.write(u"[\n")
    bacteria_list_f.close()
    refseq_df = pandas.read_csv(refseq_bacteria_table, sep="\t", dtype=str).sort_values(by="ncbi_link")
    genbank_df = pandas.read_csv(genbank_bacteria_table, sep="\t", dtype=str).sort_values(by="ncbi_link")
    n = 1
    i = 0
    j = 0
    params_list = list()
    while i < refseq_df.shape[0] or j < genbank_df.shape[0]:
        if first_bact and n < first_bact:
            n += 1
            continue
        if last_bact and n > last_bact: break
        if i >= refseq_df.shape[0] or j >= genbank_df.shape[0]:
            if i >= refseq_df.shape[0]:
                params_list.append({
                    'function': get_bacterium,
                    'prepared_bacteria': prepared_bacteria,
                    'logger_name': eagle_logger.name,
                    'ncbi_db_link': genbank_df.iloc[j]["ncbi_link"],
                    'bacterium_name': genbank_df.iloc[j]["org_name"],
                    'is_repr': bool_from_str(genbank_df.iloc[j]["repr"]),
                    'db_dir': bactdb_dir,
                    'source_db': "genbank",
                    'try_err_message': "%s is not prepared: " % genbank_df.iloc[j]["org_name"],
                })
                j += 1
            else:
                params_list.append({
                    'function': get_bacterium,
                    'prepared_bacteria': prepared_bacteria,
                    'logger_name': eagle_logger.name,
                    'ncbi_db_link': refseq_df.iloc[i]["ncbi_link"],
                    'bacterium_name': refseq_df.iloc[i]["org_name"],
                    'is_repr': bool_from_str(refseq_df.iloc[i]["repr"]),
                    'db_dir': bactdb_dir,
                    'source_db': "refseq",
                    'try_err_message': "%s is not prepared: " % refseq_df.iloc[i]["org_name"],
                })
                i += 1
        else:
            if genbank_df.iloc[j]["ncbi_link"].replace("GCA", "GCF") < refseq_df.iloc[i]["ncbi_link"]:
                params_list.append({
                    'function': get_bacterium,
                    'prepared_bacteria': prepared_bacteria,
                    'logger_name': eagle_logger.name,
                    'ncbi_db_link': genbank_df.iloc[j]["ncbi_link"],
                    'bacterium_name': genbank_df.iloc[j]["org_name"],
                    'is_repr': bool_from_str(genbank_df.iloc[j]["repr"]),
                    'db_dir': bactdb_dir,
                    'source_db': "genbank",
                    'try_err_message': "%s is not prepared: " % genbank_df.iloc[j]["org_name"],
                })
                j += 1
            else:
                params_list.append({
                    'function': get_bacterium,
                    'prepared_bacteria': prepared_bacteria,
                    'logger_name': eagle_logger.name,
                    'ncbi_db_link': refseq_df.iloc[i]["ncbi_link"],
                    'bacterium_name': refseq_df.iloc[i]["org_name"],
                    'is_repr': bool_from_str(refseq_df.iloc[i]["repr"]),
                    'db_dir': bactdb_dir,
                    'source_db': "refseq",
                    'try_err_message': "%s is not prepared: " % refseq_df.iloc[i]["org_name"],
                })
                i += 1
            if genbank_df.iloc[j]["ncbi_link"].replace("GCA", "GCF") == refseq_df.iloc[i-1]["ncbi_link"]:
                j += 1
        n += 1
    eagle_logger.info("got download links for %s bacteria" % len(params_list))
    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
    pool.close()
    pool.join()
    prepared_bacteria_f = open(os.path.join(bactdb_dir, PREPARED_BACTERIA_F_NAME), "w")
    json.dump(dict(prepared_bacteria), prepared_bacteria_f)
    prepared_bacteria_f.close()
    bacteria_list_f = io.open(bacteria_list_f_path, 'a', newline="\n")
    bacteria_list_f.write(u"  {}\n]")
    bacteria_list_f.close()
    with open(bacteria_list_f_path) as bacteria_list_f:
        return json.load(bacteria_list_f)


def get_bacterium(ncbi_db_link, bacterium_name, is_repr, prepared_bacteria, db_dir, source_db=None, **kwargs):
    if conf_constants_db.only_repr and not is_repr:
        return
    assembly_id = ncbi_db_link.split("/")[-1]
    bacterium_info = GenomeInfo(genome_id=assembly_id,
                                ncbi_download_prefix=(ncbi_db_link+"/"+assembly_id).replace("https", "ftp"),
                                source_db=source_db,
                                is_repr=is_repr)

    eagle_logger.info("%s getting started" % bacterium_name)
    eagle_logger.info('bacterium link: %s' % bacterium_info.ncbi_download_prefix)
    if prepared_bacteria.get(bacterium_info.genome_id, None):
        eagle_logger.info("%s is in prepared bacteria" % bacterium_name)
        return
    download_organism_files(bacterium_info.ncbi_download_prefix,
                            ["_wgsmaster.gbff.gz", "_rna_from_genomic.fna.gz"],
                            db_dir,
                            logger=eagle_logger)
    tax_f_name = assembly_id + "_wgsmaster.gbff.gz"
    if not os.path.exists(os.path.join(db_dir, tax_f_name)):
        tax_f_name = None
        download_organism_files(bacterium_info.ncbi_download_prefix, "_genomic.gbff.gz", db_dir, logger=eagle_logger)
        tax_f_name = assembly_id + "_genomic.gbff.gz"
    bacterium_info.taxonomy, bacterium_info.org_name = get_taxonomy(tax_f_name, db_dir)
    eagle_logger.info("got %s taxonomy" % bacterium_info.org_name)

    # TODO: no need to obtain 16S rRNA during this stage
    bacterium_info.btc_seqs_fasta, seq_id_list = get_16S_fasta(assembly_id + "_rna_from_genomic.fna.gz",
                                                               db_dir,
                                                               bacterium_info.genome_id)
    bacterium_info.btc_seqs_id = {seq_id: "16S_rRNA" for seq_id in seq_id_list}
    eagle_logger.info("got %s 16S rRNA" % bacterium_info.org_name)
    f = io.open(os.path.join(db_dir, BACTERIA_LIST_F_NAME), 'a', newline="\n")
    f.write(unicode("  "+json.dumps(bacterium_info.get_json())+",\n"))
    f.close()
    prepared_bacteria[bacterium_info.genome_id] = True


def get_taxonomy(f_name, f_dir, remove_tax_f=True):
    prepared_taxonomy = list()
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
            prepared_taxonomy = prepare_tax_list(tax_list, genus, species, strain)
            break
        if line[:8] == "ORGANISM":
            org = True
            line_list = line.split()
            genus = line_list[1]
            species = genus + "_" + line_list[2]
            strain = "_".join(line_list[1:])
        elif org:
            tax_list.extend(list(prepare_tax_line(line)))
    f.close()
    if remove_tax_f:
        os.remove(f_path)
    return prepared_taxonomy, strain


def prepare_tax_list(tax_list, genus, species, strain):
    prepared_tax_list = list()
    genus_match = False
    species_match = False
    strain_match = False
    _aceae_seen = False
    after_aceae = False
    for tax in tax_list:
        if tax.replace(" ", "_") == strain:
            strain_match = True
            break
        elif tax.replace(" ", "_") == species:
            species_match = True
            if not prepared_tax_list:
                prepared_tax_list = ["Unclassified", genus, species]
            elif _aceae_seen and not after_aceae:
                prepared_tax_list.extend([genus, species])
                after_aceae = True
            else:
                prepared_tax_list.append(species)
            break
        elif tax == genus:
            genus_match = True
            prepared_tax_list.extend([genus, species])
            if _aceae_seen:
                after_aceae = True
            break
        else:
            prepared_tax_list.append(tax)
            if _aceae_seen:
                after_aceae = True
                prepared_tax_list.append(species)
                break
            elif tax[-5:] == "aceae":
                _aceae_seen = True

    if _aceae_seen and not after_aceae:
        prepared_tax_list.extend([genus, species])
    elif len(prepared_tax_list) < 2 or (not genus_match and not species_match and not strain_match and not _aceae_seen):
        prepared_tax_list = ["Unclassified", genus, species]
    return prepared_tax_list


def get_family(tax_list, g, sp, st):  # currently is not used
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
        if elm:
            yield elm.replace(" ", "_")


def get_16S_fasta(f_name, f_dir, genome_id, remove_rna_f=True, max_l=3000):
    fasta_path = os.path.join(f_dir, genome_id + "_16S_rRNA.fasta")
    fasta_f = open(fasta_path, 'w')
    f_path = os.path.join(f_dir, f_name)
    title = None
    seq_list = []
    seq_id_list = list()
    f = gzip.open(f_path, 'rb')
    for line_ in f:
        line = None
        line = line_.decode("utf-8").strip()
        if not line: continue
        if line[0] == ">":
            if title is not None:
                seq = "".join(seq_list)
                if len(seq) <= max_l:
                    fasta_f.write(">" + title + "\n")
                    fasta_f.write(seq + "\n")
                else:
                    del seq_id_list[-1]
                seq = None
                title = None
            if "[product=16S ribosomal RNA]" in line:
                title = line[1:]
                seq_id_list.append(title)
        elif title is not None:
            seq_list.append(line)
    if title is not None:
        seq = "".join(seq_list)
        if len(seq) <= max_l:
            fasta_f.write(">" + title + "\n")
            fasta_f.write(seq + "\n")
        else:
            del seq_id_list[-1]
    fasta_f.close()
    f.close()
    if remove_rna_f:
        os.remove(f_path)
    return fasta_path, seq_id_list


def get_btax_dict(genomes_list,
                  btax_level,
                  btc_profiles,
                  db_dir,
                  num_threads=None,
                  build_tree=False,
                  config_path=None,
                  **kwargs):

    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)
    if not num_threads:
        num_threads = conf_constants.num_threads
    else:
        conf_constants.num_threads = num_threads

    btax_dict = defaultdict(BtaxInfo)
    btc_fasta_dict = defaultdict(dict)
    seq_ids_to_orgs = dict()
    for genome_dict in genomes_list:
        if not genome_dict:
            continue
        genome_info = GenomeInfo.load_from_dict(genome_dict)
        if not genome_info.btc_seqs_id:
            continue
        btax_name = None
        try:
            btax_name = genome_info.taxonomy[-btax_level]
        except IndexError:
            btax_name = genome_info.taxonomy[0]
        btax_dict[btax_name].genomes.append(genome_info.get_json())
        if btax_dict[btax_name].name is None:
            btax_dict[btax_name].name = btax_name
        btc_seqs_fasta_dict = load_fasta_to_dict(genome_info.btc_seqs_fasta)
        for btc_seq_id in genome_info.btc_seqs_id:
            seq_ids_to_orgs[btc_seq_id] = genome_info.org_name
            btc_fasta_dict[genome_info.btc_seqs_id[btc_seq_id]][btc_seq_id] = btc_seqs_fasta_dict[btc_seq_id]

    btc_profile_types = dict()
    for btc_profile_dict in btc_profiles:
        btc_profile_info = SeqProfileInfo.load_from_dict(btc_profile_dict)
        btc_profile_types[btc_profile_info.name] = btc_profile_info.seq_type
    btc_dist_dict = defaultdict(pandas.DataFrame)
    short_to_full_seq_names = dict()
    for btc_profile_name in btc_fasta_dict:
        btc_mult_aln = construct_mult_aln(seq_dict=btc_fasta_dict[btc_profile_name],
                                          aln_type=btc_profile_types[btc_profile_name],
                                          aln_name=btc_profile_name+"_aln",
                                          tmp_dir=kwargs.get("aln_tmp_dir", "mult_aln_tmp"),
                                          method=conf_constants_db.btc_profile_aln_method,
                                          num_threads=num_threads)
        btc_mult_aln.short_to_full_seq_names = short_to_full_seq_names.copy()
        btc_mult_aln.remove_paralogs(seq_ids_to_orgs=seq_ids_to_orgs, inplace=True)
        btc_mult_aln.improve_aln(inplace=True)
        btc_dist_dict[btc_profile_name] = btc_mult_aln.get_distance_matrix()  # TODO: implement specific positions method
        short_to_full_seq_names.update(btc_mult_aln.short_to_full_seq_names)
        if kwargs.get("save_alignments", False):
            btc_mult_aln.dump_alignment(aln_fasta_path=os.path.join(db_dir, btc_mult_aln.aln_name+".fasta"))

    global_dist_matr = get_global_dist(btc_dist_dict, btc_profiles, seq_ids_to_orgs)
    global_dist_matr_path = os.path.join(db_dir, BACTERIA_GLOBAL_DIST_MATRIX)
    short_to_full_seq_names_path = os.path.join(db_dir, BACTERIA_SHORT_TO_FULL_ORG_NAMES)
    short_to_full_seq_names = global_dist_matr.dump(matrix_path=global_dist_matr_path, matr_format="phylip")
    with open(short_to_full_seq_names_path, "w") as short_to_full_org_names_f:
        json.dump(short_to_full_seq_names, short_to_full_org_names_f, indent=2)

    # Steps follows can be parallel
    # for btax_name in btax_dict:
    #     btax_dict[btax_name] = filter_btax(btax_dict[btax_name], ...)

    # while btax_to_merge and n_it < conf_constants_db.:
    #    for btax_name in btax_to_merge:
    #
    #        btax_dict[]
    #    n_it += 1

    full_to_short_seq_names = {v: k for k, v in short_to_full_seq_names.items()}
    for btax_name in btax_dict:
        btax_orgs = set(GenomeInfo.load_from_dict(genome).org_name for genome in btax_dict[btax_name].genomes)
        if build_tree:
            btax_dict[btax_name].mean_d = global_dist_matr[btax_orgs].mean_dist
            btax_dict[btax_name].median_d = global_dist_matr[btax_orgs].median_dist
            if len(btax_orgs) > 2:
                btax_dict[btax_name].ref_tree_newick = build_tree_by_dist(global_dist_matr[btax_orgs],
                                                                          tree_name=btax_name+"_tree").newick
            # btax_dict[btax_name].repr_profiles =
        btax_dict[btax_name].ref_tree_full_names = \
            {full_to_short_seq_names[btax_org]: btax_org for btax_org in btax_orgs}
        btax_dict[btax_name] = btax_dict[btax_name].get_json()
    return btax_dict


def get_global_dist(btc_dist_dict, btc_profiles, seq_ids_to_orgs):
    seqs_order = {org_name: i for i, org_name in enumerate(set(seq_ids_to_orgs.values()))}
    nseqs = len(seqs_order)
    global_dist_matrix = DistanceMatrix(seqs_order=seqs_order,
                                        matr=np.zeros((nseqs, nseqs)),
                                        aln_type="btc_global")
    print(global_dist_matrix.matr.shape)###
    sumw = 0.0
    for btc_profile in btc_profiles:
        btc_profile_info = SeqProfileInfo.load_from_dict(btc_profile)
        btc_profile_matr = btc_dist_dict[btc_profile_info.name]
        btc_profile_orgs = {seq_ids_to_orgs[btc_seq_name]: btc_seq_name for btc_seq_name in btc_profile_matr.seq_names}
        dist_0 = btc_profile_matr.mean_dist
        eagle_logger.info("%s distance mean %s" % (btc_profile_info.name, dist_0))
        for i, seq_name in enumerate(global_dist_matrix.seq_names):
            if seq_name in btc_profile_orgs:
                btc_seq_dist = btc_profile_matr[btc_profile_orgs[seq_name]]
                seq_dists = pandas.Series(([0.0] * i + [btc_seq_dist.get(btc_profile_orgs[seq_name_], dist_0)
                                           for seq_name_ in global_dist_matrix.seq_names[i:]]),
                                          index=global_dist_matrix.seq_names)
            else:
                seq_dists = pandas.Series([0.0] * i + [dist_0] * (nseqs-i), index=global_dist_matrix.seq_names)
                seq_dists[seq_name] = 0.0
            global_dist_matrix[seq_name] += seq_dists * float(btc_profile_info.weight)
        sumw += float(btc_profile_info.weight)
    global_dist_matrix.matr = global_dist_matrix.matr / sumw
    return global_dist_matrix


def get_btax_blastdb(btax_dict, db_dir, btr_profiles=None, num_threads=None, config_path=None):
    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)
    if not num_threads:
        num_threads = conf_constants.num_threads
    else:
        conf_constants.num_threads = num_threads

    # Can be parallel
    for btax_name in btax_dict:
        btax_info = BtaxInfo.load_from_dict(btax_dict[btax_name])
        btax_info.btax_fna, btax_info.fna_id, downloaded_fna = get_btax_fna(btax_genomes=btax_info.genomes,
                                                                            btax_name=btax_info.name,
                                                                            db_dir=db_dir)
        for i, btax_genome in enumerate(btax_info.genomes):
            genome_info = GenomeInfo.load_from_dict(btax_genome)
            if genome_info.genome_id in downloaded_fna:
                genome_info.fna_path = downloaded_fna[genome_info.genome_id]
                btax_info.genomes[i] = genome_info.get_json()
        btax_info.blastdb = create_btax_blastdb(btax_fna_path=btax_info.btax_fna,
                                                btax_name=btax_info.name,
                                                db_dir=db_dir,
                                                blast_inst_dir=conf_constants.blast_inst_dir,
                                                logger=eagle_logger)
        if btr_profiles is not None:
            # create repr profile
            pass
        btax_dict[btax_name] = btax_info.get_json()
    return btax_dict


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
        if not bacterium:
            continue
        if not os.path.exists(bacterium["16S_rRNA_file"]):
            continue
        bacterium_data = {"download_prefix": bacterium["download_prefix"],
                          "16S_rRNA_file": bacterium["16S_rRNA_file"],
                          "fna_file": None,
                          "source_db": bacterium["source_db"],
                          "repr": bacterium['repr']}
        if only_repr and not bacterium['repr']: continue
        if families_dict.get(bacterium['family'], None):
            if families_dict[bacterium['family']].get(bacterium['genus'], None):
                if families_dict[bacterium['family']][bacterium['genus']].get(bacterium['species'], None):
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
            families_dict[bacterium['family']] = \
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
    prepare_families(families_dict, db_dir, bact_fam_f_path, num_threads=num_threads)

    return json.load(open(bact_fam_f_path))


def prepare_families(families_dict, db_dir, bact_fam_f_path, num_threads=4):
    bact_fam_f = io.open(bact_fam_f_path, 'w', newline="\n")
    bact_fam_f.write(u"{\n")
    bact_fam_f.close()

    params_list = list()
    for family in families_dict.keys():
        params_list.append({'function': prepare_family,
                            'family_name': family,
                            'family_data': families_dict[family],
                            'bact_fam_f_path': bact_fam_f_path,
                            'db_dir': db_dir,
                            'logger_name': eagle_logger.name,
                            'try_err_message': "%s is not prepared: " % family})

    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
    pool.close()
    pool.join()

    bact_fam_f = io.open(bact_fam_f_path, 'a', newline="\n")
    bact_fam_f.write(u'  "db_dir": "%s",\n  "db_repr_profiles": "%s"\n}' %
                     (db_dir, os.path.join(db_dir, PROFILES_DB_NAME)))
    bact_fam_f.close()


def prepare_family(family_name, family_data, bact_fam_f_path, db_dir, **kwargs):
    # TODO: refactor it
    special_keys = ("16S_rRNA_tree", "16S_rRNA_tsv", "16S_rRNA_fasta", "blastdb", "repr_profile")
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
    red, reduced_orgs = reduce_seq_names(fasta_dict=dict(map(lambda x: (x, True), set(ids_to_org_dict.values()))),
                                         num_letters=7,
                                         num_words=2)
    rev_reduced_orgs = dict(map(lambda x: (x[1], x[0]), reduced_orgs.items()))
    comp_seq_id_dict = defaultdict(int)
    short_ids_dict = dict()
    for seq_id in rRNA_seqs_dict.keys():
        short_seq_id = None
        short_seq_id = rev_reduced_orgs[ids_to_org_dict[seq_id]]+"x"+\
                       str(get_un_fix(un_num=comp_seq_id_dict[rev_reduced_orgs[ids_to_org_dict[seq_id]]], fix_len=2))
        comp_seq_id_dict[rev_reduced_orgs[ids_to_org_dict[seq_id]]] += 1
        short_ids_dict[short_seq_id] = rev_reduced_orgs[ids_to_org_dict[seq_id]]+"x"+seq_id
        ids_to_org_dict[short_ids_dict[short_seq_id]] = {"organism_name": ids_to_org_dict.pop(seq_id)}
        rRNA_seqs_dict[short_seq_id] = rRNA_seqs_dict.pop(seq_id)
    eagle_logger.info("%s rRNA loaded" % family_name)
    # TODO: follows
    ### This section will be upgraded with my own alignment method but now MUSCLE and hmmer 3 are used
    tmp_fam_dir = os.path.join(db_dir, family_name+"_tmp")
    rRNA_aln = construct_mult_aln(seq_dict=rRNA_seqs_dict,
                                  method="MUSCLE",
                                  aln_type="nucl",
                                  aln_name=family_name+"_rRNA_aln",
                                  tmp_dir=tmp_fam_dir,
                                  muscle_exec_path=conf_constants.muscle_exec_path,
                                  emboss_inst_dir=conf_constants.emboss_inst_dir,
                                  hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                                  logger=eagle_logger)
    eagle_logger.info("%s rRNA alignment constructed" % family_name)
    rRNA_aln.short_to_full_seq_names = short_ids_dict
    # May be errors: not tested
    rRNA_aln.remove_paralogs(ids_to_org_dict, method="min_dist", inplace=True)  # If I use my own alignment method: method="spec_pos"
    rRNA_tree = build_tree_by_dist(rRNA_aln.get_distance_matrix(),
                                   full_seq_names=rRNA_aln.full_to_short_seq_names,
                                   tmp_dir=tmp_fam_dir,
                                   logger=eagle_logger)
    # TODO: write it for not only repr bacteria usage
    # fam_tax = {family_name: get_tree_from_dict(family_data, stop_level=3, special_keys=special_keys)}
    # rRNA_tree, removed_seqs = rRNA_tree.according_to_taxonomy(taxonomy=fam_tax)
    # rRNA_aln.remove_seqs(seqs_list=removed_seqs)
    ###
    family_data["16S_rRNA_tree"] = {"newick": rRNA_tree.newick,
                                    "full_seq_names": rRNA_tree.full_seq_names}
    family_data["16S_rRNA_tsv"] = os.path.join(db_dir, family_name+"_16S_rRNA.tsv")
    family_data["16S_rRNA_fasta"] = os.path.join(db_dir, family_name+"_16S_rRNA.fasta")
    rRNA_aln.get_blocks_tsv(tsv_path=family_data["16S_rRNA_tsv"],
                            fasta_path=family_data["16S_rRNA_fasta"],
                            meta_dict=ids_to_org_dict)
    remained_orgs = list(map(lambda seq_id: ids_to_org_dict[seq_id]["organism_name"], rRNA_aln.seqs()))
    family_data = clean_btax_data(family_data, remained_orgs, stop_level=3, special_keys=special_keys)
    family_data = download_btax_files(key_prefix_pairs={"fna_file": "_genomic.fna.gz"},
                                      btax_data=family_data,
                                      download_dir=db_dir,
                                      logger=eagle_logger)
    family_data["fam_fna"], family_data["chr_id"] = get_btax_fna(fna_key="fna_file",
                                                                 btax_genomes=family_data,
                                                                 btax_name=family_name,
                                                                 db_dir=db_dir)
    family_data["blastdb"] = create_btax_blastdb(btax_fna_path=family_data["fam_fna"],
                                                 btax_name=family_name,
                                                 db_dir=db_dir,
                                                 blast_inst_dir=conf_constants.blast_inst_dir,
                                                 logger=eagle_logger)
    # repr_alns = <function that builds alignments for set of representative genes (returns dict = {aln_name: MultAln object})>
    family_data["repr_profile"] = generate_btax_profile(source={"16S_rRNA": rRNA_aln},
                                                        db_dir=db_dir,
                                                        btax_name=family_name,
                                                        method="hmmer")  # TODO: the source should be repr_alns
    # family_data["codon_usage"] = get_btax_cu(family_data)
    bact_fam_json_f = open(bact_fam_f_path, 'a')
    bact_fam_json_f.write('  "'+family_name+'": '+json.dumps(family_data)+",\n")
    bact_fam_json_f.close()
    eagle_logger.info("%s prepared" % family_name)


def create_bactdb(input_table_refseq=None,
                  input_table_genbank=None,
                  input_table_custom=None,
                  btax_level=int(),
                  btax_class_profile=None,
                  btax_rep_profile=None,
                  db_dir=DEFAULT_BACTDB_DIR,
                  num_threads=None,
                  prepared_genomes=PREPARED_BACTERIA_F_NAME,
                  prepared_genomes_info=None,
                  config_path=None,
                  **kwargs):

    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)
    if not btax_level:
        btax_level = conf_constants_db.btax_level
    else:
        conf_constants_db.btax_level = btax_level
    if not db_dir:
        db_dir = DEFAULT_BACTDB_DIR
    if num_threads:
        int_num_threads = int(num_threads)
        num_threads = None
        num_threads = int_num_threads
        conf_constants.num_threads = num_threads
    else:
        num_threads = conf_constants.num_threads
    if not prepared_genomes:
        prepared_genomes = PREPARED_BACTERIA_F_NAME

    if btax_class_profile is not None:
        # TODO: implement loading btc_profiles from custom profiles
        eagle_logger.warning("custom btax classification profiles are not implemented currently - default will be used")
    # else:
    btc_profiles = [SeqProfileInfo(name="16S_rRNA", seq_type="nucl").get_json()]  # TODO: include it to 'else' bock

    if btax_rep_profile is not None:
        # TODO: implement loading btr_profiles from custom profiles
        eagle_logger.warning("custom btax representative profiles are not implemented currently - default will be used")
    # else:
    btr_profiles = None  # TODO: include it to 'else' bock

    # TODO: this code should not get the btax classification sequence (16S rRNA)
    if input_table_custom is None and input_table_refseq is None and input_table_genbank is None:
        input_table_refseq = DEFAULT_REFSEQ_BACTERIA_TABLE
        input_table_genbank = DEFAULT_GENBANK_BACTERIA_TABLE
    bacteria_list = list()
    if input_table_refseq is not None or input_table_genbank is not None:
        bacteria_list = get_bacteria_from_ncbi(refseq_bacteria_table=input_table_refseq,
                                               genbank_bacteria_table=input_table_genbank,
                                               bactdb_dir=db_dir,
                                               num_threads=num_threads,
                                               prepared_bacteria_f_path=prepared_genomes)
    if input_table_custom is not None:
        eagle_logger.warning("custom genomes input is not implemented yet")
        # TODO: implement custom genomes input
        # bacteria_list.extend()
    if prepared_genomes_info:
        with open(prepared_genomes_info) as prep_genomes_info_f:
            bacteria_list = join_genomes_lists(genomes_list_1=bacteria_list,
                                               genomes_list_2=json.load(prep_genomes_info_f))

    # TODO: implement code to obtain btax classification sequence from fna with hmm profile
    # profiles input should be a list of SeqProfilesInfo objects
    # result - btc_seqs_path field of GenomeInfo objects in bacteria_list filled
    # currently it is filled during get_bacteria_from_ncbi run - not good

    btax_dict = get_btax_dict(genomes_list=bacteria_list,
                              btax_level=btax_level,
                              btc_profiles=btc_profiles,
                              db_dir=db_dir,
                              num_threads=num_threads,
                              build_tree=not bool(btr_profiles))

    btax_dict = get_btax_blastdb(btax_dict,
                                 db_dir=db_dir,
                                 btr_profiles=btr_profiles,
                                 num_threads=num_threads)

    repr_profiles_path = create_profiles_db(btax_dict,
                                            db_dir=db_dir,
                                            profiles_db_name=PROFILES_DB_NAME,
                                            method="hmmer",
                                            hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                                            config_path=config_path,
                                            logger=eagle_logger)
    with open(os.path.join(db_dir, BTAX_JSON_NAME), "w") as btax_json_f:
        json.dump(btax_dict, btax_json_f, indent=2)  # maybe btax_dict will be dumped in get_btax_dict
    db_info = DBInfo(all_genomes=os.path.join(db_dir, BACTERIA_LIST_F_NAME),
                     btax_json=os.path.join(db_dir, BTAX_JSON_NAME),
                     repr_profiles=repr_profiles_path,
                     global_dist_matrix=os.path.join(db_dir, BACTERIA_GLOBAL_DIST_MATRIX),
                     all_org_full_names=os.path.join(db_dir, BACTERIA_SHORT_TO_FULL_ORG_NAMES)).get_json()
    with open(os.path.join(db_dir, DB_INFO_NAME), "w") as db_info_f:
        json.dump(db_info, db_info_f, indent=2)
    return db_info
