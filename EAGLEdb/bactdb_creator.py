import gzip
import io
import json
import multiprocessing as mp
import os
from collections import defaultdict

import pandas

from EAGLE.constants import EAGLE_logger, conf_constants
from EAGLE.lib.alignment import construct_mult_aln
from EAGLE.lib.general import worker, load_fasta_to_dict, reduce_seq_names, get_un_fix, bool_from_str
from EAGLE.lib.phylo import build_tree_by_dist
from EAGLEdb import join_bacteria_lists
from EAGLEdb.constants import BACTERIA_LIST_F_NAME, ANALYZED_BACTERIA_F_NAME, BACT_FAM_F_NAME, conf_constants_db, \
    DEFAULT_REFSEQ_BACTERIA_TABLE, DEFAULT_GENBANK_BACTERIA_TABLE, DEFAULT_BACTDB_DIR
from EAGLEdb.lib.db_creator import download_organism_files, clean_btax_data, download_btax_files, create_btax_blastdb, \
    generate_btax_profile


def get_bacteria_from_ncbi(refseq_bacteria_table=None,
                           genbank_bacteria_table=None,
                           bactdb_dir=DEFAULT_BACTDB_DIR,
                           num_threads=None,
                           first_bact=None,
                           last_bact=None,
                           analyzed_bacteria_f_path=ANALYZED_BACTERIA_F_NAME,
                           remove_bact_list_f=False,
                           config_path=None):

    if not refseq_bacteria_table and not genbank_bacteria_table:
        refseq_bacteria_table = DEFAULT_REFSEQ_BACTERIA_TABLE
        genbank_bacteria_table = DEFAULT_GENBANK_BACTERIA_TABLE
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
        EAGLE_logger.warning("bactdb directory exists")
    analyzed_bacteria = mp.Manager().dict()
    if os.path.exists(analyzed_bacteria_f_path):
        EAGLE_logger.info("loading analyzed bacteria from '%s'" % analyzed_bacteria_f_path)
        analyzed_bacteria_f = open(analyzed_bacteria_f_path)
        analyzed_bacteria.update(json.load(analyzed_bacteria_f))
        analyzed_bacteria_f.close()
        EAGLE_logger.info("analyzed bacteria loaded")
    bacteria_list_f_path = os.path.join(bactdb_dir, BACTERIA_LIST_F_NAME)
    bacteria_list_f = io.open(bacteria_list_f_path, 'w', newline="\n")
    bacteria_list_f.write(u"[\n")
    bacteria_list_f.close()
    refseq_df = pandas.read_csv(refseq_bacteria_table, sep="\t", dtype=str)
    genbank_df = pandas.read_csv(genbank_bacteria_table, sep="\t", dtype=str)
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
                params_list.append({'function': get_bacterium,
                                    'analyzed_bacteria': analyzed_bacteria,
                                    'logger_name': EAGLE_logger.name,
                                    'ncbi_db_link': genbank_df.loc[j]["ncbi_link"],
                                    'bacterium_name': genbank_df.loc[j]["org_name"],
                                    'repr': bool_from_str(genbank_df.loc[j]["repr"]),
                                    'db_dir': bactdb_dir,
                                    'source_db': "genbank",
                                    'try_err_message': "%s is not prepared: " % genbank_df.loc[j]["org_name"]})
                j += 1
            else:
                params_list.append({'function': get_bacterium,
                                    'analyzed_bacteria': analyzed_bacteria,
                                    'logger_name': EAGLE_logger.name,
                                    'ncbi_db_link': refseq_df.loc[i]["ncbi_link"],
                                    'bacterium_name': refseq_df.loc[i]["org_name"],
                                    'repr': bool_from_str(refseq_df.loc[i]["repr"]),
                                    'db_dir': bactdb_dir,
                                    'source_db': "refseq",
                                    'try_err_message': "%s is not prepared: " % refseq_df.loc[i]["org_name"]})
                i += 1
        else:
            if genbank_df.loc[j]["org_name"] < refseq_df.loc[i]["org_name"]:
                params_list.append({'function': get_bacterium,
                                    'analyzed_bacteria': analyzed_bacteria,
                                    'logger_name': EAGLE_logger.name,
                                    'ncbi_db_link': genbank_df.loc[j]["ncbi_link"],
                                    'bacterium_name': genbank_df.loc[j]["org_name"],
                                    'repr': bool_from_str(genbank_df.loc[j]["repr"]),
                                    'db_dir': bactdb_dir,
                                    'source_db': "genbank",
                                    'try_err_message': "%s is not prepared: " % genbank_df.loc[j]["org_name"]})
                j += 1
            else:
                params_list.append({'function': get_bacterium,
                                    'analyzed_bacteria': analyzed_bacteria,
                                    'logger_name': EAGLE_logger.name,
                                    'ncbi_db_link': refseq_df.loc[i]["ncbi_link"],
                                    'bacterium_name': refseq_df.loc[i]["org_name"],
                                    'repr': bool_from_str(refseq_df.loc[i]["repr"]),
                                    'db_dir': bactdb_dir,
                                    'source_db': "refseq",
                                    'try_err_message': "%s is not prepared: " % refseq_df.loc[i]["org_name"]})
                i += 1
            if genbank_df.loc[j]["org_name"] == refseq_df.loc[i-1]["org_name"]:
                j += 1
        n += 1
    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
    pool.close()
    pool.join()
    analyzed_bacteria_f = open(os.path.join(bactdb_dir, ANALYZED_BACTERIA_F_NAME), "w")
    json.dump(dict(analyzed_bacteria), analyzed_bacteria_f)
    analyzed_bacteria_f.close()
    bacteria_list_f = io.open(bacteria_list_f_path, 'a', newline="\n")
    bacteria_list_f.write(u"  {}\n]")
    bacteria_list_f.close()
    return json.load(open(bacteria_list_f_path))


def get_bacterium(ncbi_db_link, bacterium_name, repr, analyzed_bacteria, db_dir, source_db=None, **kwargs):
    assembly_id = ncbi_db_link.split("/")[-1]
    bacterium_info = {"family": None,
                      "genus": None,
                      "species": None,
                      "strain": None,
                      "download_prefix": (ncbi_db_link+"/"+assembly_id).replace("https", "ftp"),
                      "16S_rRNA_file": None,
                      "source_db": source_db,
                      "repr": repr}
    EAGLE_logger.info("%s getting started" % bacterium_name)
    EAGLE_logger.info('bacterium link: %s' % bacterium_info["download_prefix"])
    if analyzed_bacteria.get(bacterium_name, None):
        EAGLE_logger.info("%s is in analyzed bacteria" % bacterium_name)
        return 0
    download_organism_files(bacterium_info["download_prefix"],
                            ["_wgsmaster.gbff.gz", "_rna_from_genomic.fna.gz"],
                            db_dir,
                            logger=EAGLE_logger)
    tax_f_name = assembly_id + "_wgsmaster.gbff.gz"
    if not os.path.exists(os.path.join(db_dir, tax_f_name)):
        tax_f_name = None
        download_organism_files(bacterium_info["download_prefix"], "_genomic.gbff.gz", db_dir, logger=EAGLE_logger)
        tax_f_name = assembly_id + "_genomic.gbff.gz"
    bacterium_info["family"], bacterium_info["genus"], bacterium_info["species"], bacterium_info["strain"] = \
        get_taxonomy(tax_f_name, db_dir)
    EAGLE_logger.info("got %s taxonomy" % bacterium_info["strain"])
    #if not os.path.exists():
    #
    bacterium_info["16S_rRNA_file"] = get_16S_fasta(assembly_id + "_rna_from_genomic.fna.gz",
                                                    db_dir,
                                                    bacterium_info["strain"])
    EAGLE_logger.info("got %s 16S rRNA" % bacterium_info["strain"])
    f = io.open(os.path.join(db_dir, BACTERIA_LIST_F_NAME), 'a', newline="\n")
    f.write(unicode("  "+json.dumps(bacterium_info)+",\n"))
    f.close()
    analyzed_bacteria[bacterium_name] = True


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
                            'logger_name': EAGLE_logger.name,
                            'try_err_message': "%s is not prepared: " % family})

    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
    pool.close()
    pool.join()

    bact_fam_f = io.open(bact_fam_f_path, 'a', newline="\n")
    bact_fam_f.write(u"  {}\n}")
    bact_fam_f.close()


def prepare_family(family_name, family_data, bact_fam_f_path, db_dir):
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
    EAGLE_logger.info("%s rRNA loaded" % family_name)
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
                                  logger=EAGLE_logger)
    EAGLE_logger.info("%s rRNA alignment constructed" % family_name)
    rRNA_aln.short_to_full_seq_names = short_ids_dict
    rRNA_aln.full_to_short_seq_names = None
    rRNA_aln.mult_aln_dict = dict(map(
        lambda x: (rRNA_aln.short_to_full_seq_names[x[0]], rRNA_aln.mult_aln_dict_short_id[x[0]]),
        rRNA_aln.mult_aln_dict_short_id.items()))
    rRNA_aln.remove_paralogs(ids_to_org_dict, method="min_dist", inplace=True)  # If I use my own alignment method: method="spec_pos"
    rRNA_tree = build_tree_by_dist(rRNA_aln.get_distance_matrix(), tmp_dir=tmp_fam_dir)
    rRNA_tree.full_seq_names = dict(map(lambda x: (x, ids_to_org_dict[short_ids_dict[x]]),
                                        rRNA_aln.mult_aln_dict_short_id.keys()))
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
    remained_orgs = map(lambda seq_id: ids_to_org_dict[seq_id]["organism_name"], rRNA_aln.seqs())
    family_data = clean_btax_data(family_data, remained_orgs, stop_level=3, special_keys=special_keys)
    family_data = download_btax_files(key_prefix_pairs={"fna_file": "_genomic.fna.gz"},
                                      btax_data=family_data,
                                      download_dir=db_dir,
                                      logger=EAGLE_logger)
    family_data["blastdb"] = create_btax_blastdb(family_data, family_name, db_dir,
                                                 blast_inst_dir=conf_constants.blast_inst_dir,
                                                 logger=EAGLE_logger)
    # repr_alns = <function that builds alignments for set of representative genes (returns dict = {aln_name: MultAln object})>
    family_data["repr_profile"] = generate_btax_profile(source={"16S_rRNA": rRNA_aln},
                                                        db_dir=db_dir,
                                                        btax_name=family_name,
                                                        method="hmmer")  # TODO: the source should be repr_alns
    # family_data["codon_usage"] = get_btax_cu(family_data)
    bact_fam_json_f = open(bact_fam_f_path, 'a')
    bact_fam_json_f.write('  "'+family_name+'": '+json.dumps(family_data)+",\n")
    bact_fam_json_f.close()
    EAGLE_logger.info("%s prepared" % family_name)


def create_bactdb(input_table_refseq=None,
                  input_table_genbank=None,
                  db_dir=DEFAULT_BACTDB_DIR,
                  num_threads=None,
                  analyzed_organisms=ANALYZED_BACTERIA_F_NAME,
                  analyzed_organisms_info=None,
                  config_path=None,
                  **kwargs):

    if config_path:
        conf_constants.update_by_config(config_path=config_path)
        conf_constants_db.update_by_config(config_path=config_path)
    if not db_dir:
        db_dir = DEFAULT_BACTDB_DIR
    if num_threads:
        int_num_threads = int(num_threads)
        num_threads = None
        num_threads = int_num_threads
        conf_constants.num_threads = num_threads
    else:
        num_threads = conf_constants.num_threads
    if not analyzed_organisms:
        analyzed_organisms = ANALYZED_BACTERIA_F_NAME

    bacteria_list = get_bacteria_from_ncbi(refseq_bacteria_table=input_table_refseq,
                                           genbank_bacteria_table=input_table_genbank,
                                           bactdb_dir=db_dir,
                                           num_threads=num_threads,
                                           analyzed_bacteria_f_path=analyzed_organisms)
    if analyzed_organisms_info:
        bacteria_list = join_bacteria_lists(bacteria_list_1=bacteria_list,
                                            bacteria_list_2=json.load(open(analyzed_organisms_info)))
    families_dict = get_families_dict(bacteria_list=bacteria_list,
                                      db_dir=db_dir,
                                      num_threads=num_threads,
                                      only_repr=True)
