# WARNING: this module cannot be imported in any eagle.lib module!
import os
import multiprocessing as mp

import numpy as np
import pandas as pd

from eagle.constants import conf_constants
from eagle.lib.general import send_log_message, join_files, worker
from eagle.lib.seqs import SeqsDict, read_blast_out
from eagle.lib.alignment import MultAln, BlastHandler, search_profile
from eagle.lib.phylo import PhyloTree, compare_trees, build_tree_by_dist
from eagledb.scheme import GenomeInfo, SeqProfileInfo


def hom_search_blast(seqs_dict, genomes_list, blastdb_path=None, work_dir="./blast_search", low_memory="auto",
                     num_threads=None, blast_inst_dir=None, autoremove=True):
    assert isinstance(seqs_dict, SeqsDict), "ERROR: the value for 'seqs_dict' argument should be an object " \
                                            "of eagle.lib.seqs.SeqsDict class"

    if num_threads is None:
        num_threads = conf_constants.num_threads
    if blast_inst_dir is None:
        blast_inst_dir = conf_constants.blast_inst_dir

    genome_id2org = dict()
    genome_id2fna = dict()
    for genome_dict in genomes_list:
        genome_info = GenomeInfo.load_from_dict(in_dict=genome_dict)
        genome_id2org[genome_info.genome_id] = genome_info.org_name
        genome_id2fna[genome_info.genome_id] = genome_info.fna_path
    blast_handler = BlastHandler(inst_dir=blast_inst_dir)
    if blastdb_path is None:
        blastdb_path = os.path.join(work_dir, "genomes_blastdb", "genomes")
        os.makedirs(os.path.dirname(blastdb_path))
        genomes_fasta = join_files(list(genome_id2fna.values()), out_file_path=blastdb_path+".fasta")
        blast_handler.make_blastdb(in_fasta=genomes_fasta, dbtype="nucl", db_name=blastdb_path)
    if seqs_dict.seqs_type in seqs_dict.prot_types:
        blast_type = "tblastn"
    else:
        blast_type = "blastn"
    blast_out_path = os.path.join(work_dir, "blast.out")
    blast_handler.run_blast_search(blast_type=blast_type, query=seqs_dict, db=blastdb_path, out=blast_out_path,
                                   num_threads=num_threads)
    blast_res_dict = read_blast_out(blast_out_path=blast_out_path)

    for seq_name in seqs_dict:
        pass
    return


def hom_search_profile(alns_or_profiles, genomes_list, work_dir="./profile_search", low_memory="auto",
                       num_threads=None, autoremove=True):
    if num_threads is None:
        num_threads = conf_constants.num_threads

    seq_types = mp.Manager().dict()
    seq_types["prot"] = False
    seq_types["nucl"] = False
    seq_ids_to_orgs = mp.Manager().dict()
    protdb_path = os.path.join(work_dir, "protdb.fasta")
    nucldb_path = os.path.join(work_dir, "nucldb.fasta")

    pool = mp.Pool(num_threads)
    profiles_list = pool.map(worker, map(lambda a_or_p: {"function": _prepare_profile,
                                                         "a_or_p": a_or_p,
                                                         "work_dir": work_dir,
                                                         "seq_types": seq_types,
                                                         "protdb_path": protdb_path,
                                                         "nucldb_path": nucldb_path,
                                                         "seq_ids_to_orgs": seq_ids_to_orgs},
                                         alns_or_profiles))
    pool.close()
    pool.join()

    for genome_dict in genomes_list:
        genome_info = GenomeInfo.load_from_dict(in_dict=genome_dict)
        if seq_types["prot"]:
            pass
        if seq_types["nucl"]:
            pass

    pool = mp.Pool(num_threads)
    hom_alns = pool.map(worker, profiles_list)
    pool.close()
    pool.join()
    return hom_alns


def _prepare_profile(a_or_p, work_dir, seq_types, protdb_path, nucldb_path, seq_ids_to_orgs):
    if isinstance(a_or_p, MultAln):
        profile_path = None
        profile_path = os.path.join(work_dir, a_or_p.aln_name, "profile.hmm")
        a_or_p.get_hmm_profile(profile_path=profile_path, method="hmmer")
        p = SeqProfileInfo(name=a_or_p.aln_name + "_profile", path=profile_path, seq_type=a_or_p.aln_type)
    else:
        p = SeqProfileInfo.load_from_dict(in_dict=a_or_p)
    if p.seq_type in SeqsDict.prot_types:
        seq_types["prot"] = True
        seqdb = protdb_path
    else:
        seq_types["nucl"] = True
        seqdb = nucldb_path
    return {"function": search_profile,
            "profile_dict": p.get_json(),
            "seqdb": seqdb,
            "seq_ids_to_orgs": seq_ids_to_orgs}


def explore_ortho_group(homologs_mult_aln, remove_paralogs=True, ref_tree_newick=None, **kwargs):
    assert isinstance(homologs_mult_aln, MultAln), \
        "ERROR: the value for 'homologs_mult_aln' argument is not an object of eagle.lib.alignment.MultAln class"
    logger = kwargs.get("logger", None)
    work_dir = kwargs.get("work_dir", "./")
    ref_tree_full_names = kwargs.get("ref_tree_full_names", None)
    kaks_for_only_first = kwargs.get("kaks_for_only_first", False)
    orgs_tree = kwargs.get("orgs_tree", False)
    # TODO: store parameters for paralogs removing, paralogs removing, tree construction etc.

    stats_dict = {
        "uniformity_std": -1.0,
        "phylo_diff": -1.0,
        "Ka/Ks": -1.0,
    }

    if remove_paralogs and homologs_mult_aln.seq_ids_to_orgs:
        homologs_mult_aln.remove_paralogs(method="min_dist", inplace=True)
    if len(homologs_mult_aln.seq_names) < 4:
        send_log_message(message="A few (<4) homologs number in the alignment '%s'" % homologs_mult_aln.aln_name,
                         mes_type="w", logger=logger)
        return homologs_mult_aln, stats_dict, None
    homologs_mult_aln.improve_aln(inplace=True)

    # Uniformity
    stats_dict["uniformity_std"] = homologs_mult_aln.estimate_uniformity(
        cons_thr=conf_constants.cons_thr,
        window_l=conf_constants.unif_window_l,
        windows_step=conf_constants.unif_windows_step
    )
    if np.isnan(stats_dict["uniformity_std"]):
        stats_dict["uniformity_std"] = -1.0
    send_log_message(message="got uniformity_std for the alignment '%s'" % homologs_mult_aln.aln_name,
                     mes_type="info", logger=logger)

    # Ka/Ks
    if homologs_mult_aln.aln_type.lower() in MultAln.prot_types and homologs_mult_aln.nucl_seqs_dict:
        aln_kaks = homologs_mult_aln.calculate_KaKs_windows(only_first_seq=kaks_for_only_first)
        if not pd.isna(aln_kaks):
            stats_dict["Ka/Ks"] = aln_kaks
            send_log_message(message="got Ka/Ks for the alignment '%s'" % homologs_mult_aln.aln_name,
                             mes_type="info", logger=logger)

    # Phylo
    if orgs_tree:
        org_names_dict = dict((short_id, homologs_mult_aln.seq_ids_to_orgs[full_id])
                              for short_id, full_id in homologs_mult_aln.short_to_full_seq_names.items())
    phylo_tmp_dir = os.path.join(work_dir, homologs_mult_aln.aln_name.replace("|:", "_")+"_phylo_tmp")
    homologs_tree = build_tree_by_dist(
        dist_matrix=homologs_mult_aln.get_distance_matrix(),
        method="FastME",
        full_seq_names=org_names_dict if orgs_tree else homologs_mult_aln.short_to_full_seq_names,
        tree_name=homologs_mult_aln.aln_name.replace("|:", "_")+"_tree",
        tmp_dir=phylo_tmp_dir,
        logger=logger
    )
    homologs_tree.set_full_names(inplace=True)
    if ref_tree_newick is not None:
        ref_tree = PhyloTree.load_tree_from_str(
            tree_str=ref_tree_newick,
            full_seq_names=ref_tree_full_names,
            tree_name="ref_tree",
            tmp_dir=phylo_tmp_dir,
            logger=logger
        )
        ref_tree.set_full_names(inplace=True)
        phylo_diff = compare_trees(phylo_tree1=homologs_tree,
                                   phylo_tree2=ref_tree,
                                   method="Robinson-Foulds")
        if not pd.isna(phylo_diff):
            stats_dict["phylo_diff"] = phylo_diff
        send_log_message(message="got phylo_diff for the alignment '%s'" % homologs_mult_aln.aln_name,
                         mes_type="info", logger=logger)

    return homologs_mult_aln, stats_dict, homologs_tree
