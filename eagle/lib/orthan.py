# WARNING: this module cannot be imported in any eagle.lib module!
import os
import multiprocessing as mp

import numpy as np
import pandas as pd

from eagle.constants import conf_constants
from eagle.lib.general import join_files
from eagle.lib.workers import process_worker
from eagle.lib.logging import send_log_message
from eagle.lib.seqs import SeqsDict, read_blast_out, get_orfs
from eagle.lib.alignment import MultAln
from eagle.lib.phylo import PhyloTree, compare_trees
from eagledb.scheme import GenomeInfo, SeqProfileInfo


def hom_search_blast(seqs_dict, genomes_list, blastdb_path=None, work_dir="./blast_search", low_memory="auto",
                     num_threads=None, blast_inst_dir=None, autoremove=True):
    assert isinstance(seqs_dict, SeqsDict), "ERROR: the value for 'seqs_dict' argument should be an object " \
                                            "of eagle.lib.seqs.SeqsDict class"

    if num_threads is None:
        num_threads = conf_constants.num_threads
    if blast_inst_dir is None:
        blast_inst_dir = conf_constants.blast_inst_dir

    # TODO: this part can be parallel
    fna_id2org = dict()
    fna_paths = set()
    for genome_dict in genomes_list:
        genome_info = GenomeInfo.load_from_dict(in_dict=genome_dict)
        fna_id2org.update(_link_fna_id2org(genome_info))
        fna_paths.add(genome_info.fna_path)

    blast_handler = BlastHandler(inst_dir=blast_inst_dir)
    if blastdb_path is None:
        blastdb_path = os.path.join(work_dir, "genomes_blastdb", "genomes")
        os.makedirs(os.path.dirname(blastdb_path))
        genomes_fasta = join_files(fna_paths, out_file_path=blastdb_path+".fasta")
        blast_handler.make_blastdb(in_fasta=genomes_fasta, dbtype="nucl", db_name=blastdb_path)
    if seqs_dict.seqs_type in seqs_dict.prot_types:
        blast_type = "tblastn"
    else:
        blast_type = "blastn"
    blast_out_path = os.path.join(work_dir, "blast.out")
    blast_handler.run_blast_search(blast_type=blast_type, query=seqs_dict, db=blastdb_path, out=blast_out_path,
                                   num_threads=num_threads)
    blast_res_dict = read_blast_out(blast_out_path=blast_out_path)  # TODO: consider to reduce memory of 'blast_res_dict'

    for seq_name in seqs_dict:
        pass
    return


def hom_search_profile(alns_or_profiles, genomes_list, work_dir="./profile_search", low_memory="auto",
                       num_threads=None, autoremove=True):
    if num_threads is None:
        num_threads = conf_constants.num_threads
    try:
        cpu = max(num_threads//len(genomes_list), 1)
    except TypeError:
        cpu = 1

    seq_types = mp.Manager().dict()
    seq_types["prot"] = False
    seq_types["nucl"] = False
    fna_id2org = mp.Manager().dict()
    protdb_path = os.path.join(work_dir, "protdb.fasta")
    nucldb_path = os.path.join(work_dir, "nucldb.fasta")

    pool = mp.Pool(num_threads)
    profiles_list = pool.map(process_worker, map(lambda a_or_p: {"function": _prepare_profile,
                                                         "a_or_p": a_or_p,
                                                         "work_dir": work_dir,
                                                         "seq_types": seq_types,
                                                         "protdb_path": protdb_path,
                                                         "nucldb_path": nucldb_path,
                                                         "fna_id2org": fna_id2org,
                                                         "cpu": cpu},
                                                 alns_or_profiles))
    pool.close()
    pool.join()

    pool = mp.Pool(num_threads)
    db_part_paths = pool.map(process_worker, map(lambda genome_dict: {
                                             "function": _prepare_genome_db,
                                             "genome_dict": genome_dict,
                                             "fna_id2org": fna_id2org,
                                             "work_dir": work_dir,
                                             "seq_prot": seq_types["prot"],
                                             "seq_nucl": seq_types["nucl"]
                                         },
                                                 genomes_list))
    pool.close()
    pool.join()

    orfs_fasta_paths = set()
    nucl_fasta_paths = set()
    for db_part_path in db_part_paths:
        if db_part_path["prot"] is not None:
            orfs_fasta_paths.add(db_part_path["prot"])
        if db_part_path["nucl"] is not None:
            nucl_fasta_paths.add(db_part_path["nucl"])
    join_files(in_files_list=orfs_fasta_paths, out_file_path=protdb_path, remove_infiles=True)
    join_files(in_files_list=nucl_fasta_paths, out_file_path=nucldb_path, remove_infiles=True)

    pool = mp.Pool(num_threads)
    hom_alns = pool.map(process_worker, profiles_list)
    pool.close()
    pool.join()
    return hom_alns


def _link_fna_id2org(genome_info):
    assert isinstance(genome_info, GenomeInfo), "ERROR: the value for 'genome_info' argument should be \
                                                 an object of GenomeInfo class"
    if genome_info.fna_id_list:
        return {fna_id: genome_info.org_name for fna_id in genome_info.fna_id_list}
    else:
        return {
            fna_id: genome_info.org_name for fna_id in SeqsDict.load_from_file(genome_info.fna_path).keys()
        }    


def _prepare_profile(a_or_p, work_dir, seq_types, protdb_path, nucldb_path, fna_id2org, cpu):
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
            "seqdb_id2org": fna_id2org,
            "shred_seqdb": True if p.seq_type.lower() in SeqsDict.nucl_types else False,
            "read_hsr_shred": True if p.seq_type.lower() in SeqsDict.prot_types else False,
            "hit_coords_transform": lambda c: c*3 if p.seq_type.lower() in SeqsDict.prot_types else None,
            "cpu": cpu}


def _prepare_genome_db(genome_dict, fna_id2org, work_dir, seq_prot, seq_nucl, **kwargs):
    genome_db_path = {"prot": None, "nucl": None}
    genome_info = GenomeInfo.load_from_dict(in_dict=genome_dict)
    if seq_prot:
        genome_db_path["prot"] = os.path.join(work_dir, os.path.splitext(os.path.basename(genome_info.fna_path))[0]) + \
                                 "_orfs.fasta"
        orf_ids = get_orfs(in_fasta_path=genome_info.fna_path, 
                           out_fasta_path=genome_db_path["prot"], 
                           minsize=90).keys()
        if genome_info.fna_id_list:
            fna_id2org.update({orf_id: genome_info.org_name for orf_id in 
                               filter(lambda orfid: orfid in genome_info.fna_id_list, orf_ids)})
        else:      
            fna_id2org.update({orf_id: genome_info.org_name for orf_id in orf_ids})
    if seq_nucl:
        genome_db_path["nucl"] = genome_info.fna_path
        fna_id2org.update(_link_fna_id2org(genome_info))
    return genome_db_path


def explore_ortho_group(homologs_mult_aln, remove_paralogs=True, ref_tree_newick=None, phylo_params=None, **kwargs):
    assert isinstance(homologs_mult_aln, MultAln), \
        "ERROR: the value for 'homologs_mult_aln' argument is not an object of eagle.lib.alignment.MultAln class"
    work_dir = kwargs.get("work_dir", "./")
    ref_tree_full_names = kwargs.get("ref_tree_full_names", None)
    kaks_for_only_first = kwargs.get("kaks_for_only_first", False)
    orgs_tree = kwargs.get("orgs_tree", False)
    get_conflicting_taxons = kwargs.get("get_conflicting_taxons", True)
    # TODO: store parameters for paralogs removing, paralogs removing, tree construction etc.

    stats_dict = {
        "uniformity_std90": -1.0,
        "uniformity_std95": -1.0,
        "uniformity_std98": -1.0,
        "uniformity_std100": -1.0,
        "phylo_diff": -1.0,
        "conflicting_taxons": list(),
        "Ka/Ks": list(),
        "Ks": list(),
    }

    if remove_paralogs and homologs_mult_aln.seq_ids_to_orgs:
        homologs_mult_aln.remove_paralogs(method="min_dist", inplace=True)
    if len(homologs_mult_aln.seq_names) < 4:
        send_log_message(message="A few (<4) homologs number in the alignment '%s'" % homologs_mult_aln.aln_name,
                         mes_type="w", logger=homologs_mult_aln.logger)
        return homologs_mult_aln, stats_dict, None
    homologs_mult_aln.improve_aln(inplace=True)

    # Uniformity
    std90 = homologs_mult_aln.estimate_uniformity(
        cons_thr=90,
        window_l=conf_constants.unif_window_l,
        windows_step=conf_constants.unif_windows_step
    )
    std95 = homologs_mult_aln.estimate_uniformity(
        cons_thr=95,
        window_l=conf_constants.unif_window_l,
        windows_step=conf_constants.unif_windows_step
    )
    std98 = homologs_mult_aln.estimate_uniformity(
        cons_thr=98,
        window_l=conf_constants.unif_window_l,
        windows_step=conf_constants.unif_windows_step
    )
    std100 = homologs_mult_aln.estimate_uniformity(
        cons_thr=100,
        window_l=conf_constants.unif_window_l,
        windows_step=conf_constants.unif_windows_step
    )
    if not np.isnan(std90):
        stats_dict["uniformity_std90"] = std90
    if not np.isnan(std95):
        stats_dict["uniformity_std95"] = std95
    if not np.isnan(std98):
        stats_dict["uniformity_std98"] = std98
    if not np.isnan(std100):
        stats_dict["uniformity_std100"] = std100
    send_log_message(message="got uniformity_std for the alignment '%s'" % homologs_mult_aln.aln_name,
                     mes_type="info", logger=homologs_mult_aln.logger)

    # Ka/Ks
    if homologs_mult_aln.aln_type.lower() in MultAln.prot_types and homologs_mult_aln.nucl_seqs_dict:
        aln_kaks = homologs_mult_aln.calculate_KaKs_windows(only_first_seq=kaks_for_only_first)
        stats_dict["Ka/Ks"], stats_dict["Ks"] = aln_kaks["Ka/Ks"], aln_kaks["Ks"]
        send_log_message(message="got Ka/Ks for the alignment '%s'" % homologs_mult_aln.aln_name,
                         mes_type="info", logger=homologs_mult_aln.logger)

    # Phylo
    if phylo_params is None:
        phylo_params = dict()
    build_tree_method = phylo_params.get("method", "FastME")
    options = phylo_params.get("options", {})
    fastme_exec_path = phylo_params.get("fastme_exec_path", None)

    phylo_tmp_dir = os.path.join(work_dir, homologs_mult_aln.aln_name.replace("|:", "_")+"_phylo_tmp")
    homologs_tree = homologs_mult_aln.build_tree(tree_name=homologs_mult_aln.aln_name.replace("|:", "_")+"_tree",
                                                 tmp_dir=phylo_tmp_dir,
                                                 method=build_tree_method,
                                                 options=options,
                                                 fastme_exec_path=fastme_exec_path)
    if orgs_tree:
        homologs_tree.full_seq_names = dict((short_id, homologs_mult_aln.seq_ids_to_orgs[full_id])
                                            for short_id, full_id in homologs_mult_aln.short_to_full_seq_names.items())
    homologs_tree.set_full_names(inplace=True)
    if ref_tree_newick is not None:
        ref_tree = PhyloTree.load_tree_from_str(
            tree_str=ref_tree_newick,
            full_seq_names=ref_tree_full_names,
            tree_name="ref_tree",
            tmp_dir=phylo_tmp_dir,
            logger=homologs_mult_aln.logger
        )
        ref_tree.set_full_names(inplace=True)
        phylo_diff = compare_trees(phylo_tree1=homologs_tree,
                                   phylo_tree2=ref_tree,
                                   method="Robinson-Foulds")
        if not pd.isna(phylo_diff):
            stats_dict["phylo_diff"] = phylo_diff
        send_log_message(message="got phylo_diff for the alignment '%s'" % homologs_mult_aln.aln_name,
                         mes_type="info", logger=homologs_mult_aln.logger)

    return homologs_mult_aln, stats_dict, homologs_tree
