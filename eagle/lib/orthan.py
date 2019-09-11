# WARNING: this module cannot be imported in any eagle.lib module!
import os

import numpy as np
import pandas as pd

from eagle.constants import conf_constants
from eagle.lib.general import send_log_message
from eagle.lib.alignment import MultAln
from eagle.lib.phylo import PhyloTree, compare_trees, build_tree_by_dist


def explore_ortho_group(homologs_mult_aln, ref_tree_newick=None, **kwargs):
    assert isinstance(homologs_mult_aln, MultAln), \
        "ERROR: the value for 'homologs_mult_aln' argument is not an object of eagle.lib.alignment.MultAln class"
    logger = kwargs.get("logger", None)
    work_dir = kwargs.get("work_dir", "./")
    ref_tree_full_names = kwargs.get("ref_tree_full_names", None)
    kaks_for_only_first = kwargs.get("kaks_for_only_first", False)

    stats_dict = {
        "uniformity_std": -1.0,
        "phylo_diff": -1.0,
        "Ka/Ks": -1.0,
    }

    homologs_mult_aln.remove_paralogs(seq_ids_to_orgs=seq_ids_to_orgs, method="min_dist", inplace=True)
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
    if homologs_mult_aln.aln_type.lower() in MultAln.prot_types:
        aln_kaks = homologs_mult_aln.calculate_KaKs_windows(only_first_seq=kaks_for_only_first)
        if not pd.isna(aln_kaks):
            stats_dict["Ka/Ks"] = aln_kaks
            send_log_message(message="got Ka/Ks for the alignment '%s'" % homologs_mult_aln.aln_name,
                             mes_type="info", logger=logger)

    # Phylo
    phylo_tmp_dir = os.path.join(work_dir, homologs_mult_aln.aln_name.replace("|:", "_")+"_phylo_tmp")
    homologs_tree = build_tree_by_dist(
        dist_matrix=homologs_mult_aln.get_distance_matrix(),
        method="FastME",
        full_seq_names=dict((short_id, seq_ids_to_orgs[full_id])
                            for short_id, full_id in homologs_mult_aln.short_to_full_seq_names.items()),
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
