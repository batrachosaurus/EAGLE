import os
import shutil
import subprocess

import pandas

from EAGLE.lib.general import ConfBase, dump_phylip_dist_matrix, load_newick


class PhyloTree(ConfBase):

    def __init__(self, newick, full_seq_names=None, config_path=None, logger=None):
        self.newick = newick
        self.full_seq_names = full_seq_names
        self.logger = logger

        super(PhyloTree, self).__init__(config_path=config_path)

    def according_to_taxonomy(self, taxonomy):
        # NOT inplace method!
        pass


def build_tree_by_dist(dist_matrix=None,
                       dist_matrix_f=None,
                       full_seq_names=None,
                       tmp_dir="tmp",
                       method="FastME",
                       fastme_exec_path="fastme",
                       config_path=None,
                       logger=None):

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if type(dist_matrix) is not pandas.DataFrame and not dist_matrix_f:
        if logger:
            logger.warning("No distance matrix input")
        else:
            print "No distance matrix input"
        return 1
    elif dist_matrix.empty and dist_matrix_f:
        if logger:
            logger.warning("No distance matrix input")
        else:
            print "No distance matrix input"
        return 1

    if method.lower() == "fastme":
        if not dist_matrix_f:
            dist_matrix_f = os.path.join(tmp_dir, "dist_matr.ph")
            dump_phylip_dist_matrix(dist_matrix=dist_matrix, matrix_path=dist_matrix_f)
        tree_f = os.path.join(tmp_dir, "tree.nwk")
        fastme_cmd = fastme_exec_path + " -i " + dist_matrix_f + " -o " + tree_f
        subprocess.call(fastme_cmd, shell=True)
        tree_newick = load_newick(newick_f_path=tree_f)
    else:
        return 1
    shutil.rmtree(tmp_dir)
    return PhyloTree(newick=tree_newick, full_seq_names=full_seq_names, config_path=config_path, logger=logger)


def compare_trees(newick1, newick2):
    pass
