import os
import shutil
import subprocess

from EAGLE.lib.general import ConfBase, dump_phylip_dist_matrix


class PhyloTree(ConfBase):

    def __init__(self, newick, config_path=None, logger=None):
        self.newick = newick
        self.logger = logger
        super(PhyloTree, self).__init__(config_path=config_path)

    def update_by_config(self, config_path):
        pass


def build_tree_by_dist(dist_matrix=None,
                       dist_matrix_f=None,
                       tmp_dir="tmp",
                       method="FastME",
                       fastme_exec_path="",
                       config_path=None,
                       logger=None):

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if not dist_matrix and not dist_matrix_f:
        if logger:
            logger.warning("No distance matrix input")
        else:
            print "No distance matrix input"
        return 1
    elif not dist_matrix_f:
        dist_matrix_f = os.path.join(tmp_dir, "dist_matr.ph")
        dump_phylip_dist_matrix(dist_matrix=dist_matrix, matrix_path=dist_matrix_f)
    if method.lower() == "fastme":
        fastme_cmd = fastme_exec_path + " -i " + dist_matrix_f  ###

        subprocess.call(fastme_cmd, shell=True)

    shutil.rmtree(tmp_dir)
    return PhyloTree(newick=, config_path=config_path, logger=logger)


def compare_trees(newick1, newick2):
    pass
