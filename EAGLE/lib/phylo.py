from EAGLE.lib.general import ConfBase


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
    pass


def compare_trees(newick1, newick2):
    pass
