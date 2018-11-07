import os
import shutil
import subprocess
from collections import OrderedDict

import pandas

from EAGLE.lib.general import ConfBase, filter_list


class PhyloTree(ConfBase):

    def __init__(self, newick, full_seq_names=None, config_path=None, logger=None):
        self.newick = newick
        self.full_seq_names = full_seq_names
        self.logger = logger

        super(PhyloTree, self).__init__(config_path=config_path)

    def dump_tree(self):
        pass

    @classmethod
    def load_tree(cls, newick_path, full_seq_names=None, config_path=None, logger=None):
        return cls(newick=load_newick(newick_path=newick_path),
                   full_seq_names=full_seq_names,
                   config_path=config_path,
                   logger=logger)

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
            print("No distance matrix input")
        return 1
    elif dist_matrix.empty and dist_matrix_f:
        if logger:
            logger.warning("No distance matrix input")
        else:
            print("No distance matrix input")
        return 1

    if method.lower() == "fastme":
        if not dist_matrix_f:
            dist_matrix_f = os.path.join(tmp_dir, "dist_matr.ph")
            dump_phylip_dist_matrix(dist_matrix=dist_matrix, matrix_path=dist_matrix_f)
        tree_path = os.path.join(tmp_dir, "tree.nwk")
        fastme_cmd = fastme_exec_path + " -i " + dist_matrix_f + " -o " + tree_path
        subprocess.call(fastme_cmd, shell=True)
        phylo_tree = PhyloTree.load_tree(newick_path=tree_path,
                                         full_seq_names=full_seq_names,
                                         config_path=config_path,
                                         logger=logger)
        if logger:
            logger.info("Phylogenetic tree built with FastME")
        else:
            print("Phylogenetic tree built with FastME")
    else:
        return 1
    shutil.rmtree(tmp_dir)
    return phylo_tree


def compare_trees(newick1, newick2):
    pass


def load_phylip_dist_matrix(matrix_path):
    matr_f = open(matrix_path)
    lines_dict = OrderedDict()
    seqs_list = list()
    matrix_started = False
    seq_dists_list = list()
    num_seqs = 0
    got_seqs = 0
    for line_ in matr_f:
        line = None
        line = line_.strip()
        if not line:
            continue
        line_list = filter_list(line.split())
        if len(line_list) == 1 and not matrix_started:
            num_seqs = int(line_list[0])
            continue
        if not matrix_started:
            matrix_started = True
        if got_seqs == 0:
            seqs_list.append(line_list[0])
            seq_dists_list.__iadd__(line_list[1:])
            got_seqs += len(line_list[1:])
        else:
            seq_dists_list.__iadd__(line_list)
            got_seqs += len(line_list)
        if got_seqs == num_seqs:
            lines_dict[seqs_list[-1]] = seq_dists_list
            seq_dists_list = list()
            got_seqs = 0
    dist_matrix = pandas.DataFrame.from_dict(data=lines_dict, orient='index')
    dist_matrix.columns = seqs_list
    return dist_matrix


def dump_phylip_dist_matrix(dist_matrix, matrix_path):
    matr_f = open(matrix_path, 'w')
    matr_f.write("    %s\n" % len(dist_matrix.columns))
    for seq in dist_matrix.index:
        num_spaces_to_add = 10 - len(seq)
        spaces_to_add = [" " for i in range(num_spaces_to_add)]
        matr_f.write("%s %s\n" % (seq+"".join(spaces_to_add), " ".join(dist_matrix.loc[seq].tolist())))
    matr_f.close()


def load_newick(newick_path):
    newick_f = open(newick_path)
    tree_list = list()
    for line_ in newick_f:
        line = None
        line = line_.strip()
        tree_list.append(line)
        if line[-1] == ";":
            break
    return "".join(tree_list)


def dump_tree_newick(tree_newick, newick_f_path):
    newick_f = open(newick_f_path, "w")
    newick_f.write(tree_newick)
    newick_f.close()