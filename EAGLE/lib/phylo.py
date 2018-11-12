import os
import shutil
from copy import deepcopy
import subprocess
from collections import OrderedDict

import pandas
from Bio import Phylo

from EAGLE.constants import conf_constants, EAGLE_logger
from EAGLE.lib.general import ConfBase, filter_list


class PhyloTree(ConfBase):

    def __init__(self, tree, full_seq_names=None, tree_name="phylo_tree", tmp_dir="tmp", config_path=None, logger=None):
        self.tree = tree
        self.full_seq_names = full_seq_names
        self.tree_name = tree_name
        self.tmp_dir = tmp_dir
        self.logger = logger

        super(PhyloTree, self).__init__(config_path=config_path)

    @property
    def names(self):
        return map(lambda c: c.name, self.tree.get_terminals())

    @property
    def newick(self):
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        newick_path = os.path.join(self.tmp_dir, self.tree_name+".nwk")
        self.dump_tree(tree_path=newick_path, format="newick")
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
        return load_newick(newick_path=newick_path)

    def dump_tree(self, tree_path, format="newick"):
        Phylo.write(trees=[self.tree], file=tree_path, format=format)

    @classmethod
    def load_tree(cls,
                  tree_file,
                  tree_name="phylo_tree",
                  tmp_dir="tmp",
                  format="newick",
                  full_seq_names=None,
                  config_path=None,
                  logger=None):

        return cls(tree=Phylo.read(tree_file, format=format),
                   full_seq_names=full_seq_names,
                   tree_name=tree_name,
                   tmp_dir=tmp_dir,
                   config_path=config_path,
                   logger=logger)

    def set_full_names(self, inplace=False):
        if inplace:
            if not self.tree.rooted:
                self.tree.root_at_midpoint()
            to_remove = set(map(self._full_name_clade, self.tree.get_terminals()))
            while True:
                try:
                    self.tree.collapse(lambda c: c.name in to_remove)
                except ValueError:
                    break
        else:
            logger = self.logger
            self.logger = None
            full_names_pht = deepcopy(self)
            full_names_pht.logger = logger
            self.logger = logger
            if not full_names_pht.tree.rooted:
                full_names_pht.tree.root_at_midpoint()
            to_remove = set(map(full_names_pht._full_name_clade, full_names_pht.tree.get_terminals()))
            while True:
                try:
                    full_names_pht.tree.collapse(lambda c: c.name in to_remove)
                except ValueError:
                    break
            return full_names_pht

    def remove_names(self, names_to_remove):
        while True:
            try:
                self.tree.collapse(lambda c: c.name in names_to_remove)
            except ValueError:
                break

    def _full_name_clade(self, c):
        try:
            c.name = self.full_seq_names[c.name]
        except KeyError:
            return c.name

    def according_to_taxonomy(self, taxonomy):
        # NOT inplace method!
        pass


def build_tree_by_dist(dist_matrix=None,
                       dist_matrix_path=None,
                       full_seq_names=None,
                       tree_name="phylo_tree",
                       tmp_dir="tmp",
                       method="FastME",
                       fastme_exec_path=conf_constants.fastme_exec_path,
                       config_path=None,
                       logger=None):

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if type(dist_matrix) is not pandas.DataFrame and not dist_matrix_path:
        if logger:
            logger.warning("No distance matrix input")
        else:
            print("No distance matrix input")
        return 1
    elif dist_matrix.empty and dist_matrix_path:
        if logger:
            logger.warning("No distance matrix input")
        else:
            print("No distance matrix input")
        return 1

    if method.lower() == "fastme":
        if not dist_matrix_path:
            dist_matrix_path = os.path.join(tmp_dir, "dist_matr.ph")
            dump_phylip_dist_matrix(dist_matrix=dist_matrix, matrix_path=dist_matrix_path)
        tree_path = os.path.join(tmp_dir, "tree.nwk")
        fastme_cmd = fastme_exec_path + " -i " + dist_matrix_path + " -o " + tree_path
        subprocess.call(fastme_cmd, shell=True)
        phylo_tree = PhyloTree.load_tree(tree_file=tree_path,
                                         format="newick",
                                         full_seq_names=full_seq_names,
                                         config_path=config_path,
                                         logger=logger)
        if logger:
            logger.info("Phylogenetic tree built with FastME")
        else:
            print("Phylogenetic tree built with FastME")
    else:
        return 1
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return phylo_tree


def compare_trees(phylo_tree1,
                  phylo_tree2,
                  method="Robinson-Foulds",
                  tmp_dir="tmp",
                  emboss_inst_dir=conf_constants.emboss_inst_dir):

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    names_to_remain = set(phylo_tree1.names) & set(phylo_tree2.names)
    phylo_tree1.remove_names(set(phylo_tree1.names) - names_to_remain)
    phylo_tree2.remove_names(set(phylo_tree2.names) - names_to_remain)
    to_comp_path = os.path.join(tmp_dir, "to_comp.nwk")
    comp_res_path = os.path.join(tmp_dir, "comp.res")
    Phylo.write([phylo_tree1.tree, phylo_tree2.tree], file=to_comp_path, format="newick")
    if method.lower() in ("robinson-foulds", "rf", "s"):
        method = None
        method = "s"
    else:
        method = None
        method = "b"
    subprocess.call(os.path.join(emboss_inst_dir, "ftreedist") + " " + to_comp_path + " " +
                    comp_res_path + " -dtype " + method, shell=True)
    with open(comp_res_path) as comp_res_f:
        try:
            trees_diff = float(list(filter(lambda s: s.strip(), list(filter(lambda l: l.strip(),
                                                                            comp_res_f.readlines()))[-1].split()))[-1])
        except (ValueError, TypeError):
            trees_diff = None
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return trees_diff


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
        elif got_seqs <= num_seqs:
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
