import os
import shutil
from copy import deepcopy
import subprocess
from collections import OrderedDict

import pandas
import dendropy

from eagle.constants import conf_constants
from eagle.lib.general import ConfBase, filter_list, generate_random_string
from eagle.lib.alignment import DistanceMatrix


class PhyloTree(ConfBase):

    def __init__(self,
                 tree,
                 full_seq_names=None,
                 tree_name="phylo_tree",
                 tmp_dir=None,
                 config_path=None,
                 logger=None):

        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_phylo_tree_tmp"
        if full_seq_names is None:
            full_seq_names = dict()

        self.tree = tree
        self.full_seq_names = full_seq_names
        self.tree_name = tree_name
        self.tmp_dir = tmp_dir
        self.logger = logger

        super(PhyloTree, self).__init__(config_path=config_path)

    def copy(self):
        return self._sub_phylo_tree(tree=self.tree,
                                    tree_name=self.tree_name)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self._sub_phylo_tree(tree=deepcopy(self.tree),
                                    full_seq_names=self.full_seq_names.copy(),
                                    tree_name=self.tree_name)

    def _sub_phylo_tree(self,
                        tree,
                        full_seq_names=None,
                        tree_name="phylo_tree",
                        tmp_dir=None,
                        config_path=None,
                        logger=None):

        if full_seq_names is None:
            full_seq_names = self.full_seq_names
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_phylo_tree_tmp"
        if config_path is None:
            config_path = self.config_path
        if logger is None:
            logger = self.logger

        return PhyloTree(tree=tree,
                         full_seq_names=full_seq_names,
                         tree_name=tree_name,
                         tmp_dir=tmp_dir,
                         config_path=config_path,
                         logger=logger)

    @property
    def names(self):
        return filter(None, map(lambda node: node.taxon.label if node.taxon else None, self.tree.nodes()))

    @property
    def newick(self):
        return self.tree.as_string(schema="newick").replace(" ", "_").strip()

    def dump_tree(self, tree_path, tree_format="newick"):
        tree_str = self.tree.as_string(schema=tree_format).replace(" ", "_")
        if tree_format == "newick":
            dump_tree_newick(tree_newick=tree_str, newick_f_path=tree_path)

    @classmethod
    def load_tree(cls,
                  tree_path,
                  tree_name="phylo_tree",
                  tmp_dir="tmp",
                  tree_format="newick",
                  full_seq_names=None,
                  config_path=None,
                  logger=None):

        if tree_format == "newick":
            tree = load_newick(tree_path)
        else:
            print("loading trees from %s not implemented yet - use load_tree_from_str method" % tree_format)
            return

        return cls(tree=dendropy.Tree.get_from_string(tree, schema=tree_format),
                   full_seq_names=full_seq_names,
                   tree_name=tree_name,
                   tmp_dir=tmp_dir,
                   config_path=config_path,
                   logger=logger)

    @classmethod
    def load_tree_from_str(cls,
                           tree_str,
                           tree_name="phylo_tree",
                           tmp_dir="tmp",
                           tree_format="newick",
                           full_seq_names=None,
                           config_path=None,
                           logger=None):

        return cls(tree=dendropy.Tree.get_from_string(tree_str, schema=tree_format),
                   full_seq_names=full_seq_names,
                   tree_name=tree_name,
                   tmp_dir=tmp_dir,
                   config_path=config_path,
                   logger=logger)

    def set_full_names(self, inplace=False):
        if inplace:
            if not self.full_seq_names:
                return
            names_to_remove = list()
            for name in self.names:
                node_to_rename = self.tree.find_node_with_taxon_label(name)
                try:
                    node_to_rename.taxon.label = self.full_seq_names[name.replace(" ", "_")]
                except KeyError:
                    names_to_remove.append(name)
            self.tree.prune_taxa_with_labels(names_to_remove)
        else:
            full_names_pht = deepcopy(self)
            full_names_pht.set_full_names(inplace=True)
            return full_names_pht

    def remove_names(self, names_to_remove, inplace=False):
        if inplace:
            self.tree.prune_taxa_with_labels(names_to_remove)
        else:
            rem_names_pht = deepcopy(self)
            rem_names_pht.remove_names(names_to_remove, inplace=True)
            return rem_names_pht

    def according_to_taxonomy(self, taxonomy):
        # NOT inplace method!
        pass


def build_tree_by_dist(dist_matrix=None,
                       dist_matrix_path=None,
                       full_seq_names=None,
                       tree_name="phylo_tree",
                       tmp_dir="phylo_tree_tmp",
                       method="FastME",
                       fastme_exec_path=None,
                       config_path=None,
                       logger=None):

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if fastme_exec_path is None:
        fastme_exec_path = conf_constants.fastme_exec_path
    if not isinstance(dist_matrix, DistanceMatrix) and not dist_matrix_path:
        if logger:
            logger.warning("No distance matrix input or "
                           "value for dist_matrix argument is not eagle.lib.alignment.DistanceMatrix object")
        else:
            print("No distance matrix input or"
                  "value for dist_matrix argument is not eagle.lib.alignment.DistanceMatrix object")
        return 1

    if method.lower() == "fastme":
        if not dist_matrix_path:
            dist_matrix_path = os.path.join(tmp_dir, "dist_matr.ph")
            dist_matrix.replace_negative(inplace=False).dump(matrix_path=dist_matrix_path, matr_format="phylip")
        tree_path = os.path.join(tmp_dir, "tree.nwk")
        fastme_cmd = fastme_exec_path + " -i " + dist_matrix_path + " -o " + tree_path
        subprocess.call(fastme_cmd, shell=True)
        phylo_tree = PhyloTree.load_tree(tree_path=tree_path,
                                         tree_format="newick",
                                         full_seq_names=full_seq_names,
                                         tree_name=tree_name,
                                         tmp_dir=tmp_dir,
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
                  method="Robinson-Foulds"):

    taxa = dendropy.TaxonNamespace()
    names_to_remain = set(phylo_tree1.names) & set(phylo_tree2.names)
    pht1 = dendropy.Tree.get_from_string(phylo_tree1.remove_names(list(set(phylo_tree1.names)-names_to_remain)).newick,
                                         schema='newick',
                                         taxon_namespace=taxa)
    pht2 = dendropy.Tree.get_from_string(phylo_tree2.remove_names(list(set(phylo_tree2.names)-names_to_remain)).newick,
                                         schema='newick',
                                         taxon_namespace=taxa)
    if method.lower() in ("robinson-foulds", "rf") and len(pht1.leaf_nodes()) > 3:
        pht1.encode_bipartitions()
        pht2.encode_bipartitions()
        trees_diff = float(dendropy.treecalc.treecompare.symmetric_difference(
            tree1=pht1,
            tree2=pht2,
            is_bipartitions_updated=True,
        )) / float(2*len(pht1.leaf_nodes())-6)
    else:
        return
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
