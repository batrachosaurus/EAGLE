import os
import shutil
from copy import deepcopy
import subprocess
from collections import OrderedDict
import numbers

import numpy as np
import pandas
import dendropy

from eagle.constants import conf_constants
from eagle.lib.general import ConfBase, filter_list, generate_random_string, send_log_message, fullmatch_regexp_list
from eagle.lib.seqs import reduce_seq_names, dump_fasta_dict


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
            return dump_tree_newick(tree_newick=tree_str, newick_f_path=tree_path)

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


class DistanceMatrix(object):
    # implement _sub_distance_matrx method

    default_format = "phylip"
    default_method = "FastME"

    def __init__(self, seqs_order, matr, short_to_full_seq_names=None, **kwargs):
        self.seqs_order = seqs_order
        assert isinstance(matr, np.ndarray), "ERROR: value for 'matr' argument should be numpy.ndarray"
        self.matr = matr
        self.short_to_full_seq_names = short_to_full_seq_names
        if self.short_to_full_seq_names is None:
            self.short_to_full_seq_names = dict()

        self.aln_name = kwargs.get("aln_name", None)
        self.aln_type = kwargs.get("aln_type", None)
        self.calc_method = kwargs.get("calc_method", None)

        self.emboss_inst_dir = kwargs.get("emboss_inst_dir", conf_constants.emboss_inst_dir)
        self.fastme_exec_path = kwargs.get("fastme_exec_path", conf_constants.fastme_exec_path)
        self.tmp_dir = kwargs.get("tmp_dir", generate_random_string(10) + "_dm_tmp")
        self.logger = kwargs.get("logger", None)
        self.config_path = kwargs.get("config_path", None)

    @property
    def seq_names(self):
        return [seq_name for seq_name, seq_i in sorted(self.seqs_order.items(), key=lambda si: si[1])]

    @property
    def full_to_short_seq_names(self):
        if not self.short_to_full_seq_names:
            self.short_to_full_seq_names = reduce_seq_names({seq_name: "" for seq_name in self.seq_names},
                                                            num_letters=10)[1]
        return {full_name: short_name for short_name, full_name in self.short_to_full_seq_names.items()}

    def __getitem__(self, item):
        if type(item) in (list, set):
            seq_names = sorted(item, key=lambda i: self.seqs_order[i])
            matr = list()
            short_to_full_seq_names = dict()
            full_to_short_seq_names = self.full_to_short_seq_names
            for seq_name in seq_names:
                short_to_full_seq_names[full_to_short_seq_names[seq_name]] = seq_name
                matr.append(list(self[seq_name][seq_names]))

            return DistanceMatrix(seqs_order={seq_name: i for i, seq_name in enumerate(seq_names)},
                                  matr=np.array(matr),
                                  short_to_full_seq_names=short_to_full_seq_names,
                                  aln_type=self.aln_type,
                                  calc_method=self.calc_method,
                                  emboss_inst_dir=self.emboss_inst_dir,
                                  fastme_exec_path=self.fastme_exec_path,
                                  logger=self.logger,
                                  config_path=self.config_path)
        else:
            return pandas.Series(self.matr[self.seqs_order[item]], index=self.seq_names)

    def __setitem__(self, key, value):
        # WARNING: new key cannot be set (key must be in self.seq_names)
        assert isinstance(value, pandas.Series), \
            "Error: value must be a pandas.Series object with org_names as indices"
        for seq_name in value.index:
            if seq_name == key:
                self.matr[self.seqs_order[key]] = np.array([value[seq_name_] for seq_name_ in self.seq_names])
            else:
                self.matr[self.seqs_order[seq_name]][self.seqs_order[key]] = value[seq_name]

    def __iter__(self):
        for seq_name in self.seq_names:
            yield seq_name

    @property
    def mean_dist(self):
        return np.mean(self.matr[self.matr >= 0.0])

    @property
    def median_dist(self):
        return np.median(self.matr[self.matr >= 0.0])

    @property
    def nseqs(self):
        return len(self.seq_names)

    @classmethod
    def calculate(cls, mult_aln, method=default_method, options=None, emboss_inst_dir=None, fastme_exec_path=None, **kwargs):
        # if method == "FastME":
        #     only '--dna' or '--protein' (full parameter names) can be used as keys in "options" dict to set the model
        if emboss_inst_dir is None:
            emboss_inst_dir = conf_constants.emboss_inst_dir
        if fastme_exec_path is None:
            fastme_exec_path = conf_constants.fastme_exec_path
        if options is None:
            options = dict()

        from eagle.lib.alignment import MultAln
        assert isinstance(mult_aln, MultAln), \
            "Error: the value for mult_aln argument is not eagle.lib.alignment.MultAln object"
        if not os.path.exists(mult_aln.tmp_dir):
            os.makedirs(mult_aln.tmp_dir)
        if method.lower() == "phylip":
            aln_fasta_path = os.path.join(mult_aln.tmp_dir, mult_aln.aln_name+".fasta")
            phylip_matrix_path = os.path.join(mult_aln.tmp_dir, mult_aln.aln_name+".phylip")
            if not mult_aln.mult_aln_dict:
                if mult_aln.logger:
                    mult_aln.logger.warning("No sequences in alignment")
                else:
                    print("No sequences in alignment")
                return
            dump_fasta_dict(fasta_dict=mult_aln.mult_aln_dict_short_id, fasta_path=aln_fasta_path)
            if mult_aln.aln_type.lower() in mult_aln.prot_types:
                if mult_aln.logger:
                    mult_aln.logger.info("protdist is starting")
                else:
                    print("protdist is starting")
                phylip_cmd = os.path.join(emboss_inst_dir, "fprotdist") + " -sequence " + aln_fasta_path + \
                             " -method d -outfile " + phylip_matrix_path
            else:
                if mult_aln.logger:
                    mult_aln.logger.info("dnadist is starting")
                else:
                    print("dnadist is starting")
                phylip_cmd = os.path.join(emboss_inst_dir, "fdnadist") + " -sequence " + aln_fasta_path + \
                             " -method f -outfile " + phylip_matrix_path
            subprocess.call(phylip_cmd, shell=True)
            if mult_aln.logger:
                mult_aln.logger.info("distance calculations finished")
            else:
                print("distance calculations finished")
            distance_matrix = cls.load(matrix_path=phylip_matrix_path,
                                       matr_format="phylip",
                                       short_to_full_seq_names=mult_aln.short_to_full_seq_names,
                                       emboss_inst_dir=emboss_inst_dir,
                                       fastme_exec_path=fastme_exec_path,
                                       **kwargs)
        if method.lower() == "fastme":
            aln_phylip_path = os.path.join(mult_aln.tmp_dir, mult_aln.aln_name+".phylip")
            phylip_matrix_path = os.path.join(mult_aln.tmp_dir, mult_aln.aln_name+"_dm.phylip")
            if not mult_aln.mult_aln_dict:
                send_log_message(message="no sequences in alignment", mes_type="e", logger=mult_aln.logger)
                return
            mult_aln.dump_alignment(aln_path=aln_phylip_path, aln_format="phylip")
            for option in ("-O", "--output_matrix"):
                options.pop(option, None)
            if mult_aln.aln_type.lower() in mult_aln.prot_types:
                fastme_options = {"--protein": "L", "-O": phylip_matrix_path}
            else:
                fastme_options = {"--dna": 4, "-O": phylip_matrix_path}
            fastme_options.update(options)
            send_log_message(message="distances calculation started", mes_type="i", logger=mult_aln.logger)
            run_fastme(input_data=aln_phylip_path, options=fastme_options, fastme_exec_path=fastme_exec_path)
            send_log_message(message="distances calculation finished", mes_type="i", logger=mult_aln.logger)

            distance_matrix = cls.load(matrix_path=phylip_matrix_path,
                                       matr_format="phylip",
                                       short_to_full_seq_names=mult_aln.short_to_full_seq_names,
                                       emboss_inst_dir=emboss_inst_dir,
                                       fastme_exec_path=fastme_exec_path,
                                       **kwargs)

        shutil.rmtree(mult_aln.tmp_dir)
        distance_matrix.aln_name = mult_aln.aln_name
        distance_matrix.aln_type = mult_aln.aln_type
        return distance_matrix

    @classmethod
    def load(cls,
             matrix_path,
             matr_format=default_format,
             short_to_full_seq_names=None,
             **kwargs):

        if matr_format == "phylip":
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
                    lines_dict[seqs_list[-1]] = list(map(float, seq_dists_list))
                    seq_dists_list = list()
                    got_seqs = 0
            matr_f.close()

            seqs_order = dict()
            matr = list()
            i = 0
            for seq_id in lines_dict:
                matr.append(lines_dict[seq_id])
                if short_to_full_seq_names:
                    seqs_order[short_to_full_seq_names[seq_id]] = i
                else:
                    seqs_order[seq_id] = i
                i += 1

            distance_matrix = cls(seqs_order=seqs_order,
                                  matr=np.array(matr),
                                  short_to_full_seq_names=short_to_full_seq_names,
                                  **kwargs)

        return distance_matrix

    def dump(self, matrix_path, matr_format=default_format):
        float_to_str = lambda x: "%f" % x if isinstance(x, numbers.Number) else str(x)
        if matr_format == "phylip":
            matr_f = open(matrix_path, 'w')
            matr_f.write("    %s\n" % self.nseqs)
            full_to_short_seq_names = self.full_to_short_seq_names
            for seq_name in self.seq_names:
                num_spaces_to_add = 10 - len(full_to_short_seq_names[seq_name])
                spaces_to_add = [" " for i in range(num_spaces_to_add)]
                matr_f.write("%s %s\n" % (full_to_short_seq_names[seq_name] + "".join(spaces_to_add),
                                          " ".join(map(float_to_str, self[seq_name]))))
            matr_f.close()
        return self.short_to_full_seq_names

    def set_new_seq_names(self, old_to_new_seq_names):
        if self.short_to_full_seq_names:
            full_to_short_seq_names = self.full_to_short_seq_names
        for seq_old_name in old_to_new_seq_names:
            self.seqs_order[old_to_new_seq_names[seq_old_name]] = self.seqs_order.pop(seq_old_name)
            if self.short_to_full_seq_names:
                self.short_to_full_seq_names[full_to_short_seq_names[seq_old_name]] = old_to_new_seq_names[seq_old_name]

    def replace_negative(self, value=None, inplace=False):
        if value is None:
            value = self.mean_dist
        if inplace:
            self.matr[self.matr < 0.0] = value
        else:
            matr = self.matr.copy()
            matr[matr < 0.0] = value
            return DistanceMatrix(seqs_order=self.seqs_order,
                                  matr=matr,
                                  short_to_full_seq_names=self.short_to_full_seq_names,
                                  aln_name=self.aln_name,
                                  aln_type=self.aln_type,
                                  calc_method=self.calc_method,
                                  emboss_inst_dir=self.emboss_inst_dir,
                                  logger=self.logger,
                                  config_path=self.config_path)

    def build_tree(self,
                   tree_name="phylo_tree",
                   tmp_dir="phylo_tree_tmp",
                   method="FastME",
                   options=None,
                   fastme_exec_path=None):

        if fastme_exec_path is None:
            fastme_exec_path = self.fastme_exec_path

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        if method.lower() == "fastme":
            dist_matrix_path = os.path.join(tmp_dir, "dist_matr.ph")
            self.replace_negative(inplace=False).dump(matrix_path=dist_matrix_path, matr_format="phylip")
            tree_path = os.path.join(tmp_dir, "tree.nwk")
            phylo_tree = run_fastme(input_data=dist_matrix_path,
                                    output_tree=tree_path,
                                    options=options,
                                    fastme_exec_path=fastme_exec_path)[0]
            phylo_tree.tree_name = tree_name
            phylo_tree.full_seq_names = self.short_to_full_seq_names
            phylo_tree.tmp_dir = tmp_dir
            phylo_tree.config_path = self.config_path
            phylo_tree.logger = self.logger
            send_log_message(message="phylogenetic tree built with FastME", mes_type="i", logger=self.logger)
        else:
            return
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return phylo_tree


def build_tree_by_aln(mult_aln=None,
                      mult_aln_fasta=None,
                      tree_name="phylo_tree",
                      tmp_dir="phylo_tree_tmp",
                      method="FastME",
                      options=None,
                      fastme_exec_path=None,
                      config_path=None,
                      logger=None):

    from eagle.lib.alignment import MultAln

    if mult_aln is None and mult_aln_fasta is not None:
        mult_aln = MultAln.load_alignment(aln_path=mult_aln_fasta, config_path=config_path, logger=logger)
    if not isinstance(mult_aln, MultAln):
        send_log_message(message="the value for mult_aln argument is not eagle.lib.alignment.MultAln object",
                         mes_type='e', logger=logger)
        return

    return mult_aln.build_tree(tree_name=tree_name, tmp_dir=tmp_dir, method=method, options=options,
                               fastme_exec_path=fastme_exec_path)


def build_tree_by_dist(dist_matrix=None,
                       dist_matrix_path=None,  # should be in phylip format
                       full_seq_names=None,
                       tree_name="phylo_tree",
                       tmp_dir="phylo_tree_tmp",
                       method="FastME",
                       options=None,
                       fastme_exec_path=None,
                       config_path=None,
                       logger=None):

    if dist_matrix is None and dist_matrix_path is not None:
        dist_matrix = DistanceMatrix.load(matrix_path=dist_matrix_path, short_to_full_seq_names=full_seq_names)
    if not isinstance(dist_matrix, DistanceMatrix):
        send_log_message(message="no distance matrix input or "
                                 "value for dist_matrix argument is not eagle.lib.alignment.DistanceMatrix object",
                         mes_type="e", logger=logger)
        return

    return dist_matrix.build_tree(tree_name=tree_name, tmp_dir=tmp_dir, method=method, options=options,
                                  fastme_exec_path=fastme_exec_path)


def run_fastme(input_data, output_tree=None, options=None, fastme_exec_path=None, **kwargs):
    if fastme_exec_path is None:
        fastme_exec_path = conf_constants.fastme_exec_path

    if options is None:
        options = dict()
    for option in ("-i", "--input_data", "-o", "--output_tree", "-n", "--nni", "-c"):
        options.pop(option, None)
    options["-i"] = input_data
    if output_tree is not None:
        options["-o"] = output_tree
    else:
        options["-c"] = True
    output_matrix = options.get("-O", options.get("--output_matrix", None))
    if not kwargs.get("no_nni", False):
        if not list(filter(None, fullmatch_regexp_list("-n.*", options))) \
               and not list(filter(None, fullmatch_regexp_list("--nni.*", options))):
            options["-n"] = True

    filter_opt = lambda opt: [opt[0]] if opt[1] is True else [opt[0], str(opt[1])]
    prepare_opt = lambda opt: " " + "=".join(opt) if "--" in opt[0] else " " + " ".join(opt)
    fastme_cmd = fastme_exec_path + "".join(prepare_opt(filter_opt(opt)) for opt in options.items())
    subprocess.call(fastme_cmd, shell=True)

    phylo_tree = None
    dist_matrix = None
    if output_tree is not None:
        phylo_tree = PhyloTree.load_tree(tree_path=output_tree, tree_format="newick")
    if output_matrix is not None:
        dist_matrix = DistanceMatrix.load(matrix_path=output_matrix, matr_format="phylip")  # actually it's not phylip
    return phylo_tree, dist_matrix


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
    return newick_f_path
