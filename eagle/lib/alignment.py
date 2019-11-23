# Consider create package alignment and place the objects into separate modules
import os
import re
import shutil
import subprocess
import multiprocessing as mp
from copy import copy, deepcopy
from collections import defaultdict, OrderedDict

import numpy as np
from scipy.stats import gmean
import pandas
from Bio.Seq import Seq

from eagledb.scheme import SeqProfileInfo
from eagle.constants import conf_constants
from eagle.lib.general import ConfBase, join_files, generate_random_string, send_log_message, fullmatch_regexp_list
from eagle.lib.phylo import DistanceMatrix, run_fastme
from eagle.lib.seqs import load_fasta_to_dict, dump_fasta_dict, reduce_seq_names, SeqsDict, detect_seqs_type, \
    shred_fasta


class MultAln(ConfBase):
    # Consider inherit it from SeqsDict and making ConfBase the mixin

    prot_type = SeqsDict.prot_type
    nucl_type = SeqsDict.nucl_type
    prot_types = SeqsDict.prot_types
    nucl_types = SeqsDict.nucl_types

    dist_matr_method = "FastME"
    low_memory = False

    def __init__(self,
                 mult_aln_dict=None,
                 seq_info_dict=None,
                 aln_type=None,
                 states_seq=None,
                 aln_name="mult_aln",
                 tmp_dir=None,
                 emboss_inst_dir=None,
                 hmmer_inst_dir=None,
                 fastme_exec_path=None,
                 kaks_calculator_exec_path=None,
                 config_path=None,
                 logger=None):

        if mult_aln_dict is None:
            mult_aln_dict = dict()
        if seq_info_dict is None:
            seq_info_dict = defaultdict(lambda: defaultdict(lambda: None))
        if states_seq is None:
            states_seq = list()
        if emboss_inst_dir is None:
            emboss_inst_dir = conf_constants.emboss_inst_dir
        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if fastme_exec_path is None:
            fastme_exec_path = conf_constants.fastme_exec_path
        if kaks_calculator_exec_path is None:
            kaks_calculator_exec_path = conf_constants.kaks_calculator_exec_path

        self._nucl_seqs_dict = dict()
        self._short_to_full_seq_names = dict()
        self._seq_ids_to_orgs = dict()

        self.mult_aln_dict = mult_aln_dict
        self.seq_info_dict = seq_info_dict
        self.states_seq = states_seq
        self._aln_type = aln_type
        self._distance_matrix = None
        self.aln_name = aln_name
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_mult_aln_tmp"
        self.tmp_dir = tmp_dir
        self.emboss_inst_dir = emboss_inst_dir
        self.hmmer_inst_dir = hmmer_inst_dir
        self.fastme_exec_path = fastme_exec_path
        self.kaks_calculator_exec_path=kaks_calculator_exec_path
        self.logger = logger
        super(MultAln, self).__init__(config_path=config_path)

    @property
    def nucl_seqs_dict(self):
        return self._nucl_seqs_dict

    @nucl_seqs_dict.setter
    def nucl_seqs_dict(self, nucl_seqs_dict):
        if self.aln_type.lower() not in self.prot_types:
            send_log_message(message="The alignment seems to be nucleotide - not applicable",
                             mes_type="w", logger=self.logger)
        else:
            if not isinstance(nucl_seqs_dict, SeqsDict):
                nucl_seqs_dict = SeqsDict.load_from_dict(nucl_seqs_dict)
            self._nucl_seqs_dict = nucl_seqs_dict

    @property
    def short_to_full_seq_names(self):
        return self._short_to_full_seq_names

    @short_to_full_seq_names.setter
    def short_to_full_seq_names(self, stfsn_dict):
        if isinstance(stfsn_dict, dict):  # not exhaustive condition
            self._short_to_full_seq_names = stfsn_dict
        else:
            send_log_message(message="the value to assign should be a dict", mes_type="e", logger=self.logger)

    @property
    def seq_ids_to_orgs(self):
        return self._seq_ids_to_orgs

    @seq_ids_to_orgs.setter
    def seq_ids_to_orgs(self, sito_dict):
        if isinstance(sito_dict, dict):  # not exhaustive condition
            self._seq_ids_to_orgs = sito_dict
        else:
            send_log_message(message="the value to assign should be a dict", mes_type="e", logger=self.logger)

    def __len__(self):
        return len(self.mult_aln_dict)

    def __getitem__(self, item):
        if type(item) is slice:
            states_seq = self.states_seq
            if states_seq:
                states_seq = states_seq[item]
            return self._sub_mult_aln(mult_aln_dict=SeqsDict.load_from_dict({seq: self[seq][item] for seq in self}),
                                      aln_type=self.aln_type,
                                      aln_name=self.aln_name+"[%s:%s]" % (item.start, item.stop),
                                      states_seq=states_seq)
        elif type(item) in (list, set) and isinstance(self.mult_aln_dict, SeqsDict):
            return self._sub_mult_aln(self.mult_aln_dict.get_sample(seqs=item, low_memory=self.low_memory),
                                      aln_type=self.aln_type)
        else:
            return self.mult_aln_dict[item]

    def __setitem__(self, seq_id, seq):
        self.mult_aln_dict[seq_id] = seq

    def __delitem__(self, seq_id):
        del self.short_to_full_seq_names[self.full_to_short_seq_names[seq_id]]
        del self.mult_aln_dict[seq_id]
        if self._distance_matrix:
            self._distance_matrix = None

    def pop(self, seq_id):
        del self.short_to_full_seq_names[self.full_to_short_seq_names[seq_id]]
        seq = self.mult_aln_dict.pop(seq_id)
        if self._distance_matrix:
            self._distance_matrix = None
        return seq

    def __iter__(self):
        for seq in self.seq_names:
            yield seq

    def _sub_mult_aln(self,
                      mult_aln_dict,
                      seq_info_dict=None,
                      aln_type=None,
                      states_seq=None,
                      aln_name="mult_aln",
                      tmp_dir=None,
                      short_to_full_seq_names=None,
                      emboss_inst_dir=None,
                      fastme_exec_path=None,
                      hmmer_inst_dir=None,
                      kaks_calculator_exec_path=None,
                      config_path=None,
                      logger=None):

        if seq_info_dict is None:
            seq_info_dict = self.seq_info_dict
        if aln_type is None:
            aln_type = self.aln_type
        if states_seq is None:
            states_seq = list()
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_mult_aln_tmp"
        if short_to_full_seq_names is None:
            short_to_full_seq_names = self.short_to_full_seq_names
        if emboss_inst_dir is None:
            emboss_inst_dir = self.emboss_inst_dir
        if fastme_exec_path is None:
            fastme_exec_path = self.fastme_exec_path
        if hmmer_inst_dir is None:
            hmmer_inst_dir = self.hmmer_inst_dir
        if kaks_calculator_exec_path is None:
            kaks_calculator_exec_path = self.kaks_calculator_exec_path
        if config_path is None:
            config_path = self.config_path
        if logger is None:
            logger = self.logger

        mult_aln = MultAln(mult_aln_dict=mult_aln_dict,
                           seq_info_dict=seq_info_dict,
                           aln_type=aln_type,
                           aln_name=aln_name,
                           tmp_dir=tmp_dir,
                           states_seq=states_seq,
                           emboss_inst_dir=emboss_inst_dir,
                           fastme_exec_path=fastme_exec_path,
                           hmmer_inst_dir=hmmer_inst_dir,
                           kaks_calculator_exec_path=kaks_calculator_exec_path,
                           config_path=config_path,
                           logger=logger)
        mult_aln.short_to_full_seq_names = short_to_full_seq_names
        mult_aln._seq_ids_to_orgs = self._seq_ids_to_orgs
        return mult_aln

    def copy(self):
        return self._sub_mult_aln(mult_aln_dict=self.mult_aln_dict,
                                  aln_name=self.aln_name,
                                  states_seq=self.states_seq)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self._sub_mult_aln(mult_aln_dict=deepcopy(self.mult_aln_dict),
                                  aln_name=self.aln_name,
                                  short_to_full_seq_names=self.short_to_full_seq_names.copy(),
                                  states_seq=copy(self.states_seq))  # may be deepcopy is needed

    def __contains__(self, item):
        return item in self.mult_aln_dict

    @property
    def seq_names(self):
        if self.mult_aln_dict:
            return list(self.mult_aln_dict.keys())
        else:
            return list()

    @property
    def name(self):
        return self.aln_name

    @property
    def aln_type(self):
        if self._aln_type is None:
            self._aln_type = detect_seqs_type(seqs_dict=self.mult_aln_dict)
        return self._aln_type

    @aln_type.setter
    def aln_type(self, aln_type):
        self._aln_type = aln_type

    def rename_seqs(self, old_to_new_dict):
        if isinstance(self.mult_aln_dict, SeqsDict):
            self.mult_aln_dict.rename_seqs(old_to_new_dict)
        else:
            for old_seq in old_to_new_dict:
                if old_seq in self.mult_aln_dict:
                    self.mult_aln_dict[old_to_new_dict[old_seq]] = self.mult_aln_dict.pop(old_seq)
        for short_name in list(self.short_to_full_seq_names.keys()):
            self.short_to_full_seq_names[short_name] = old_to_new_dict.get(self.short_to_full_seq_names[short_name],
                                                                           self.short_to_full_seq_names[short_name])
        self._distance_matrix = None

    @property
    def length(self):
        if self.mult_aln_dict:
            return len(self.mult_aln_dict[self.seq_names[0]])
        else:
            return 0

    @property
    def num_seqs(self):
        return len(self.mult_aln_dict)

    @property
    def shape(self):
        return self.num_seqs, self.length

    @property
    def mult_aln_dict_short_id(self):
        if self.mult_aln_dict:
            if self.short_to_full_seq_names:
                full_to_short_seq_names = self.full_to_short_seq_names
                return dict((full_to_short_seq_names[seq_id], seq) for seq_id, seq in self.mult_aln_dict.items())
            else:
                mult_aln_dict_short_ids, short_to_full_seq_names = reduce_seq_names(fasta_dict=self.mult_aln_dict,
                                                                                    num_letters=10)
                self.short_to_full_seq_names = short_to_full_seq_names
            return mult_aln_dict_short_ids
        else:
            return dict()

    @property
    def full_to_short_seq_names(self):
        if self.short_to_full_seq_names:
            seq_names = self.seq_names
            return dict(filter(None, map(lambda x: (x[1], x[0]) if x[1] in seq_names else None,
                                         self.short_to_full_seq_names.items())))
        else:
            return dict()

    def dump_alignment(self, aln_path, aln_format="fasta"):
        if aln_format.lower() == "fasta":
            dump_fasta_dict(fasta_dict=self.mult_aln_dict, fasta_path=aln_path)
        if aln_format.lower() == "phylip":
            mult_aln_dict_short_id = self.mult_aln_dict_short_id
            with open(aln_path, "w") as aln_f:
                aln_f.write("    %s    %s\n" % self.shape)
                for seq_name in mult_aln_dict_short_id:
                    num_spaces_to_add = 10 - len(seq_name)
                    spaces_to_add = [" " for i in range(num_spaces_to_add)]
                    aln_f.write("%s %s\n" % (seq_name+"".join(spaces_to_add), mult_aln_dict_short_id[seq_name]))
        return self.short_to_full_seq_names

    @classmethod
    def load_alignment(cls, aln_path=None, aln_format="fasta", aln_type=None, aln_name=None,
                       config_path=None, logger=None, **kwargs):
        # For backward compatibility
        if aln_path is None:
            aln_path = kwargs.get("aln_fasta_path", None)

        if aln_name is None:
            aln_name = os.path.splitext(aln_path)[0]
        if aln_format.lower() in ("fasta", "fas"):
            mult_aln_dict = load_fasta_to_dict(fasta_path=aln_path, **kwargs)
        else:
            print("ERROR: only fasta supported yet")
            return
        return cls(mult_aln_dict=mult_aln_dict,
                   aln_type=aln_type,
                   aln_name=aln_name,
                   config_path=config_path,
                   logger=logger)

    def improve_aln(self,
                    max_gap_fract=0.75,  # maximal fraction of gaps in a column to keep it in alignment
                    max_mismatch_fract=1.0,  # maximal fraction of mismatches in a column to keep it in alignment
                    remove_seq=False,  # if True can remove sequences for cleaning gaps between aln blocks
                    dist_filt=False,  # if True returns a list of alignments (or blocks coordinates) with different strictness of sequences removing (different thresholds for clustering)
                    output='alignment',
                    inplace=False,
                    **kwargs):
        # TODO: write methods for removing gaps between blocks
        if self.logger:
            self.logger.info("'%s' multiple alignment improvement started" % self.aln_name)
        else:
            print("'%s' multiple alignment improvement started" % self.aln_name)
        
        if dist_filt:
            # Constructs a list of rarefied by different distance thresholds alignments and runs improve_aln on it with dist_filt=False
            pass
        coords = list()
        seq_id_list = list(self.mult_aln_dict.keys())
        for i in range(self.length):
            if self.states_seq:
                if self.states_seq[i] in ("match", "specific"):
                    coords = self._update_coords(i+1, coords)
                    continue
            let_counts = defaultdict(int)
            let_counts["-"] = 0
            for seq_id in seq_id_list:
                let_counts[self.mult_aln_dict[seq_id][i]] += 1
            if float(let_counts.pop("-"))/float(len(seq_id_list)) <= max_gap_fract:
                if sorted(let_counts.values(), reverse=True)[0] >= 1.0-max_mismatch_fract:
                    coords = self._update_coords(i+1, coords)
            elif remove_seq:
                # TODO: write detection of seqs to remove
                pass
        if len(coords[-1]) == 1:
            coords[-1] = (coords[-1][0], coords[-1][0])
        if not remove_seq and len(coords) >= 1:
            coords = [(coords[0][0], coords[-1][1])]
        if output.lower() in ("coordinates", "coords", "coord", "c"):
            seq_coord_list = list()
            for seq_id in self.seq_names:
                seq_c_dict = {"seq_id": seq_id}
                coords_for_seq = self._get_coords_for_seq(coords, self.mult_aln_dict[seq_id])
                for k in range(len(coords_for_seq)):
                    seq_c_dict["c%s" % ((k+1)*2-1)] = coords_for_seq[k][0]
                    seq_c_dict["c%s" % ((k+1)*2)] = coords_for_seq[k][1]
                seq_coord_list.append(seq_c_dict)
            return pandas.DataFrame(seq_coord_list)
        else:
            if inplace:
                seqs_ids = self.seq_names
                for seq_id in seqs_ids:
                    self.mult_aln_dict[seq_id] = "".join(self.mult_aln_dict[seq_id][c1: c2] for c1, c2 in coords)
                self._distance_matrix = None
                if self.logger:
                    self.logger.info("'%s' multiple alignment improved" % self.aln_name)
                else:
                    print("'%s' multiple alignment improved" % self.aln_name)
            else:
                impr_mult_aln_dict = {seq_id: "".join(self[seq_id][c1:c2] for c1, c2 in coords)
                                      for seq_id in self.seq_names}
                impr_mult_aln = self._sub_mult_aln(mult_aln_dict=SeqsDict.load_from_dict(impr_mult_aln_dict),
                                                   aln_type=self.aln_type,
                                                   aln_name="impr_"+self.aln_name)
                if self.logger:
                    self.logger.info("'%s' multiple alignment improved" % self.aln_name)
                else:
                    print("'%s' multiple alignment improved" % self.aln_name)
                return impr_mult_aln

    @staticmethod
    def _update_coords(i, coords):
        if not coords:
            coords.append([i])
        elif len(coords[-1]) == 1:
            coords[-1].append(i)
        elif coords[-1][1] == i - 1:
            coords[-1][1] = i
        else:
            coords.append([i])
        return coords

    @staticmethod
    def _get_coords_for_seq(coords, gapped_seq):
        no_gaps_coords = list()
        for coord_pair in coords:
            c1_shift = gapped_seq[:coord_pair[0]].count("-")
            c2_shift = gapped_seq[:coord_pair[1]].count("-")
            no_gaps_coords.append((coord_pair[0]-c1_shift, coord_pair[1]-c2_shift))
        return no_gaps_coords

    @property
    def distance_matrix(self):
        if self._distance_matrix is None:
            self._distance_matrix = DistanceMatrix.calculate(mult_aln=self,
                                                             method=self.dist_matr_method,
                                                             emboss_inst_dir=self.emboss_inst_dir,
                                                             fastme_exec_path=self.fastme_exec_path,
                                                             logger=self.logger,
                                                             config_path=self.config_path)
        return self._distance_matrix

    def get_distance_matrix(self, method=None):
        if method is None:
            method = self.dist_matr_method
        self._distance_matrix = DistanceMatrix.calculate(mult_aln=self,
                                                         method=method,
                                                         emboss_inst_dir=self.emboss_inst_dir,
                                                         fastme_exec_path=self.fastme_exec_path,
                                                         logger=self.logger,
                                                         config_path=self.config_path)
        return self._distance_matrix

    def remove_paralogs(self, seq_ids_to_orgs=None, method="min_dist", inplace=False):
        """

        :param seq_ids_to_orgs: {seq_id: org_name}
        :param method:
        :param inplace:
        :return:
        """
        send_log_message(message="paralogs removing started", mes_type='i', logger=self.logger)

        if seq_ids_to_orgs is not None:
            self.seq_ids_to_orgs = seq_ids_to_orgs

        org_dist_dict = defaultdict(dict)
        for seq_name in self.distance_matrix.seq_names:
            org_dist_dict[self.seq_ids_to_orgs[seq_name]][seq_name] = sum(map(float, self.distance_matrix[seq_name]))
        org_names = org_dist_dict.keys()
        for org in org_names:
            org_dist_dict[org] = sorted(org_dist_dict[org].items(), key=lambda x: x[1])
        mult_aln_dict_filt = dict()
        short_to_full_seq_names_filt = dict()
        if method.lower() in ("minimal_distance", "min_dist", "md"):
            for seq_name in self.seq_names:
                if seq_name == org_dist_dict[self.seq_ids_to_orgs[seq_name]][0][0]:
                    mult_aln_dict_filt[seq_name] = self.mult_aln_dict[seq_name]
        else:
            # TODO: write spec_pos method
            pass
        for seq_short_name in self.short_to_full_seq_names:
            if self.short_to_full_seq_names[seq_short_name] in mult_aln_dict_filt:
                short_to_full_seq_names_filt[seq_short_name] = self.short_to_full_seq_names[seq_short_name]
        mult_aln_dict_filt = SeqsDict.load_from_dict(in_dict=mult_aln_dict_filt)
        if inplace:
            self.mult_aln_dict = mult_aln_dict_filt
            self.short_to_full_seq_names = short_to_full_seq_names_filt
            self._distance_matrix = None
            if self.logger:
                self.logger.info("paralogs removed")
            else:
                print("paralogs removed")
        else:
            filtered_aln = self._sub_mult_aln(mult_aln_dict=mult_aln_dict_filt,
                                              aln_type=self.aln_type,
                                              aln_name="no_par_"+self.aln_name,
                                              short_to_full_seq_names=short_to_full_seq_names_filt)
            if self.logger:
                self.logger.info("paralogs removed")
            else:
                print("paralogs removed")
            return filtered_aln

    def get_blocks_tsv(self, tsv_path, fasta_path, split_into_blocks=False, meta_dict=None):
        # Three block types: CB - conservetive block, SB - specific block, MB - major block. UB - unclassified block (not any of three types => MB=NA)
        cut_ends_df = self.improve_aln(output='coords')
        if split_into_blocks:
            blocks_df = self.define_aln_blocks(cut_ends_coord_dict=cut_ends_df)
        else:
            blocks_df = cut_ends_df
            blocks_df["block_type"] = pandas.Series(["UB"]*cut_ends_df.shape[0])
            blocks_df["block_descr"] = pandas.Series(['major_block_id "NA"; block_id "UB1"; block_name "NA"']*
                                                     cut_ends_df.shape[0])
            blocks_df = blocks_df.rename(index=str, columns={"c1": "start", "c2": "end"})
            blocks_df = blocks_df[["seq_id", "block_type", "start", "end", "block_descr"]]
        tsv_f = open(tsv_path, 'w')
        seqs_dict = dict(list(blocks_df.apply(self._write_blocks_tsv, axis=1, args=(tsv_f, meta_dict))))
        tsv_f.close()
        dump_fasta_dict(fasta_dict=seqs_dict, fasta_path=fasta_path)
        if self.logger:
            self.logger.info("Blocks tsv %s ready" % tsv_path)
        else:
            print("Blocks tsv %s ready" % tsv_path)

    def _write_blocks_tsv(self, block_info, tsv_f, meta_dict=None):
        tsv_f.write("\t".join(map(str, list(block_info)))+"; %s\n" %
                    "; ".join('%s "%s"' % (meta[0], meta[1]) for meta in meta_dict[block_info["seq_id"]].items()))
        return block_info["seq_id"], self.mult_aln_dict[block_info["seq_id"]].replace("-", "")

    def remove_seqs(self, seqs_list):
        pass

    def define_aln_blocks(self, cut_ends_coord_dict):
        blocks_info = list()  # list of dicts
        pass
        return pandas.DataFrame(blocks_info)

    def get_hmm_profile(self, profile_path, method="hmmer"):
        if method.lower() == "hmmer":
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)
            aln_fasta_path = os.path.join(self.tmp_dir, self.aln_name+".fasta")
            dump_fasta_dict(self.mult_aln_dict, aln_fasta_path)
            hmmer_handler = HmmerHandler(inst_dir=self.hmmer_inst_dir, config_path=self.config_path, logger=self.logger)
            hmmer_handler.build_hmm_profile(profile_path=profile_path, in_aln_path=aln_fasta_path)
            shutil.rmtree(self.tmp_dir)
            return profile_path

    def nucl_by_prot_aln(self, nucl_fasta_dict=None, nucl_fasta_path=None):
        if self.aln_type.lower() not in self.prot_types:
            if self.logger:
                self.logger.error("reference alignment type is not protein")
            else:
                print("ERROR: reference alignment type is not protein")
            return

        if nucl_fasta_dict is not None:
            self.nucl_seqs_dict = nucl_fasta_dict
        elif nucl_fasta_path is not None:
            self.nucl_seqs_dict = load_fasta_to_dict(fasta_path=nucl_fasta_path)
        if not self.nucl_seqs_dict:
            send_log_message(message="no nucleotide sequences are input", mes_type="e", logger=self.logger)
            return

        nucl_aln_dict = dict()
        for seq_name in self.seq_names:
            match = re.search(re.sub("[-\.]", "", self[seq_name]).replace("*", "."),
                              str(Seq(self.nucl_seqs_dict[seq_name]).translate()))
            nucl_aln_dict[seq_name] = nucl_accord_prot(self[seq_name],
                                                       self.nucl_seqs_dict[seq_name][match.start()*3: match.end()*3])
        nucl_aln = self._sub_mult_aln(SeqsDict.load_from_dict(nucl_aln_dict),
                                      aln_type=self.nucl_type,
                                      aln_name="nucl_"+self.aln_name)
        return nucl_aln

    def estimate_uniformity(self, cons_thr=None, window_l=None, windows_step=None):
        if cons_thr is None:
            cons_thr = conf_constants.cons_thr
        if window_l is None:
            window_l = conf_constants.unif_window_l
        if windows_step is None:
            windows_step = conf_constants.unif_windows_step

        windows_list = self.split_into_windows(window_l=window_l, windows_step=windows_step)
        cons_cols_by_windows = np.array([w.cons_cols_num(cons_thr=cons_thr) for w in windows_list])
        return np.std(cons_cols_by_windows)

    def split_into_windows(self, window_l, windows_step=None):
        if windows_step is None:
            windows_step = window_l

        windows_list = list()
        i = 0
        while i < self.length:
            windows_list.append(self[i: i + window_l])
            i += windows_step
        return windows_list

    def cons_cols_num(self, cons_thr=None):
        if cons_thr is None:
            cons_thr = conf_constants.cons_thr
        if cons_thr > 1.0:
            cons_thr = float(cons_thr) / 100.0

        cln = 0
        for i in range(self.length):
            s_num_dict = defaultdict(int)
            for seq_id in self.mult_aln_dict:
                s_num_dict[self.mult_aln_dict[seq_id][i].lower()] += 1
            all_s_num = sum(s_num_dict.values())
            if float(s_num_dict.get("-", 0)) / float(all_s_num) <= 1.0 - cons_thr:
                if float(sorted(s_num_dict.values(), reverse=True)[0]) / float(all_s_num) >= cons_thr:
                    cln += 1
        return cln

    def rarefy(self, seqs_to_remain=100):
        seqs_ids = self.mult_aln_dict.keys()
        if len(seqs_ids) <= seqs_to_remain:
            return self
        rarefied_aln_dict = dict()
        for i in range(seqs_to_remain):
            seq_id = None
            seq_id = seqs_ids.pop(np.random.randint(len(seqs_ids)))
            rarefied_aln_dict[seq_id] = self.mult_aln_dict[seq_id]
        return self._sub_mult_aln(SeqsDict.load_from_dict(rarefied_aln_dict),
                                  aln_name="raref_"+self.aln_name)

    def calculate_KaKs_windows(self,
                               nucl_seqs_dict=None,
                               window_l=None,
                               top_fract=None,
                               method="KaKs_Calculator",
                               only_first_seq=False,
                               **kwargs):
        if nucl_seqs_dict is not None:
            self.nucl_seqs_dict = nucl_seqs_dict
        if window_l is None:
            window_l = conf_constants.kaks_window_l
        if top_fract is None:
            top_fract = conf_constants.kaks_top_fract
        if self.aln_type is None:
            self.aln_type = detect_seqs_type(seqs_dict=self.mult_aln_dict)

        if self.aln_type.lower() in self.prot_types:
            if not self.nucl_seqs_dict:
                if self.logger:
                    self.logger.error("protein alignment but no value input for argument 'nucl_seqs_dict'")
                else:
                    print("ERROR: protein alignment but no value input for argument 'nucl_seqs_dict'")
                return
            nucl_aln = self.nucl_by_prot_aln()
            windows_list = nucl_aln.split_into_windows(window_l=window_l)
        else:
            windows_list = self.split_into_windows(window_l=window_l)

        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        kaks_list = list()
        ks_list = list()
        for w in windows_list:
            w.tmp_dir = self.tmp_dir
            w_kaks = w.calculate_KaKs(method=method,
                                      remove_tmp=False,
                                      only_first_seq=only_first_seq,
                                      raref_base=kwargs.pop("raref_base", 10.0),
                                      **kwargs)
            if type(w_kaks) is float or isinstance(w_kaks, np.floating):
                if w_kaks["Ka/Ks"] >= 0.0 and not pandas.isna(w_kaks["Ka/Ks"]):
                    kaks_list.append(w_kaks["Ka/Ks"])
                    ks_list.append(w_kaks["Ks"])
                else:
                    kaks_list.append(-1.0)
                    ks_list.append(-1.0)

        shutil.rmtree(self.tmp_dir)
        return {"Ka/Ks": kaks_list, "Ks": ks_list}

    def calculate_KaKs(self,
                       nucl_seqs_dict=None,
                       method="KaKs_Calculator",
                       max_kaks=10.0,
                       only_first_seq=False,
                       **kwargs):

        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if nucl_seqs_dict is not None:
            self.nucl_seqs_dict = nucl_seqs_dict

        if self.aln_type.lower() in self.prot_types:
            if not self.nucl_seqs_dict:
                send_log_message(message="protein alignment but no value input for argument 'nucl_seqs_dict'",
                                 mes_type="e",
                                 logger=self.logger)
                return
            nucl_aln = self.nucl_by_prot_aln()
            return nucl_aln.calculate_KaKs(method=method, max_kaks=max_kaks, only_first_seq=only_first_seq, **kwargs)

        seq_pairs = self.generate_seqs_pairs(only_first_seq=only_first_seq, raref_base=kwargs.get("raref_base", None))
        if method.lower() == "kaks_calculator":
            pairs_axt_path = os.path.join(self.tmp_dir, self.aln_name+".axt")
            with open(pairs_axt_path, "w") as pairs_axt_f:
                for seq_pair in seq_pairs:
                    pairs_axt_f.write("vs".join(seq_pair)+"\n")
                    pairs_axt_f.write(seq_pairs[seq_pair][0]+"\n")
                    pairs_axt_f.write(seq_pairs[seq_pair][1]+"\n\n")
            out_path = os.path.join(self.tmp_dir, self.aln_name+".kaks")
            kaks_calculator_cmd = self.kaks_calculator_exec_path + \
                                  " -i " + pairs_axt_path + " -o " + out_path + " -m YN"
            if self.logger:
                self.logger.info("running command: '%s'" % kaks_calculator_cmd)
            else:
                print("INFO: running command: '%s'" % kaks_calculator_cmd)
            subprocess.call(kaks_calculator_cmd, shell=True)
            try:
                kaks_df = pandas.read_csv(out_path, sep="\t")
                kaks = gmean(kaks_df["Ka/Ks"].apply(lambda v: min(v, max_kaks)))
                ks = gmean(kaks_df["Ks"])
            except:
                kaks = -1.0
                ks = -1.0
        if kwargs.get("remove_tmp", True):
            shutil.rmtree(self.tmp_dir)
        return {"Ka/Ks": kaks, "Ks": ks}

    def generate_seqs_pairs(self, only_first_seq=False, raref_base=None):
        """
        Generates pars of aligned sequences from multiple alignment (n(n-1)/2 by deafault)
        :param only_first_seq: set it True to get only pairs containig first sequence (default False)
        :param raref_base: randomly rarefies pairs: n(n-1)/2 -> n*raref_base/2 (10.0 makes sense to use)
        :return:
        """
        seqs_pairs = dict()
        seq_names_list = list(self.seq_names)
        for i, seqi_name in enumerate(seq_names_list[:-1]):
            for seqj_name in seq_names_list[i+1:]:
                seqs_pairs[frozenset({seqi_name, seqj_name})] = [self[seqi_name], self[seqj_name]]
            if only_first_seq:
                break
        if raref_base is not None and int(raref_base) < len(seq_names_list) and not only_first_seq:
            for seqs_pair in np.random.permutation(list(seqs_pairs.keys()))[(len(seq_names_list)-1)*int(raref_base/2.0):]:
                del seqs_pairs[seqs_pair]
        return seqs_pairs

    def stop_codons_stats(self, improve_aln=False, **kwargs):
        if improve_aln:
            improved_aln = self.improve_aln(**kwargs)
        else:
            improved_aln = self
        stops_per_seq = list()
        if improved_aln.aln_type in self.prot_types:
            for seq_name in improved_aln:
                stops_per_seq.append(improved_aln[seq_name].count("*"))
        else:
            for seq_name in improved_aln:
                stops_per_seq.append(re.sub("(tag|taa|tga)", "*", improved_aln[seq_name].lower()).count("*"))
        stops_per_seq.sort()
        return {"stops_per_seq_median": np.median(stops_per_seq),
                "seqs_with_stops_fract": float(len(list(filter(lambda x: x > 0, stops_per_seq)))) / float(len(self))}

    def build_tree(self,
                   tree_name="phylo_tree",
                   tmp_dir="phylo_tree_tmp",
                   method="FastME",
                   options=None,
                   fastme_exec_path=None,
                   **kwargs):

        if fastme_exec_path is None:
            fastme_exec_path = self.fastme_exec_path

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        if method.lower() == "fastme":
            if "-b" not in options and "--bootstrap" not in options and \
                    self._distance_matrix is not None and not kwargs.get("rebuild_dist_matr", False):
                return self.distance_matrix.build_tree(tree_name=tree_name, tmp_dir=tmp_dir, method=method,
                                                       options=options, fastme_exec_path=fastme_exec_path)
            if self.aln_type.lower() in self.nucl_types:
                if not list(filter(None, fullmatch_regexp_list("-d.*", options))) \
                       and not list(filter(None, fullmatch_regexp_list("--dna.*", options))):
                    options["-d"] = True
            else:
                if not list(filter(None, fullmatch_regexp_list("-p.*", options))) \
                       and not list(filter(None, fullmatch_regexp_list("--protein.*", options))):
                    options["-p"] = True
            aln_phylip_path = os.path.join(tmp_dir, self.aln_name+".phylip")
            self.dump_alignment(aln_path=aln_phylip_path, aln_format="phylip")
            tree_path = os.path.join(tmp_dir, "tree.nwk")
            phylo_tree = run_fastme(input_data=aln_phylip_path,
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


class BlastHandler(ConfBase):

    def __init__(self,
                 inst_dir=None,
                 tmp_dir=None,
                 config_path=None,
                 logger=None):

        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_blast_tmp"

        self.inst_dir = inst_dir
        if self.inst_dir is None:
            self.inst_dir = conf_constants.blast_inst_dir
        self.tmp_dir = tmp_dir
        self.logger = logger

        super(BlastHandler, self).__init__(config_path=config_path)

    def make_blastdb(self, in_fasta, dbtype, db_name=None):
        if db_name:
            makeblastdb_cmd = os.path.join(self.inst_dir, "makeblastdb") + " -in " + in_fasta + " -dbtype " + dbtype + \
                              " -out " + db_name
        else:
            makeblastdb_cmd = os.path.join(self.inst_dir, "makeblastdb") + " -in " + in_fasta + " -dbtype " + dbtype
        subprocess.call(makeblastdb_cmd, shell=True)

    def run_blast_search(self, blast_type, query, db, out, num_threads=1, outfmt=7, max_hsps=100, **kwargs):
        if num_threads > 1 and kwargs.get("split_input", True):
            if self.logger is not None:
                self.logger.info("splitting '%s' into %s parts" % (query, num_threads))
            else:
                print("INFO: splitting '%s' into %s parts" % (query, num_threads))
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)
            if isinstance(query, SeqsDict):
                query_dict = query
                query_path = os.path.join(self.tmp_dir, "query.fasta")
            else:
                query_dict = load_fasta_to_dict(fasta_path=query, dat_path=os.path.join(self.tmp_dir, ".%s.dat" % query))
                query_path = query
            query_chunk_size = len(query_dict) // num_threads + 1
            p_list = list()
            query_seqs = list(query_dict.keys())
            i = 0
            query_chunks_list = list()
            for i in range(num_threads):
                query_chunk_path = None
                query_chunk_path = os.path.join(self.tmp_dir,
                                                ("_%s" % i).join(os.path.splitext(os.path.basename(query_path))))
                query_chunks_list.append([query_chunk_path, query_chunk_path + ".bl"])
                if len(query_dict) == i*query_chunk_size:
                    del query_chunks_list[-1]
                    continue
                elif (i+1) * query_chunk_size > len(query_dict):
                    query_dict.get_sample(query_seqs[i * query_chunk_size:]).dump(seqs_path=query_chunks_list[-1][0])
                else:
                    query_dict.get_sample(query_seqs[i*query_chunk_size:
                                                     (i+1)*query_chunk_size]).dump(seqs_path=query_chunks_list[-1][0])
                p = mp.Process(target=self.run_blast_search,
                               args=(blast_type, query_chunks_list[-1][0], db,
                                     query_chunks_list[-1][1], 1, outfmt, max_hsps),
                               kwargs=kwargs)
                p.start()
                p_list.append(p)
            for p in p_list:
                p.join()
            join_files(in_files_list=list(map(lambda p: p[1], query_chunks_list)), out_file_path=out)

            if kwargs.get("remove_tmp", True):
                shutil.rmtree(self.tmp_dir)
        else:
            if isinstance(query, SeqsDict):
                query_path = os.path.join(self.tmp_dir, "query.fasta")
            else:
                query_path = query
            blast_search_cmd = os.path.join(self.inst_dir, blast_type) + \
                               " -query " + query_path + \
                               " -db " + db + \
                               " -out " + out + \
                               " -word_size " + kwargs.get("word_size", str(3)) + \
                               " -num_threads " + str(num_threads) + \
                               " -outfmt " + str(outfmt) + \
                               " -max_hsps " + str(max_hsps)
            if self.logger is not None:
                self.logger.info("run '%s' command" % blast_search_cmd)
            else:
                print("INFO: run '%s' command" % blast_search_cmd)
            subprocess.call(blast_search_cmd, shell=True)


class HmmerHandler(ConfBase):

    def __init__(self,
                 inst_dir=None,
                 tmp_dir=None,
                 config_path=None,
                 logger=None):

        if inst_dir is None:
            inst_dir = conf_constants.hmmer_inst_dir
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_hmm_tmp"

        self.tmp_dir = tmp_dir
        self.inst_dir = inst_dir
        self.logger = logger

        super(HmmerHandler, self).__init__(config_path=config_path)

    def build_hmm_profile(self, profile_path, in_aln_path):
        hmmbuild_cmd = os.path.join(self.inst_dir, "hmmbuild") + " " + profile_path + " " + in_aln_path
        subprocess.call(hmmbuild_cmd, shell=True)
        return profile_path

    def make_profiles_db(self, profiles_list, profiles_db_path):
        join_files(profiles_list, profiles_db_path)
        hmmpress_cmd = os.path.join(self.inst_dir, "hmmpress") + " " + profiles_db_path
        subprocess.call(hmmpress_cmd, shell=True)

    def run_hmmsearch(self, in_profile_path, seqdb, out_path=None, shred_seqdb=False):
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if shred_seqdb:
            seqdb_path = shred_fasta(in_fasta=seqdb,
                                     shredded_fasta_path=os.path.join(self.tmp_dir, os.path.basename(seqdb)),
                                     part_l=50000,
                                     parts_ov=5000,)
        elif isinstance(seqdb, SeqsDict):
            seqdb_path = seqdb.dump(os.path.join(self.tmp_dir, "seqdb_to_search.fasta"))
        else:
            seqdb_path = seqdb
        if out_path is None:
            read_hsr = True
            out_path = os.path.splitext(in_profile_path) + ".hsr"
        else:
            read_hsr = False
        hmmsearch_cmd = os.path.join(self.inst_dir, "hmmsearch") + " " + in_profile_path + " " + seqdb_path
        subprocess.call(hmmsearch_cmd, shell=True)
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
        if read_hsr:
            return self.read_hsr(hsr_path=out_path, is_shred=shred_seqdb)
        else:
            return out_path

    def run_hmmscan(self, profiles_db, in_seqs, num_threads=4, out_path=None, shred_in_fasta=False):
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if shred_in_fasta:
            in_fasta_path = shred_fasta(in_fasta=in_seqs,
                                        shredded_fasta_path=os.path.join(self.tmp_dir, os.path.basename(in_seqs)),
                                        part_l=50000,
                                        parts_ov=5000)
        elif isinstance(in_seqs, SeqsDict):
            in_fasta_path = in_seqs.dump(os.path.join(self.tmp_dir, "seqfile_to_scan.fasta"))
        else:
            in_fasta_path = in_seqs
        if out_path is None:
            read_hsr = True
            if "." in in_seqs:
                out_path = ".".join(in_seqs.split(".")[:-1]) + ".hsr"
            else:
                out_path = in_seqs + ".hsr"
        else:
            read_hsr = False
        hmmscan_cmd = os.path.join(self.inst_dir, "hmmscan") + " --cpu " + str(num_threads) + " " + profiles_db + " " +\
                      in_fasta_path + " > " + out_path
        subprocess.call(hmmscan_cmd, shell=True)
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
        if read_hsr:
            return self.read_hsr(hsr_path=out_path, is_shred=shred_in_fasta)
        else:
            return out_path

    @staticmethod
    def read_hsr(hsr_path, is_shred=False):###
        return


def construct_mult_aln(seq_dict=None,
                       fasta_path=None,
                       method="MUSCLE",
                       aln_type=None,
                       muscle_exec_path=None,
                       mafft_exec_path=None,
                       msaprobs_exec_path=None,
                       emboss_inst_dir=None,
                       hmmer_inst_dir=None,
                       aln_name="mult_aln",
                       tmp_dir=None,
                       remove_tmp=True,
                       num_threads=None,
                       config_path=None,
                       logger=None,
                       **kwargs):

    if muscle_exec_path is None:
        muscle_exec_path = conf_constants.muscle_exec_path
    if mafft_exec_path is None:
        mafft_exec_path = conf_constants.mafft_exec_path
    if msaprobs_exec_path is None:
        msaprobs_exec_path = conf_constants.msaprobs_exec_path
    if emboss_inst_dir is None:
        emboss_inst_dir = conf_constants.emboss_inst_dir
    if hmmer_inst_dir is None:
        hmmer_inst_dir = conf_constants.hmmer_inst_dir
    if tmp_dir is None:
        tmp_dir = generate_random_string(10) + "_mult_aln_tmp"

    if not fasta_path and seq_dict:
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        fasta_path = os.path.join(tmp_dir, "seqs_to_aln.fasta")
        dump_fasta_dict(fasta_dict=seq_dict, fasta_path=fasta_path, replace_stops="X")  # maybe it is not good (* > X)
    if not fasta_path:
        if logger:
            logger.warning("No sequences input")
        else:
            print("No sequences input")
        return 1
    if not aln_type:
        detect_seqs_type(fasta_path)
    out_fasta_path = os.path.join(tmp_dir, aln_name+".fasta")

    if method.lower() == "muscle":
        if logger:
            logger.info("MUSCLE is starting")
        else:
            print("MUSCLE is starting")
        muscle_cmd = muscle_exec_path + " -in " + fasta_path + " -out " + out_fasta_path
        subprocess.call(muscle_cmd, shell=True)
        if logger:
            logger.info("MUSCLE finished")
        else:
            print("MUSCLE finished")

    if method.lower() == "mafft":
        if num_threads is None:
            num_threads = conf_constants.num_threads
        if logger:
            logger.info("MAFFT is starting")
        else:
            print("MAFFT is starting")
        mafft_cmd = mafft_exec_path + " --auto --op " + str(kwargs.get("op", kwargs.get("gap_open_penalty", 1.53))) + \
                    " --ep " + str(kwargs.get("ep", kwargs.get("gap_ext_penalty", 0.123))) + " --thread " + str(num_threads) + \
                    " " + fasta_path + " > " + out_fasta_path
        subprocess.call(mafft_cmd, shell=True)
        if logger:
            logger.info("MAFFT finished")
        else:
            print("MAFFT finished")

    if method.lower() == "msaprobs":
        if num_threads is None:
            num_threads = conf_constants.num_threads
        if logger:
            logger.info("MSAProbs is starting")
        else:
            print("MSAProbs is starting")
        msaprobs_cmd = msaprobs_exec_path + " -num_threads " + str(num_threads) + " -v " + \
                       fasta_path + " > " + out_fasta_path
        subprocess.call(msaprobs_cmd, shell=True)
        if logger:
            logger.info("MSAProbs finished")
        else:
            print("MSAProbs finished")

    mult_aln = MultAln.load_alignment(aln_fasta_path=out_fasta_path,
                                      aln_type=aln_type,
                                      aln_name=aln_name,
                                      config_path=config_path,
                                      logger=logger,
                                      restore_stops=True,
                                      **kwargs)
    mult_aln.emboss_inst_dir = emboss_inst_dir
    mult_aln.hmmer_inst_dir = hmmer_inst_dir
    mult_aln.tmp_dir = tmp_dir
    if remove_tmp:
        shutil.rmtree(tmp_dir)
    return mult_aln


def nucl_accord_prot(prot_seq, nucl_seq):
    nucl_seq_list = list()
    i = 0
    for aa in prot_seq:
        if aa == "-":
            nucl_seq_list.append("---")
        elif aa == ".":
            nucl_seq_list.append("...")
        else:
            nucl_seq_list.append(nucl_seq[i*3:(i+1)*3])
            i += 1
    return "".join(nucl_seq_list)


def search_profile(profile_dict, seqdb, seq_ids_to_orgs, work_dir="./", method="hmmer", hmmer_inst_dir=None,
                   shred_seqdb=False, config_path=None, logger=None, **kwargs):###
    p = SeqProfileInfo.load_from_dict(in_dict=profile_dict)
    if method.lower() == "hmmer":
        if kwargs.get("read_hsr_shred", False):
            hsr_path = os.path.join(work_dir, p.name + ".hsr")
        else:
            hsr_path = None
        hmmer_handler = HmmerHandler(inst_dir=hmmer_inst_dir,
                                     tmp_dir=os.path.join(work_dir, p.name+"_tmp"),
                                     config_path=config_path,
                                     logger=logger)
        search_res = hmmer_handler.run_hmmsearch(in_profile_path=p.path,
                                                 seqdb=seqdb,
                                                 out_path=hsr_path,
                                                 shred_seqdb=shred_seqdb)
        if hsr_path is not None:
            search_res = None
            search_res = hmmer_handler.read_hsr(hsr_path=hsr_path, is_shred=True)

        pass
    return


def get_kaks_gmean(kaks_array, ks_array=None, stype="negative", top_fract=None):
    if top_fract is None:
        top_fract = conf_constants.kaks_top_fract

    l = len(kaks_array)
    kaks_syn_array = sorted(filter(lambda p: True if p[0] > 0 else False,
                                   [(kaks_array[i], ks_array[i] if ks_array is not None else None) for i in range(l)]),
                            key=lambda p: p[0])
    kaks_list = list(map(lambda p: p[0], kaks_syn_array))
    ks_list = list(map(lambda p: p[1], kaks_syn_array))
    sep_ind = top_fract*float(l)

    ks = None
    if stype.lower() in ("n", "neg", "negative", "s", "stab", "stabilising"):
        kaks = gmean(kaks_list[: sep_ind])
        if ks_array is not None:
            ks = gmean(ks_list[: sep_ind])
    else:
        kaks = gmean(kaks_list[-sep_ind: ])
        if ks_array is not None:
            ks = gmean(ks_list[-sep_ind: ])
    return kaks, ks
