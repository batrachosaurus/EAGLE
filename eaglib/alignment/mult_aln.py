import os
import re
import shutil
import subprocess
from copy import deepcopy
from collections import defaultdict, Counter

from deprecated import deprecated
import numpy as np
import pandas as pd
from scipy.stats import gmean
from Bio.Seq import Seq

from eagle.constants import conf_constants
from eaglib._utils.strings import generate_random_string, fullmatch_regexp_list
from eaglib._utils.logging import send_log_message
from eaglib.seqs import SeqsDict, load_fasta_to_dict, nucl_accord_prot
from eaglib.phylo import DistanceMatrix, run_fastme


class MultAln(SeqsDict):

    name0 = "mult_aln"

    dist_matr_method = "FastME"
    dist_matr_options = dict()

    def __init__(self, seqs_order, seqs_array,
                 seq_info_dict=None,
                 seqs_type=None,
                 low_memory=False,
                 states_seq=None,
                 name=name0,
                 tmp_dir=None,
                 emboss_inst_dir=None,
                 hmmer_inst_dir=None,
                 infernal_inst_dir=None,
                 fastme_exec_path=None,
                 kaks_calculator_exec_path=None,
                 **kwargs):

        super(MultAln, self).__init__(seqs_order, seqs_array,
                                      seq_info_dict=seq_info_dict,
                                      seqs_type=seqs_type,
                                      low_memory=low_memory,
                                      **kwargs)

        if states_seq is None:
            states_seq = list()
        if emboss_inst_dir is None:
            emboss_inst_dir = conf_constants.emboss_inst_dir
        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir
        if fastme_exec_path is None:
            fastme_exec_path = conf_constants.fastme_exec_path
        if kaks_calculator_exec_path is None:
            kaks_calculator_exec_path = conf_constants.kaks_calculator_exec_path

        self.states_seq = states_seq
        self._distance_matrix = None
        self.name = name
        if tmp_dir is None:
            tmp_dir = self.name + "_%s_tmp" % generate_random_string(10)
        self.tmp_dir = tmp_dir
        self.emboss_inst_dir = emboss_inst_dir
        self.infernal_inst_dir = infernal_inst_dir
        self.hmmer_inst_dir = hmmer_inst_dir
        self.fastme_exec_path = fastme_exec_path
        self.kaks_calculator_exec_path = kaks_calculator_exec_path

    def _init_subset(self, seqs_order, seqs_array,
                     seq_info_dict=None,
                     seqs_type=None,
                     low_memory=None,
                     states_seq=None,
                     name=name0,
                     tmp_dir=None,
                     emboss_inst_dir=None,
                     hmmer_inst_dir=None,
                     infernal_inst_dir=None,
                     fastme_exec_path=None,
                     kaks_calculator_exec_path=None,
                     **kwargs):

        if states_seq is None:
            states_seq = deepcopy(self.states_seq)
        if name is None:
            name = self.name
        if tmp_dir is None:
            tmp_dir = self.tmp_dir
        if emboss_inst_dir is None:
            emboss_inst_dir = self.emboss_inst_dir
        if hmmer_inst_dir is None:
            hmmer_inst_dir = self.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = self.infernal_inst_dir
        if fastme_exec_path is None:
            fastme_exec_path = self.fastme_exec_path
        if kaks_calculator_exec_path is None:
            kaks_calculator_exec_path = self.kaks_calculator_exec_path

        mult_aln = super(MultAln, self)._init_subset(seqs_order, seqs_array,
                                                     seq_info_dict=seq_info_dict,
                                                     seqs_type=seqs_type,
                                                     low_memory=low_memory,
                                                     states_seq=states_seq,
                                                     name=name,
                                                     tmp_dir=tmp_dir,
                                                     emboss_inst_dir=emboss_inst_dir,
                                                     hmmer_inst_dir=hmmer_inst_dir,
                                                     infernal_inst_dir=infernal_inst_dir,
                                                     fastme_exec_path=fastme_exec_path,
                                                     kaks_calculator_exec_path=kaks_calculator_exec_path,
                                                     **kwargs)
        return mult_aln

    def _inplace(self, mult_aln):
        assert isinstance(mult_aln, MultAln), "ERROR: value for argument 'mult_aln' should be an instance of " \
                                              "eagle.eaglib.alignment.MultAln class"
        self.seqs_order = mult_aln.seqs_order
        self.seqs_array = mult_aln.seqs_array
        self._distance_matrix = None
        self.seq_info_dict = mult_aln.seq_info_dict
        self.nucl_seqs_dict = mult_aln.nucl_seqs_dict
        self.full_to_short_seq_names = mult_aln.full_to_short_seq_names

    def col(self, item):
        seqs_dict = self.load_from_dict({seq_name: seq[item] for seq_name, seq in self.items()},
                                        seqs_type=self.seqs_type)
        return self._init_subset(seqs_dict.seqs_order, seqs_dict.seqs_array)

    def __setitem__(self, key, value):
        super(MultAln, self).__setitem__(key=key, value=value)
        self._distance_matrix = None

    def __delitem__(self, key):
        self._distance_matrix = None
        super(MultAln, self).__delitem__(key=key)

    @property
    @deprecated(reason="use MultAln instance directly")
    def mult_aln_dict(self):
        return self

    @property
    @deprecated(reason="use 'name' attribute")
    def aln_name(self):
        return self.name

    @aln_name.setter
    @deprecated(reason="use 'name' attribute")
    def aln_name(self, name):
        self.name = name

    @property
    @deprecated(reason="use 'seqs_type' attribute")
    def aln_type(self):
        return self.seqs_type

    @aln_type.setter
    @deprecated(reason="use 'seqs_type' attribute")
    def aln_type(self, aln_type):
        self.seqs_type = aln_type

    def rename_seqs(self, old_to_new_dict):
        super(MultAln, self).rename_seqs(old_to_new_dict=old_to_new_dict)
        self._distance_matrix = None

    @property
    def length(self):
        if self.seqs_order:
            return len(self[self.seq_names[0]])
        else:
            return 0

    @property
    def shape(self):
        return self.num_seqs, self.length

    @property
    @deprecated(reason="use short_id attribute")
    def mult_aln_dict_short_id(self):
        return self.short_id

    @deprecated(reason="use 'dump' method")
    def dump_alignment(self, aln_path, aln_format="fasta"):
        return self.dump(aln_path, format=aln_format)

    def dump(self, fname=None, format="fasta", overwrite=True, **kwargs):
        super(MultAln, self).dump(fname, format=format, overwrite=overwrite, shape=self.shape, **kwargs)

        return self.short_to_full_seq_names  # TODO: update output (remained for backward compatibility)

    @classmethod
    @deprecated(reason="use 'load_from_file' method")
    def load_alignment(cls, aln_path=None, aln_format="fasta", aln_type=None, aln_name=None, logger=None, **kwargs):
        if "aln_fasta_path" in kwargs:
            aln_path = None
            aln_path = kwargs["aln_fasta_path"]

        return cls.load_from_file(fname=aln_path, format=aln_format, seqs_type=aln_type, name=aln_name, logger=logger,
                                  **kwargs)

    @classmethod
    def load_from_file(cls, fname=None, format="fasta", seqs_type=None, low_memory=None, **kwargs):
        if format.lower() in ("fasta", "fas", "fa"):
            with open(fname) as fasta_f:
                seq_l = 0
                for line_ in fasta_f:
                    line = None
                    line = line_.strip()
                    if not line:
                        continue
                    if line[0] == ">" and seq_l > 0:
                        break
                    else:
                        seq_l += len(line)
        if format.lower() in ("phylip", "phy", "ph"):
            print("not implemented yet")  # TODO: implement phylip format reading

        chunk_size = kwargs.pop("chunk_size", seq_l)
        return super(MultAln, cls).load_from_file(fname=fname, format=format, seqs_type=seqs_type,
                                                  low_memory=low_memory, chunk_size=chunk_size, **kwargs)

    @classmethod
    def load_from_dict(cls, in_dict, seqs_type=None, low_memory=None, **kwargs):
        chunk_size = kwargs.pop("chunk_size", len(in_dict[list(in_dict.keys())[0]]))
        return super(MultAln, cls).load_from_dict(in_dict=in_dict, seqs_type=seqs_type, low_memory=low_memory,
                                                  chunk_size=chunk_size, **kwargs)

    def improve_aln(self,
                    max_gap_fract=0.75,
                    max_mismatch_fract=1.0,
                    remove_seq=False,
                    dist_filt=False,
                    num_threads=None,
                    output='alignment',
                    inplace=False,
                    **kwargs):

        """
        Method for removing not relevant parts of alignment
        :param max_gap_fract: maximal fraction of gaps in a column to keep it in alignment
        :param max_mismatch_fract: 1.0 - 'minimal allowed frequency of most frequent letter in a column'
        :param remove_seq: if True can remove sequences for cleaning gaps between aln blocks
        :param dist_filt: if True returns a list of alignments (or blocks coordinates) with different strictness of
                          sequences removing (different thresholds for clustering)
        :param output:
        :param inplace:
        :param kwargs:
        :return:
        """
        # TODO: write methods for removing gaps between blocks
        send_log_message("'%s' multiple alignment improvement started" % self.name, mes_type='info')
        to_return = None
        if num_threads is None:
            num_threads = num_threads

        if dist_filt:
            # Constructs a list of rarefied by different distance thresholds alignments and runs improve_aln on it with dist_filt=False
            pass

        # TODO: speed up that with threading
        block_coords = list()
        for i in range(self.length):
            if self.states_seq:
                if self.states_seq[i] in ("match", "specific"):
                    block_coords = self._update_block_coords(i, block_coords)
                    continue
            let_counts = Counter("".join(self[seq_id][i].lower() for seq_id in self.seq_names))
            if float(let_counts["-"]) <= max_gap_fract * float(self.num_seqs):
                del let_counts["-"]
                if sorted(let_counts.values(), reverse=True)[0] >= (1.0 - max_mismatch_fract) * float(self.num_seqs):
                    block_coords = self._update_block_coords(i, block_coords)
            elif remove_seq:
                # TODO: write detection of seqs to remove
                pass

        if output.lower() in ("coordinates", "coords", "coord", "c"):
            to_return = block_coords
        else:
            if block_coords:
                if remove_seq:
                    impr_mult_aln_dict = self.load_from_dict(
                        {seq_id: "".join(self[seq_id][c1:c2] for c1, c2 in block_coords) for seq_id in self.seq_names}  # TODO: self.seq_names?
                    )
                else:
                    impr_mult_aln_dict = self.load_from_dict(
                        {seq_id: self[seq_id][block_coords[0][0]:block_coords[-1][1]] for seq_id in self.seq_names}
                    )
                improved_mult_aln = self._init_subset(impr_mult_aln_dict.seqs_order, impr_mult_aln_dict.seqs_array,
                                                      name="improved_" + self.name)
                if inplace:
                    self._inplace(improved_mult_aln)
                else:
                    to_return = improved_mult_aln
            else:
                send_log_message("no columns left in '%s' multiple alignment" % self.name, mes_type='warning')

        send_log_message("'%s' multiple alignment improved" % self.name, mes_type='info')
        return to_return

    @staticmethod
    def _update_block_coords(i, block_coords):
        if not block_coords or block_coords[-1][1] < i:
            block_coords.append(np.array([i, i + 1]))
        else:
            block_coords[-1][1] = i + 1
        return block_coords

    @staticmethod
    def _get_coords_for_seq(coords, gapped_seq):
        no_gaps_coords = list()
        for coord_pair in coords:
            c1_shift = gapped_seq[:coord_pair[0]].count("-")
            c2_shift = gapped_seq[:coord_pair[1]].count("-")
            no_gaps_coords.append((coord_pair[0] - c1_shift, coord_pair[1] - c2_shift))
        return no_gaps_coords

    @property
    def distance_matrix(self):
        if self._distance_matrix is None:
            self._distance_matrix = DistanceMatrix.calculate(mult_aln=self,
                                                             method=self.dist_matr_method,
                                                             options=self.dist_matr_options,
                                                             emboss_inst_dir=self.emboss_inst_dir,
                                                             fastme_exec_path=self.fastme_exec_path)
        return self._distance_matrix

    @deprecated(reason="Use 'calculate_distance_matrix' method")
    def get_distance_matrix(self, method=None, options=None):
        return self.calculate_distance_matrix(method=method, options=options)

    def calculate_distance_matrix(self, method=None, options=None, **kwargs):
        if method is None:
            method = self.dist_matr_method
        if options is None:
            options = self.dist_matr_options

        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        if method.lower() == "phylip":
            aln_fasta_path = os.path.join(self.tmp_dir, self.name+".fasta")
            phylip_matrix_path = os.path.join(self.tmp_dir, self.name+".phylip")
            if not self:
                send_log_message("no sequences in alignment", mes_type="warning")
                return
            self.dump(aln_fasta_path, format="fasta")
            if self.seqs_type.lower() in self.prot_types:
                send_log_message("protdist is starting", mes_type="info")
                phylip_cmd = os.path.join(self.emboss_inst_dir, "fprotdist") + " -sequence " + aln_fasta_path + \
                             " -method d -outfile " + phylip_matrix_path
            else:
                send_log_message("dnadist is starting", mes_type="info")
                phylip_cmd = os.path.join(self.emboss_inst_dir, "fdnadist") + " -sequence " + aln_fasta_path + \
                             " -method f -outfile " + phylip_matrix_path
            subprocess.call(phylip_cmd, shell=True)
            send_log_message("distance calculations finished", mes_type="info")
            self._distance_matrix = DistanceMatrix.load(matrix_path=phylip_matrix_path,
                                                        matr_format="phylip",
                                                        short_to_full_seq_names=self.short_to_full_seq_names,
                                                        emboss_inst_dir=self.emboss_inst_dir,
                                                        fastme_exec_path=self.fastme_exec_path,
                                                        **kwargs)

        if method.lower() == "fastme":
            aln_phylip_path = os.path.join(self.tmp_dir, self.name+".phylip")
            phylip_matrix_path = os.path.join(self.tmp_dir, self.name+"_dm.phylip")
            if not self:
                send_log_message(message="no sequences in alignment", mes_type="e")
                return
            self.dump(aln_phylip_path, format="phylip")
            for option in ("-O", "--output_matrix"):
                options.pop(option, None)
            if self.seqs_type.lower() in self.prot_types:
                fastme_options = {"--protein": "L", "-O": phylip_matrix_path}
            else:
                fastme_options = {"--dna": 4, "-O": phylip_matrix_path}
            fastme_options.update(options)
            send_log_message(message="distances calculation started", mes_type="i")
            run_fastme(input_data=aln_phylip_path, options=fastme_options, fastme_exec_path=self.fastme_exec_path)
            send_log_message(message="distances calculation finished", mes_type="i")
            self._distance_matrix = DistanceMatrix.load(matrix_path=phylip_matrix_path,
                                                        matr_format="phylip",
                                                        short_to_full_seq_names=self.short_to_full_seq_names,
                                                        emboss_inst_dir=self.emboss_inst_dir,
                                                        fastme_exec_path=self.fastme_exec_path,
                                                        **kwargs)
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
        self._distance_matrix.aln_name = self.name
        self._distance_matrix.aln_type = self.seqs_type

        return self._distance_matrix

    def remove_paralogs(self, seq_ids_to_orgs=None, method="min_dist", inplace=False):
        """

        :param seq_ids_to_orgs: {seq_id: org_name}
        :param method:
        :param inplace:
        :return:
        """
        send_log_message(message="paralogs removing started", mes_type='i')

        if seq_ids_to_orgs is not None:
            self.seq_ids_to_orgs = seq_ids_to_orgs

        org_dists = defaultdict(dict)
        for seq_name in self.distance_matrix.seq_names:
            org_dists[self.seq_ids_to_orgs[seq_name]][seq_name] = sum(map(float, self.distance_matrix[seq_name]))
        sorted_org_dists = dict()
        for org in org_dists:
            sorted_org_dists[org] = sorted(org_dists[org].items(), key=lambda x: x[1])
        mult_aln_dict_filt = dict()
        if method.lower() in ("minimal_distance", "min_dist", "md"):
            for seq_name in self.seq_names:
                if seq_name == sorted_org_dists[self.seq_ids_to_orgs[seq_name]][0][0]:
                    mult_aln_dict_filt[seq_name] = self.mult_aln_dict[seq_name]
        else:
            # TODO: write spec_pos method
            pass

        filtered_aln_seqs = self.load_from_dict(in_dict=mult_aln_dict_filt)
        filtered_aln = self._init_subset(filtered_aln_seqs.seqs_order, filtered_aln_seqs.seqs_array,
                                         name="no_par_" + self.name)

        if inplace:
            self._inplace(filtered_aln)
        send_log_message("paralogs removed", "info")
        if not inplace:
            return filtered_aln

    def build_profile(self, profile_name=None, profile_path=None, method="hmmer", **kwargs):
        from eaglib.alignment.seq_profiles import SeqsProfile

        return SeqsProfile.build(mult_aln=self, method=method, name=profile_name, path=profile_path, **kwargs)

    @deprecated(reason="use method build_profile")
    def get_hmm_profile(self, profile_path, method="hmmer"):
        return self.build_profile(profile_path=profile_path, method=method)

    def nucl_by_prot_aln(self, nucl_fasta_dict=None, nucl_fasta_path=None):
        if self.aln_type.lower() not in self.prot_types:
            send_log_message(message="reference alignment type is not protein", mes_type='e')
            return

        if nucl_fasta_dict is not None:
            self.nucl_seqs_dict = nucl_fasta_dict
        elif nucl_fasta_path is not None:
            self.nucl_seqs_dict = load_fasta_to_dict(fasta_path=nucl_fasta_path)
        if not self.nucl_seqs_dict:
            send_log_message(message="no nucleotide sequences are input", mes_type="e")
            return

        nucl_aln_dict = dict()
        for seq_name in self.seq_names:
            to_search = None
            if "-" in self[seq_name] or "." in self[seq_name]:
                to_search = re.sub("[-\.]", "", self[seq_name]).replace("*", ".")
            else:
                to_search = self[seq_name].replace("*", ".")
            match = re.search(to_search, str(Seq(self.nucl_seqs_dict[seq_name]).translate()))
            nucl_aln_dict[seq_name] = nucl_accord_prot(self[seq_name],
                                                       self.nucl_seqs_dict[seq_name][
                                                       match.start() * 3: match.end() * 3])
        nucl_aln_seqs = self.load_from_dict(nucl_aln_dict)
        nucl_aln = self._init_subset(nucl_aln_seqs.seqs_order, nucl_aln_seqs.seqs_array,
                                     seqs_type=self.nucl_type,
                                     name="nucl_" + self.name)
        if self.states_seq:
            nucl_aln.states_seq = list()  # TODO: transform to codons states sequence
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
            windows_list.append(self.col[i: i + window_l])
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
            for seq_id in self:
                s_num_dict[self[seq_id][i].lower()] += 1
            all_s_num = sum(s_num_dict.values())
            if float(s_num_dict.get("-", 0)) / float(all_s_num) <= 1.0 - cons_thr:
                if float(sorted(s_num_dict.values(), reverse=True)[0]) / float(all_s_num) >= cons_thr:
                    cln += 1
        return cln

    def rarefy(self, seqs_to_remain=100, **kwargs):
        return super(MultAln, self).rarefy(seqs_to_remain=seqs_to_remain, name="raref_" + self.name, **kwargs)

    def calculate_KaKs_windows(self,
                               nucl_seqs_dict=None,
                               window_l=None,
                               method="KaKs_Calculator",
                               only_first_seq=False,
                               **kwargs):
        if nucl_seqs_dict is not None:
            self.nucl_seqs_dict = nucl_seqs_dict
        if window_l is None:
            window_l = conf_constants.kaks_window_l

        if self.seqs_type.lower() in self.prot_types:
            if not self.nucl_seqs_dict:
                send_log_message(message="protein alignment but no value input for argument 'nucl_seqs_dict'",
                                 mes_type='e')
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
            w.storage_dir = self.tmp_dir
            w_kaks = w.calculate_KaKs(method=method,
                                      remove_tmp=False,
                                      only_first_seq=only_first_seq,
                                      raref_base=kwargs.pop("raref_base", 10.0),
                                      **kwargs)
            if type(w_kaks) is float or isinstance(w_kaks, np.floating):
                if w_kaks["Ka/Ks"] >= 0.0 and not pd.isna(w_kaks["Ka/Ks"]):
                    kaks_list.append(w_kaks["Ka/Ks"])
                    ks_list.append(w_kaks["Ks"])
                else:
                    kaks_list.append(-1.0)
                    ks_list.append(-1.0)

        shutil.rmtree(self.tmp_dir, ignore_errors=True)
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

        if self.seqs_type.lower() in self.prot_types:
            if not self.nucl_seqs_dict:
                send_log_message(message="protein alignment but no value input for argument 'nucl_seqs_dict'",
                                 mes_type="e")
                return
            nucl_aln = self.nucl_by_prot_aln()
            return nucl_aln.calculate_KaKs(method=method, max_kaks=max_kaks, only_first_seq=only_first_seq, **kwargs)

        seq_pairs = self.generate_seqs_pairs(only_first_seq=only_first_seq, raref_base=kwargs.get("raref_base", None))
        if method.lower() == "kaks_calculator":
            pairs_axt_path = os.path.join(self.tmp_dir, self.name + ".axt")
            with open(pairs_axt_path, "w") as pairs_axt_f:
                for seq_pair in seq_pairs:
                    pairs_axt_f.write("vs".join(seq_pair) + "\n")
                    pairs_axt_f.write(seq_pairs[seq_pair][0] + "\n")
                    pairs_axt_f.write(seq_pairs[seq_pair][1] + "\n\n")
            out_path = os.path.join(self.tmp_dir, self.name + ".kaks")
            kaks_calculator_cmd = self.kaks_calculator_exec_path + \
                                  " -i " + pairs_axt_path + " -o " + out_path + " -m YN"
            send_log_message(message="running command: '%s'" % kaks_calculator_cmd, mes_type='i')
            subprocess.call(kaks_calculator_cmd, shell=True)
            try:
                kaks_df = pd.read_csv(out_path, sep="\t")
                kaks = gmean(kaks_df["Ka/Ks"].apply(lambda v: min(v, max_kaks)))
                ks = gmean(kaks_df["Ks"])
            except:
                kaks = -1.0
                ks = -1.0
        if kwargs.get("remove_tmp", True):
            shutil.rmtree(self.tmp_dir, ignore_errors=True)
        return {"Ka/Ks": kaks, "Ks": ks}

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
            if self.seqs_type.lower() in self.nucl_types:
                if not list(filter(None, fullmatch_regexp_list("-d.*", options))) \
                        and not list(filter(None, fullmatch_regexp_list("--dna.*", options))):
                    options["-d"] = True
            else:
                if not list(filter(None, fullmatch_regexp_list("-p.*", options))) \
                        and not list(filter(None, fullmatch_regexp_list("--protein.*", options))):
                    options["-p"] = True
            aln_phylip_path = os.path.join(tmp_dir, self.name + ".phylip")
            self.dump(aln_phylip_path, format="phylip")
            tree_path = os.path.join(tmp_dir, "tree.nwk")
            phylo_tree = run_fastme(input_data=aln_phylip_path,
                                    output_tree=tree_path,
                                    options=options,
                                    fastme_exec_path=fastme_exec_path)[0]
            phylo_tree.tree_name = tree_name
            phylo_tree.full_seq_names = self.short_to_full_seq_names
            phylo_tree.storage_dir = tmp_dir
            send_log_message(message="phylogenetic tree built with FastME", mes_type="i")
        else:
            return
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return phylo_tree

    def stop_codons_stats(self, improve_aln=False, **kwargs):
        if improve_aln:
            return self.improve_aln(inplace=False, **kwargs).stop_codons_stats(improve_aln=False, **kwargs)
        else:
            return super(MultAln, self).stop_codons_stats(improve_aln=False, **kwargs)


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


def construct_mult_aln(seqs_dict=None,
                       fasta_path=None,
                       method="MUSCLE",
                       seqs_type=None,
                       muscle_exec_path=None,
                       mafft_exec_path=None,
                       msaprobs_exec_path=None,
                       emboss_inst_dir=None,
                       hmmer_inst_dir=None,
                       infernal_inst_dir=None,
                       fastme_exec_path=None,
                       kaks_calculator_exec_path=None,
                       aln_name="mult_aln",
                       tmp_dir=None,
                       remove_tmp=True,
                       num_threads=None,
                       config_path=None,
                       **kwargs):

    if "seq_dict" in kwargs:
        seqs_dict = kwargs["seq_dict"]
    if "aln_type" in kwargs:
        seqs_type = kwargs["aln_type"]

    if seqs_dict is None and fasta_path is not None:
        seqs_dict = SeqsDict.load_from_file(fasta_path, format="fasta", seqs_type=seqs_type, **kwargs)

    if isinstance(seqs_dict, SeqsDict):
        return seqs_dict.construct_mult_aln(method=method,
                                            muscle_exec_path=muscle_exec_path,
                                            mafft_exec_path=mafft_exec_path,
                                            msaprobs_exec_path=msaprobs_exec_path,
                                            emboss_inst_dir=emboss_inst_dir,
                                            hmmer_inst_dir=hmmer_inst_dir,
                                            infernal_inst_dir=infernal_inst_dir,
                                            fastme_exec_path=fastme_exec_path,
                                            kaks_calculator_exec_path=kaks_calculator_exec_path,
                                            aln_name=aln_name,
                                            tmp_dir=tmp_dir,
                                            remove_tmp=remove_tmp,
                                            num_threads=num_threads,
                                            **kwargs)
