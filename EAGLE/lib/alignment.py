import collections
import os
import shutil
import subprocess
from collections import defaultdict

from EAGLE.lib.general import load_fasta_to_dict, dump_fasta_dict, load_phylip_dist_matrix, reduce_seq_names


class MultAln:

    def __init__(self, mult_aln_dict=None, aln_type=None, states_seq=None, tmp_dir="tmp", logger=None):
        if mult_aln_dict:
            self.mult_aln_dict, self.full_seq_names = reduce_seq_names(fasta_dict=mult_aln_dict, num_letters=10)
        else:
            self.mult_aln_dict = None
            self.full_seq_names = None
        self.states_seq = states_seq
        self.aln_type = aln_type
        if not self.aln_type:
            self.aln_type = detect_seqs_type(self.mult_aln_dict)
        self.distance_matrix = None
        self.tmp_dir = tmp_dir
        self.logger = logger

    def __getitem__(self, seq_id):
        return self.mult_aln_dict[seq_id]

    def __setitem__(self, seq_id, seq):
        self.mult_aln_dict[seq_id] = seq

    def dump_alignment(self, aln_f_path):
        pass

    def improve_aln(self,
                    max_gap_fract=0.95,  # maximal fraction of gaps in a column to keep it in alignment
                    max_mismatch_fract=1.0,  # maximal fraction of mismatches in a column to keep it in alignment
                    only_ends=True,
                    split_into_blocks=False,
                    remove_seq=False,  # if True returns a list of alignments with different strictness of sequences removing (different thresholds for clustering)
                    output='alignment',
                    inplace=False):

        pass

    def get_distance_matrix(self, method="phylip"):
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if method.lower() == "phylip":
            aln_fasta_path = os.path.join(self.tmp_dir, "alignment.fasta")
            phylip_matrix_path = os.path.join(self.tmp_dir, "dist_matrix.phylip")
            dump_fasta_dict(fasta_dict=self.mult_aln_dict, fasta_path=aln_fasta_path)
            if self.aln_type.lower() in ("protein", "prot", "p"):
                phylip_cmd = "fprotdist -sequence " + aln_fasta_path + " -outfile " + phylip_matrix_path
            else:
                phylip_cmd = "fdnadist -sequence " + aln_fasta_path + " -method f -outfile " + phylip_matrix_path
            subprocess.call(phylip_cmd, shell=True)
            self.distance_matrix = load_phylip_dist_matrix(matrix_path=phylip_matrix_path)
        shutil.rmtree(self.tmp_dir)
        return self.distance_matrix

    def remove_paralogs(self, seq_ids_to_orgs, method="min_dist", inplace=False):
        short_seq_ids_to_org = self._check_seq_ids_to_org(seq_ids_to_orgs)
        if not self.distance_matrix:
            self.get_distance_matrix()
        org_dist_dict = defaultdict(dict)
        seq_ids = self.distance_matrix.index
        for seq_id in seq_ids:
            org_dist_dict[seq_ids_to_orgs[seq_id]][seq_id] = sum(map(float, list(self.distance_matrix.loc(seq_id))))
        for org in org_dist_dict.keys():
            org_dist_dict[org] = sorted(org_dist_dict[org].items(), key=lambda x: x[1])
        if method.lower() in ("minimal_distance", "min_dist", "md"):
            full_seq_names_filt = dict(filter(lambda x: x[0] == org_dist_dict.get(short_seq_ids_to_org.get(x[0], None),
                                                                                  ((None, None), None))[0][0],
                                              self.full_seq_names.items()))
            mult_aln_dict_filt = dict(filter(lambda x: x[0] in full_seq_names_filt.keys(), self.mult_aln_dict.items()))
        else:
            # TODO: write spec_pos method
            full_seq_names_filt = None
            mult_aln_dict_filt = None
        if inplace:
            self.mult_aln_dict = mult_aln_dict_filt
            self.full_seq_names = full_seq_names_filt
        else:
            filtered_aln = self
            filtered_aln.mult_aln_dict = mult_aln_dict_filt
            filtered_aln.full_seq_names = full_seq_names_filt
            return filtered_aln

    def _check_seq_ids_to_org(self, seq_ids_to_orgs):
        to_short_dict = dict()
        long_suc_num = 0
        short_suc_num = 0
        for key in self.full_seq_names.keys:
            if seq_ids_to_orgs.get(key, None):
                short_suc_num += 1
            if seq_ids_to_orgs.get(self.full_seq_names[key], None):
                to_short_dict[self.full_seq_names[key]] = key
                long_suc_num += 1
        if short_suc_num >= long_suc_num:
            return seq_ids_to_orgs
        else:
            return dict(map(lambda x: (to_short_dict.get(x[0], None),  x[1]), seq_ids_to_orgs.items()))

    def get_blocks_tsv(self, gtf_path, fasta_path, meta_dict):
        cut_ends_dict = self.improve_aln(output='coords')

        pass

    def get_hmm_profile(self, method, profile_path):
        pass

    def nucl_by_prot_aln(self, nucl_fasta_dict=None, nucl_fasta_path=None):
        if self.aln_type.lower() not in ("protein", "prot", "p"):
            if self.logger:
                self.logger.warning("Reference alignment type is not protein")
            else:
                print "Reference alignment type is not protein"
            return 1
        if not nucl_fasta_dict:
            if nucl_fasta_path:
                nucl_fasta_dict = load_fasta_to_dict(fasta_path=nucl_fasta_path)
            else:
                if self.logger:
                    self.logger("No nucleotide sequences are input")
                else:
                    print "No nucleotide sequences are input"
                return 1
        pass


def construct_mult_aln(seq_dict=None,
                       fasta_path=None,
                       method="MUSCLE",
                       aln_type=None,
                       aligner_inst_dir="",
                       tmp_dir="tmp",
                       hmmer_inst_dir="",
                       remove_tmp=True,
                       logger=None):

    if not fasta_path and seq_dict:
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        fasta_path = os.path.join(tmp_dir, "seqs_to_aln.fasta")
    if not fasta_path:
        if logger:
            logger.warning("No sequences input")
        else:
            print "No sequences input"
        return 1
    if not aln_type:
        detect_seqs_type(fasta_path)
    out_fasta_path = os.path.join(tmp_dir, "mult_aln.fasta")

    if method.lower() == "muscle":
        muscle_cmd = os.path.join(aligner_inst_dir, "muscle") + " -in " + fasta_path + " -out " + out_fasta_path
        subprocess.call(muscle_cmd, shell=True)

    mult_aln_dict = load_fasta_to_dict(out_fasta_path)
    if remove_tmp:
        shutil.rmtree(tmp_dir)
    return MultAln(mult_aln_dict=mult_aln_dict, aln_type=aln_type, tmp_dir=tmp_dir, logger=logger)


def detect_seqs_type(fasta_path=None, fasta_dict=None, nuc_freq_thr=0.75):
    seqs_list = list()
    summ_l = 0
    if not fasta_dict and fasta_path:
        fasta_dict = load_fasta_to_dict(fasta_path).lower().replace("-", "")
        for seq_key in fasta_dict.keys():
            seqs_list.append(fasta_dict[seq_key])
            summ_l += len(fasta_dict[seq_key])
        let_counts = collections.Counter("".join(seqs_list))
        if float(let_counts.get("a", 0)+let_counts.get("c", 0)+let_counts.get("g", 0)+let_counts.get("t", 0))/float(summ_l) >= nuc_freq_thr:
            return "nucl"
        else:
            return "prot"
    else:
        return None
