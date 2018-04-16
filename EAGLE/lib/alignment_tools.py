import os
import shutil
import subprocess
import collections

from EAGLE.constants import EAGLE_logger
from EAGLE import conf_constants
from EAGLE.lib import read_fasta_to_dict


class MultAln:

    def __init__(self, mult_aln_dict, aln_type=None, tmp_dir="tmp"):
        self.mult_aln_dict = mult_aln_dict
        self.aln_type = aln_type
        if not self.aln_type:
            self.aln_type = detect_aln_type(self.mult_aln_dict)
        self.distance_matrix = None
        self.tmp_dir = tmp_dir

    def __getitem__(self, seq_id):
        return self.mult_aln_dict[seq_id]

    def __setitem__(self, seq_id, seq):
        self.mult_aln_dict[seq_id] = seq

    def reduce_gaps(self, gap_percent=0.1, remove_seq=False):
        pass

    def get_distance_matrix(self):
        pass

    def remove_paralogs(self, seq_ids_to_orgs, only_dist=True):
        if not self.distance_matrix:
            self.get_distance_matrix()
        pass

    def to_gtf(self, gtf_path, fasta_path, meta_dict):
        pass

    def get_hmm_profile(self, method, profile_path):
        pass


def construct_mult_aln(seq_dict=None,
                       fasta_path=None,
                       method="MUSCLE",
                       aln_type=None,
                       aligner_inst_dir=conf_constants.muscle_inst_dir,
                       tmp_dir="tmp",
                       hmmer_inst_dir=conf_constants.hmmer_inst_dir,
                       remove_tmp=True):

    if not fasta_path and seq_dict:
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        fasta_path = os.path.join(tmp_dir, "seqs_to_aln.fasta")
    if not fasta_path:
        EAGLE_logger.warning("No sequences input")
        return 1
    if not aln_type:
        detect_aln_type(fasta_path)
    out_fasta_path = os.path.join(tmp_dir, "mult_aln.fasta")

    if method.lower() == "muscle":
        muscle_cmd = os.path.join(aligner_inst_dir, "muscle") + " -in " + fasta_path + " -out " + out_fasta_path
        subprocess.call(muscle_cmd, shell=True)

    mult_aln_dict = read_fasta_to_dict(out_fasta_path)
    if remove_tmp:
        shutil.rmtree(tmp_dir)
    return MultAln(mult_aln_dict=mult_aln_dict, aln_type=aln_type, tmp_dir=tmp_dir)


def detect_aln_type(fasta_path=None, fasta_dict=None, nuc_freq_thr=0.75):
    seqs_list = list()
    summ_l = 0
    if not fasta_dict and fasta_path:
        fasta_dict = read_fasta_to_dict(fasta_path).lower().replace("-", "")
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
