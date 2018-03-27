import os
import shutil
import subprocess

from EAGLE.constants import EAGLE_logger
from EAGLE import conf_constants
from EAGLE.lib import read_fasta_to_dict


class MultAln:

    def __init__(self, mult_aln_dict, tmp_dir="tmp", hmmer_inst_dir=""):
        self.mult_aln_dict = mult_aln_dict
        self.tmp_dir = tmp_dir
        self.hmmer_inst_dir = hmmer_inst_dir

    def __getitem__(self, seq_id):
        return self.mult_aln_dict[seq_id]

    def __setitem__(self, seq_id, seq):
        self.mult_aln_dict[seq_id] = seq

    def get_distance_matrix(self):
        pass

    def remove_paralogs(self, seq_ids_to_orgs, only_dist=True):
        pass

    def to_gtf(self, gtf_path, fasta_path, meta_dict):
        pass

    def get_hmm_profile(self, method, profile_path):
        pass


def construct_mult_aln(seq_dict=None,
                       fasta_path=None,
                       method="MUSCLE",
                       aligner_inst_dir=conf_constants.muscle_inst_dir,
                       tmp_dir="tmp",
                       hmmer_inst_dir="",
                       remove_tmp=True):

    if not fasta_path and seq_dict:
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        fasta_path = os.path.join(tmp_dir, "seqs_to_aln.fasta")
    if not fasta_path:
        EAGLE_logger.warning("No sequences input")
        return 1
    out_fasta_path = os.path.join(tmp_dir, "mult_aln.fasta")

    if method.lower() == "muscle":
        muscle_cmd = os.path.join(aligner_inst_dir, "muscle") + " -in " + fasta_path + " -out " + out_fasta_path
        subprocess.call(muscle_cmd, shell=True)

    mult_aln_dict = read_fasta_to_dict(out_fasta_path)
    if remove_tmp:
        shutil.rmtree(tmp_dir)
    return MultAln(mult_aln_dict=mult_aln_dict, tmp_dir=tmp_dir, hmmer_inst_dir=hmmer_inst_dir)
