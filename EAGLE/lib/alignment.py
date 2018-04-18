import collections
import os
import shutil
import subprocess

from EAGLE.lib.general import load_fasta_to_dict, dump_fasta_dict, load_phylip_dist_matrix


class MultAln:

    def __init__(self, mult_aln_dict, aln_type=None, tmp_dir="tmp", logger=None):
        self.mult_aln_dict = mult_aln_dict
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

    def reduce_gaps(self, gap_percent=0.1, remove_seq=False):
        pass

    def get_distance_matrix(self, method="phylip"):
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if method.lower() == "phylip":
            aln_fasta_path = os.path.join(self.tmp_dir, "alignment.fasta")
            phylip_matrix_path = os.path.join(self.tmp_dir, "dist_matrix.phylip")
            dump_fasta_dict(fasta_dict=self.mult_aln_dict, fasta_path=aln_fasta_path)
            if self.aln_type.lower() in ("protein", "prot", "p"):
                phylip_cmd = "protdist -in " + aln_fasta_path + " -out " + phylip_matrix_path  # check cmd!
            else:
                phylip_cmd = "dnadist -in " + aln_fasta_path + " -out " + phylip_matrix_path  # check cmd!
            subprocess.call(phylip_cmd, shell=True)
            self.distance_matrix = load_phylip_dist_matrix(matrix_path=phylip_matrix_path)
        shutil.rmtree(self.tmp_dir)
        return  self.distance_matrix

    def remove_paralogs(self, seq_ids_to_orgs, only_dist=True):
        if not self.distance_matrix:
            self.get_distance_matrix()
        pass

    def to_gtf(self, gtf_path, fasta_path, meta_dict):
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
