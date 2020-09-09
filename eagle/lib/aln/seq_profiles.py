import os
import shutil
import subprocess

from eagle.constants import conf_constants
from eagle.lib.general import generate_random_string, join_files
from eagle.lib.seqs import SeqsDict


class SeqProfiles(object):

    hmmer_key = "hmmer"
    infernal_key = "infernal"

    def __init__(self,
                 method=hmmer_key,
                 hmmer_inst_dir=None,
                 infernal_inst_dir=None,
                 seqdb=None,
                 profiles_db=None,
                 tmp_dir=None,
                 logger=None):

        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_profile_tmp"

        self.method = method
        self.hmmer_inst_dir = hmmer_inst_dir
        self.infernal_inst_dir = infernal_inst_dir
        self.seqdb = seqdb
        self.profiles_db = profiles_db
        self.logger = logger
        self.tmp_dir = tmp_dir

    def search(self):
        return

    def scan(self):
        return

    def build_profile(self):
        return

    def build_seqdb(self):
        return

    def build_profile_db(self):
        return


class HMMERHandler(object):

    def __init__(self, inst_dir=None, seqdb=None, profiles_db=None, tmp_dir=None, logger=None):
        if inst_dir is None:
            inst_dir = conf_constants.hmmer_inst_dir
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_hmmer_tmp"

        self.inst_dir = inst_dir
        self.seqdb = seqdb
        self.profiles_db = profiles_db
        self.logger = logger
        self.tmp_dir = tmp_dir

    def build_hmm_profile(self, profile_path, in_aln_path):###  in_aln - MultAln?
        hmmbuild_cmd = os.path.join(self.inst_dir, "hmmbuild") + " " + profile_path + " " + in_aln_path
        subprocess.call(hmmbuild_cmd, shell=True)
        return profile_path

    def make_profiles_db(self, profiles_list, profiles_db_path):
        join_files(profiles_list, profiles_db_path)### profiles list should contain ProfileInfo instancies
        hmmpress_cmd = os.path.join(self.inst_dir, "hmmpress") + " " + profiles_db_path
        subprocess.call(hmmpress_cmd, shell=True)

    def run_hmmsearch(self, in_profile_path, seqdb=None, out_path=None, cpu=1, shred_seqdb=False):
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
        hmmsearch_cmd = os.path.join(self.inst_dir, "hmmsearch") + \
                        " " + in_profile_path + \
                        " " + seqdb_path + \
                        " --cpu " + str(cpu)
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


class InfernalHandler(object):

    def __init__(self, inst_dir=None, seqdb=None, profiles_db=None, tmp_dir=None, logger=None):
        if inst_dir is None:
            inst_dir = conf_constants.infernal_inst_dir
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_infernal_tmp"

        self.inst_dir = inst_dir
        self.seqdb = seqdb
        self.profiles_db = profiles_db
        self.logger = logger
        self.tmp_dir = tmp_dir
