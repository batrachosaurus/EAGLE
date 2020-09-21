import os
import shutil
import subprocess
from collections import defaultdict

from Bio.SearchIO.HmmerIO.hmmer3_text import Hmmer3TextParser

from eagle.constants import conf_constants
from eagle.lib.general import generate_random_string, join_files
from eagle.lib.seqs import SeqsDict, shred_seqs
from eagle.lib.alignment.mult_aln import MultAln
from eagledb.scheme import SeqsProfileInfo


HMMER_KEY = "hmmer"
INFERNAL_KEY = "infernal"


class SeqsProfile(object):

    def __init__(self, seqs_profile_info: SeqsProfileInfo,
                 method=HMMER_KEY, hmmer_inst_dir=None, infernal_inst_dir=None, tmp_dir=None, logger=None):
        self.name = seqs_profile_info.name
        self.seq_type = seqs_profile_info.seq_type
        self.path = seqs_profile_info.path
        self.weight = seqs_profile_info.weight

        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir
        if tmp_dir is None:
            if self.path is not None:
                tmp_dir = self.path.split(".")[0] + "_%s_tmp" % generate_random_string(10)
            elif self.name is not None:
                tmp_dir = self.name.split(".")[0] + "_%s_tmp" % generate_random_string(10)
            else:
                tmp_dir = "seqs_profile_%s_tmp" % generate_random_string(10)

        self.method = method
        self.hmmer_inst_dir = hmmer_inst_dir
        self.infernal_inst_dir = infernal_inst_dir
        self.tmp_dir = tmp_dir
        self.logger = logger

    @classmethod
    def build(cls, mult_aln, name=None, path=None, weight=1.0,
              method=HMMER_KEY, hmmer_inst_dir=None, infernal_inst_dir=None, tmp_dir=None, logger=None,
              **kwargs):
        if name is None:
            if path is None:
                raise ValueError("")
            else:
                name = os.path.splitext(os.path.basename(path))[0]
        elif path is None:
            path = name
        if tmp_dir is None:
            tmp_dir = path.split(".")[0] + "_%s_tmp" % generate_random_string(10)

        if method.lower() == HMMER_KEY:
            os.makedirs(tmp_dir)
            if isinstance(mult_aln, str) and os.path.exists(mult_aln):
                mult_aln_path = mult_aln
            elif isinstance(mult_aln, MultAln):
                mult_aln_path = mult_aln.dump(os.path.join(tmp_dir, name)+".fasta", format="fasta")
            else:
                raise ValueError("the value for argument 'mult_aln' should be an instance of "
                                 "class eagle.lib.alignment.MultAln or path to a file with the alignment")
            hmmbuild_cmd = os.path.join(hmmer_inst_dir, "hmmbuild") + " " + path + " " + mult_aln_path
            subprocess.call(hmmbuild_cmd, shell=True)
            shutil.rmtree(tmp_dir, ignore_errors=True)

        if method.lower() == INFERNAL_KEY:
            pass

        return cls(SeqsProfileInfo(name=name, path=path, seq_type=mult_aln.seqs_type, weight=weight),
                   method=method, hmmer_inst_dir=hmmer_inst_dir, infernal_inst_dir=infernal_inst_dir,
                   tmp_dir=tmp_dir, logger=logger)

    def search(self, seqdb, out_path=None, threads=1, shred_seqdb=False, **kwargs):
        read_output = kwargs.get("read_output", False)
        if out_path is None:
            out_path = os.path.splitext(self.path)[0] + "_out_%s.psr" % generate_random_string(10)
            read_output = True

        if self.method.lower() == HMMER_KEY:
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)
            if shred_seqdb:  # I don't know if this this functional is needed for all methods
                prepared_seqdb = shred_seqs(seqdb, shredded_seqs_fasta=os.path.join(self.tmp_dir, "seqdb_shred.fasta"),
                                            part_l=50000, parts_ov=5000)
            else:
                prepared_seqdb = seqdb

            if isinstance(prepared_seqdb, SeqsDict):
                seqdb_path = prepared_seqdb.dump(os.path.join(self.tmp_dir, "seqdb_to_search.fasta"), format="fasta")
            else:
                seqdb_path = prepared_seqdb

            hmmsearch_cmd = os.path.join(self.hmmer_inst_dir, "hmmsearch") + \
                            " " + self.path + \
                            " " + seqdb_path + \
                            " --cpu " + str(threads)
            subprocess.call(hmmsearch_cmd, shell=True)
            shutil.rmtree(self.tmp_dir, ignore_errors=True)
            if read_output:
                with open(out_path) as out_f:
                    return Hmmer3TextParser(out_f)
            else:
                return out_path

        if self.method.lower() == INFERNAL_KEY:
            pass
        return

    @property
    def info(self):
        return SeqsProfileInfo(name=self.name, path=self.path, seq_type=self.seq_type, weight=self.weight)


class SeqProfilesDB(object):

    def __init__(self, name:str, method=HMMER_KEY, hmmer_inst_dir=None, infernal_inst_dir=None,
                 tmp_dir=None, logger=None):
        self.name = name

        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir
        if tmp_dir is None:
            tmp_dir = self.name.split(".")[0] + "_%s_tmp" % generate_random_string(10)

        self.method = method
        self.hmmer_inst_dir = hmmer_inst_dir
        self.infernal_inst_dir = infernal_inst_dir
        self.tmp_dir = tmp_dir
        self.logger = logger

    @classmethod
    def build(cls, profiles, name, method=HMMER_KEY, hmmer_inst_dir=None, infernal_inst_dir=None,
              tmp_dir=None, logger=None, **kwargs):
        profile_paths = list()
        for profile_info in profiles:
            if isinstance(profile_info, SeqsProfileInfo):
                profile_paths.append(profile_info.path)
            elif isinstance(profile_info, dict):
                profile_paths.append(SeqsProfileInfo.load_from_dict(profile_info).path)
            elif isinstance(profile_info, str) and os.path.exists(profile_info):
                profile_paths.append(profile_info)
        join_files(profile_paths, name)

        if method.lower() == HMMER_KEY:
            hmmpress_cmd = os.path.join(hmmer_inst_dir, "hmmpress") + " " + name
            subprocess.call(hmmpress_cmd, shell=True)

        if method.lower() == INFERNAL_KEY:
            pass

        return cls(name=name, method=method, hmmer_inst_dir=hmmer_inst_dir, infernal_inst_dir=infernal_inst_dir,
                   tmp_dir=tmp_dir, logger=logger)

    def scan(self, in_seqs, num_threads=4, out_path=None, shred_in_seqs=False, **kwargs):
        read_output = kwargs.get("read_output", False)
        if out_path is None:
            out_path = os.path.splitext(self.name)[0] + "_out_%s.psr" % generate_random_string(10)
            read_output = True

        if self.method.lower() == HMMER_KEY:
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)
            if shred_in_seqs:  # I don't know if this this functional is needed for all methods
                prepared_in_seqs = shred_seqs(in_seqs,
                                              shredded_seqs_fasta=os.path.join(self.tmp_dir, "in_seqs_shred.fasta"),
                                              part_l=50000, parts_ov=5000)
            else:
                prepared_in_seqs = in_seqs

            if isinstance(prepared_in_seqs, SeqsDict):
                in_seqs_path = prepared_in_seqs.dump(os.path.join(self.tmp_dir, "seqs_to_scan.fasta"), format="fasta")
            else:
                in_seqs_path = prepared_in_seqs
            hmmscan_cmd = os.path.join(self.hmmer_inst_dir, "hmmscan") + \
                          " --cpu " + str(num_threads) + " " + self.name + " " + \
                          in_seqs_path + " > " + out_path
            subprocess.call(hmmscan_cmd, shell=True)
            shutil.rmtree(self.tmp_dir, ignore_errors=True)
            if read_output:
                with open(out_path) as out_f:
                    return Hmmer3TextParser(out_f)
            else:
                return out_path

        if self.method.lower() == INFERNAL_KEY:
            pass
        return
