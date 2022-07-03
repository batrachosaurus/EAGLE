import os
import shutil
import subprocess
from collections.abc import Iterable

import pandas as pd

from jsondler import JsonEntry

from eaglib.constants import conf_constants
from eaglib._utils.files import join_files
from eaglib._utils.strings import generate_random_string
from eaglib.seqs import SeqsDict, shred_seqs
from eaglib.alignment.mult_aln import MultAln


HMMER_KEY = "hmmer"
INFERNAL_KEY = "infernal"


class SeqsProfileInfo(JsonEntry):

    # json keys
    name_key = "name"
    path_key = "path"
    seq_type_key = "type"
    weight_key = "weight"
    method_key = "method"

    # default values
    name_0 = None
    path_0 = None
    seq_type_0 = 'protein'
    weight_0 = 1.0
    method_0 = HMMER_KEY

    def __init__(self,
                 name=name_0,
                 path=path_0,
                 seq_type=seq_type_0,
                 weight=weight_0,
                 method=method_0):

        # attribute names must match keys form SeqProfileInfo.attr_scheme()
        self.name = name
        self.path = path
        self.seq_type = seq_type
        self.weight = weight
        self.method = method

    @classmethod
    def attr_scheme(cls):
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "name": (cls.name_key,),
            "path": (cls.path_key,),
            "seq_type": (cls.seq_type_key,),
            "weight": (cls.weight_key,),
            "method": (cls.method_key,),
        }


class SeqsProfile(object):

    def __init__(self, seqs_profile_info: SeqsProfileInfo,
                 hmmer_inst_dir=None, infernal_inst_dir=None, tmp_dir=None, **kwargs):
        self.name = seqs_profile_info.name
        self.seq_type = seqs_profile_info.seq_type
        self.path = seqs_profile_info.path
        self.weight = seqs_profile_info.weight
        self.method = seqs_profile_info.method

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

        self.hmmer_inst_dir = hmmer_inst_dir
        self.infernal_inst_dir = infernal_inst_dir
        self.tmp_dir = tmp_dir

    @classmethod
    def build(cls, mult_aln, name=None, path=None, weight=1.0,
              method=HMMER_KEY, hmmer_inst_dir=None, infernal_inst_dir=None, tmp_dir=None,
              **kwargs):
        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir

        if name is None:
            if path is None:
                raise ValueError("")
            else:
                name = os.path.splitext(os.path.basename(path))[0]
        elif path is None:
            if method.lower() == HMMER_KEY:
                path = name + ".hmm"
            elif method.lower() == INFERNAL_KEY:
                path = name + ".cm"
            else:
                path = name
        if tmp_dir is None:
            tmp_dir = path.split(".")[0] + "_%s_tmp" % generate_random_string(10)

        if method.lower() in (HMMER_KEY, INFERNAL_KEY):
            os.makedirs(tmp_dir)
            if isinstance(mult_aln, str) and os.path.exists(mult_aln):
                mult_aln_path = mult_aln
                seqs_type = kwargs.get("seqs_type", SeqsProfileInfo.seq_type_0)
            elif isinstance(mult_aln, MultAln):
                mult_aln_path = mult_aln.dump(os.path.join(tmp_dir, name)+".fasta", format="fasta")
                seqs_type = kwargs.get("seqs_type", mult_aln.seqs_type)
            else:
                raise ValueError("the value for argument 'mult_aln' should be an instance of "
                                 "class eagle.eaglib.alignment.MultAln or path to a file with the alignment")

            if method.lower() == HMMER_KEY:
                build_cmd = os.path.join(hmmer_inst_dir, "hmmbuild ") + path + " " + mult_aln_path
            if method.lower() == INFERNAL_KEY:
                if kwargs.get("noss", False):
                    build_cmd = os.path.join(infernal_inst_dir, "cmbuild --noss ") + path + " " + mult_aln_path
                else:
                    build_cmd = os.path.join(infernal_inst_dir, "cmbuild ") + path + " " + mult_aln_path + "&&" + \
                                os.path.join(infernal_inst_dir, "cmcalibrate ") + path
            subprocess.call(build_cmd, shell=True)

            shutil.rmtree(tmp_dir, ignore_errors=True)

        return cls(SeqsProfileInfo(name=name, path=path, seq_type=seqs_type, weight=weight, method=method),
                   hmmer_inst_dir=hmmer_inst_dir, infernal_inst_dir=infernal_inst_dir, tmp_dir=tmp_dir)

    def search(self, seqdb, out_path=None, threads=1, **kwargs):
        # All preparations of seqdb like translation or shredding should be done outside the method
        read_output = kwargs.get("read_output", False)
        if out_path is None:
            out_path = os.path.splitext(self.path)[0] + "_out_%s.psr" % generate_random_string(10)
            read_output = True

        if self.method.lower() in (HMMER_KEY, INFERNAL_KEY):
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)

            # The example of shred_seqs function usage, remove this after reimplementation in other modules
            # prepared_seqdb = shred_seqs(seqdb, shredded_seqs_fasta=os.path.join(self.tmp_dir, "seqdb_shred.fasta"),
            #                             part_l=50000, parts_ov=5000)

            if isinstance(seqdb, SeqsDict):
                seqdb_path = seqdb.dump(os.path.join(self.tmp_dir, "seqdb_to_search.fasta"), format="fasta")
            else:
                seqdb_path = seqdb

            if self.method.lower() == HMMER_KEY:
                search_cmd = os.path.join(self.hmmer_inst_dir, "hmmsearch --cpu ") + \
                             str(threads) + " --domtblout " + out_path + \
                             " " + self.path + " " + seqdb_path
                read_out_func = read_hmmer_domtblout
            if self.method.lower() == INFERNAL_KEY:
                search_cmd = os.path.join(self.infernal_inst_dir, "cmsearch --cpu ") + \
                             str(threads) + " --tblout " + out_path + \
                             " " + self.path + " " + seqdb_path
                read_out_func = read_infernal_tblout
            subprocess.call(search_cmd, shell=True)
            shutil.rmtree(self.tmp_dir, ignore_errors=True)
            if read_output:
                with open(out_path) as out_f:
                    return read_out_func(out_path)
            else:
                return out_path
        return

    @property
    def info(self):
        return SeqsProfileInfo(name=self.name, path=self.path, seq_type=self.seq_type, weight=self.weight,
                               method=self.method)


class SeqProfilesDB(object):

    def __init__(self, name:str, method=HMMER_KEY, hmmer_inst_dir=None, infernal_inst_dir=None, tmp_dir=None, **kwargs):
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

    @classmethod
    def build(cls, profiles, name, method=HMMER_KEY, hmmer_inst_dir=None, infernal_inst_dir=None,
              tmp_dir=None, **kwargs):
        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir

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
            hmmpress_cmd = os.path.join(hmmer_inst_dir, "hmmpress ") + name
            subprocess.call(hmmpress_cmd, shell=True)
        if method.lower() == INFERNAL_KEY:
            cmpress_cmd = os.path.join(infernal_inst_dir, "cmmpress ") + name
            subprocess.call(cmpress_cmd, shell=True)

        return cls(name=name, method=method, hmmer_inst_dir=hmmer_inst_dir, infernal_inst_dir=infernal_inst_dir,
                   tmp_dir=tmp_dir)

    def scan(self, in_seqs, num_threads=4, out_path=None, **kwargs):
        read_output = kwargs.get("read_output", False)
        if out_path is None:
            out_path = os.path.splitext(self.name)[0] + "_out_%s.psr" % generate_random_string(10)
            read_output = True

        if self.method.lower() in (HMMER_KEY, INFERNAL_KEY):
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)

            if isinstance(in_seqs, SeqsDict):
                in_seqs_path = in_seqs.dump(os.path.join(self.tmp_dir, "seqs_to_scan.fasta"), format="fasta")
            else:
                in_seqs_path = in_seqs

            if self.method.lower() == HMMER_KEY:
                scan_cmd = os.path.join(self.hmmer_inst_dir, "hmmscan --cpu ") + \
                           str(num_threads) + " --domtblout " + out_path + \
                           " " + self.name + " " + in_seqs_path
                read_out_func = read_hmmer_domtblout
            if self.method.lower() == INFERNAL_KEY:
                scan_cmd = os.path.join(self.infernal_inst_dir, "cmscan --cpu ") + \
                           str(num_threads) + " --tblout " + out_path + \
                           " " + self.name + " " + in_seqs_path
                read_out_func = read_infernal_tblout
            subprocess.call(scan_cmd, shell=True)
            shutil.rmtree(self.tmp_dir, ignore_errors=True)
            if read_output:
                with open(out_path) as out_f:
                    return read_out_func(out_path)
            else:
                return out_path
        return


def read_infernal_tblout(tblout_path):
    with open(tblout_path) as tblout_f:
        columns_markup = list()
        hits = list()
        for i, line in enumerate(tblout_f):
            if i == 0:
                columns_line = " " + line[1:]
            if i == 1:
                columns_markup = get_columns_markup(line)
                columns = [columns_line[col_cs[0]: col_cs[1]].strip() for col_cs in columns_markup]
            if line[0] != "#" and columns_markup:
                hits.append([line[col_cs[0]: col_cs[1]].strip() for col_cs in columns_markup])
    return pd.DataFrame(hits, columns=columns)


def read_hmmer_domtblout(domtblout_path):
    with open(domtblout_path) as domtblout_f:
        columns_markup = list()
        hits = list()
        for i, line in enumerate(domtblout_f):
            if i == 1:
                columns_line = " " + line[1:]
            if i == 2:
                columns_markup = get_columns_markup(line)
                columns = _fix_hmmer_columns(columns_line[col_cs[0]: col_cs[1]].strip() for col_cs in columns_markup)
            if line[0] != "#" and columns_markup:
                hits.append([line[col_cs[0]: col_cs[1]].strip() for col_cs in columns_markup])
    return pd.DataFrame(hits, columns=columns)


def _fix_hmmer_columns(columns: Iterable):
    fixed_cols = list()
    score_met = 0
    bias_met = 0
    from_met = 0
    to_met = 0
    for col in columns:
        if col == "score":
            if score_met == 0:
                fixed_cols.append("seq_score")
                score_met += 1
            elif score_met == 1:
                fixed_cols.append("domain_score")
        elif col == "bias":
            if bias_met == 0:
                fixed_cols.append("seq_bias")
                bias_met += 1
            elif bias_met == 1:
                fixed_cols.append("domain_bias")
        elif col == "from":
            if from_met == 0:
                fixed_cols.append("hmm_from")
                from_met += 1
            elif from_met == 1:
                fixed_cols.append("ali_from")
                from_met += 1
            elif from_met == 2:
                fixed_cols.append("env_from")
        elif col == "to":
            if to_met == 0:
                fixed_cols.append("hmm_to")
                to_met += 1
            elif to_met == 1:
                fixed_cols.append("ali_to")
                to_met += 1
            elif to_met == 2:
                fixed_cols.append("env_to")
        else:
            fixed_cols.append(col)
    return fixed_cols


def get_columns_markup(markup_line, line_max_len=1000):
    prev_s = None
    opened = True
    columns_markup = [[0, 1]]
    for i, s in enumerate(markup_line):
        if s == " " and prev_s == "-":
            columns_markup[-1][1] = i+1
            opened = False
        if not opened and s == "-" and prev_s == " ":
            columns_markup.append([i, i+1])
        prev_s = s
    columns_markup[-1][1] = max(columns_markup[-1][1], line_max_len)
    return columns_markup
