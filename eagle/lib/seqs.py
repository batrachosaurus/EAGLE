import os
import shutil
from collections import defaultdict, Counter
import subprocess

import psutil
import numpy as np
from Bio.Seq import Seq

from eagle.constants import conf_constants
from eagle.lib.general import generate_random_string, np_memmap_astype
from eagle.lib.general import filter_list, get_un_fix


class SeqsDict(object):

    def __init__(self, seqs_order, seqs_array, seqs_type=None, low_memory=False):
        self.seqs_order = seqs_order
        self.seqs_array = seqs_array
        self.seqs_type = seqs_type
        self.low_memory = low_memory

        self._empty_rows = list()

    def __getitem__(self, item):
        return self.seqs_array[self.seqs_order[item]]

    def get_sample(self, seqs, low_memory='auto', **kwargs):
        if low_memory == "auto":
            if self.low_memory:
                low_memory = True
            else:
                low_memory = False
        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=self.seqs_array.dtype, mode='w+', shape=len(seqs))
        else:
            seqs_array = np.zeros(len(seqs), dtype=self.seqs_array.dtype)
        seqs_order = dict()
        for i, seq in enumerate(seqs):
            seqs_order[seq] = i
            seqs_array[i] = self[seq]
        return SeqsDict(seqs_order=seqs_order, seqs_array=seqs_array, low_memory=low_memory)

    def __setitem__(self, key, value):
        if len(value) > self.seqs_array.itemsize:
            if self.low_memory:
                self.seqs_array = np_memmap_astype(dat_path=self.seqs_array.filename,
                                                   old_dtype=self.seqs_array.dtype,
                                                   new_dtype=np.dtype("S%s" % len(value)),
                                                   shape=self.seqs_array.shape)
            else:
                self.seqs_array = self.seqs_array.astype(np.dtype("S%s" % len(value)))
        if self._empty_rows:
            self.seqs_order[key] = self._empty_rows.pop(0)
            self.seqs_array[self.seqs_order[key]] = value
        elif key in self.seqs_order:
            self.seqs_array[self.seqs_order[key]] = value
        else:
            self.seqs_order[key] = self.seqs_array.shape[0]
            if self.low_memory:
                self.seqs_array = np.memmap(self.seqs_array.filename,
                                            dtype=self.seqs_array.dtype,
                                            mode='r+',
                                            shape=self.seqs_array.shape[0]+1,
                                            order='C')
                self.seqs_array[-1] = value
            else:
                self.seqs_array = np.concatenate((self.seqs_array, np.array([value], dtype=self.seqs_array.dtype)))

    def __len__(self):
        return len(self.seqs_order)

    def __iter__(self):
        for key in self.keys():
            yield key

    def pop(self, key):
        seq = self[key]
        del self[key]
        return seq

    def __delitem__(self, key):
        if self.low_memory:
            self.seqs_array[self.seqs_order[key]] = ""
            self._empty_rows.append(self.seqs_order.pop(key))
        else:
            self.seqs_array = np.concatenate((self.seqs_array[:self.seqs_order[key]],
                                              self.seqs_array[self.seqs_order[key]+1:]))
            i = 0
            for seq_id in self.keys():
                if seq_id == key:
                    del self.seqs_order[seq_id]
                else:
                    self.seqs_order[seq_id] = i
                    i += 1

    def keys(self):
        return map(lambda si: si[0], sorted(self.seqs_order.items(), key=lambda si: si[1]))

    def values(self):
        return self.seqs_array

    def items(self):
        return map(lambda seq_id: (seq_id, self[seq_id]), self.keys())

    def rename_seqs(self, old_to_new_dict):
        for old_seq in old_to_new_dict:
            self.seqs_order[old_to_new_dict[old_seq]] = self.seqs_order.pop(old_seq)

    @classmethod
    def load_from_file(cls, seqs_path, seqs_format="fasta", low_memory='auto', **kwargs):
        read_chunk_size = kwargs.get("read_chink_size", 100)
        exp_seq_len = kwargs.get("exp_seq_len", 1000)
        if low_memory == 'auto':
            low_memory = cls._check_low_memory(seqs_path=seqs_path)
        seq_i = 0
        seqs_order = dict()
        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=np.dtype("S%s" % exp_seq_len), mode='w+', shape=read_chunk_size)
        else:
            seqs_array = np.zeros(read_chunk_size, dtype=np.dtype("S%s" % exp_seq_len))
        if seqs_format.lower() == "fasta":
            seq_list = list()
            title = None
            fasta_f = open(seqs_path)
            for line_ in fasta_f:
                line = None
                line = line_.strip()
                if not line:
                    continue
                if line[0] == ">":
                    if title:
                        seqs_order[title] = seq_i
                        new_seq = None
                        if kwargs.get("restore_stops", False):
                            new_seq = "".join(seq_list).replace("X", "*")
                        else:
                            new_seq = "".join(seq_list)
                        if len(new_seq) > seqs_array.itemsize:
                            if low_memory:
                                seqs_array = np_memmap_astype(dat_path=seqs_array.filename,
                                                              old_dtype=seqs_array.dtype,
                                                              new_dtype=np.dtype("S%s" % len(new_seq)),
                                                              shape=seqs_array.shape)
                            else:
                                seqs_array = seqs_array.astype(np.dtype("S%s" % len(new_seq)))
                        seqs_array[seq_i] = new_seq
                        title = None
                        seq_i += 1
                    if seq_i > 0 and float(seq_i) % float(read_chunk_size) == 0:
                        if low_memory:
                            seqs_array = np.memmap(seqs_array.filename,
                                                   dtype=seqs_array.dtype,
                                                   mode='r+',
                                                   shape=seqs_array.shape[0]+read_chunk_size,
                                                   order='C')
                        else:
                            seqs_array = np.concatenate((seqs_array, np.zeros(read_chunk_size, dtype=seqs_array.dtype)))
                    title = line[1:]
                    seq_list = list()
                else:
                    seq_list.append(line)
            if title:
                seqs_order[title] = seq_i
                new_seq = None
                new_seq = None
                if kwargs.get("restore_stops", False):
                    new_seq = "".join(seq_list).replace("X", "*")
                else:
                    new_seq = "".join(seq_list)
                if len(new_seq) > seqs_array.itemsize:
                    seqs_array = seqs_array.astype(np.dtype("S%s" % len(new_seq)))
                seqs_array[seq_i] = new_seq
                seq_list = list()
                title = None
                seq_i += 1
            fasta_f.close()
        if seqs_array.shape[0] > seq_i:
            seqs_array = seqs_array[:seq_i]
        return cls(seqs_order=seqs_order, seqs_array=seqs_array, low_memory=low_memory)

    @classmethod
    def load_from_dict(cls, in_dict, low_memory=False, **kwargs):
        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=np.dtype("S1000"), mode='w+', shape=len(in_dict))
        else:
            seqs_array = np.zeros(len(in_dict), dtype=np.dtype("S1000"))
        seqs_order = dict()
        for seq_i, seq_id in enumerate(in_dict):
            seqs_order[seq_id] = seq_i
            if len(in_dict[seq_id]) > seqs_array.itemsize:
                seqs_array = seqs_array.astype(dtype=np.dtype("S%s" % len(in_dict[seq_id])))
            seqs_array[seq_i] = in_dict[seq_id]
        return cls(seqs_order=seqs_order, seqs_array=seqs_array, low_memory=low_memory)

    @staticmethod
    def _check_low_memory(seqs_path):
        avail_mem = psutil.virtual_memory().total
        seqs_size = os.path.getsize(seqs_path)
        if seqs_size == 0:
            seqs_size = 1
        if float(avail_mem) / float(seqs_size) < 0.5:
            low_memory = True
        else:
            low_memory = False
        return low_memory

    def dump(self, seqs_path, seqs_format="fasta", overwrite=True, **kwargs):
        if seqs_format.lower() == "fasta":
            if overwrite:
                fasta_f = open(seqs_path, 'w')
            else:
                fasta_f = open(seqs_path, 'a')
            if kwargs.get("replace_stops", False):
                for seq_id in self.seqs_order:
                    fasta_f.write(">" + seq_id + "\n")
                    fasta_f.write(self[seq_id].replace("*", "X") + "\n")
            else:
                for seq_id in self.seqs_order:
                    fasta_f.write(">" + seq_id + "\n")
                    fasta_f.write(self[seq_id] + "\n")
            fasta_f.close()

    def detect_seqs_type(self, nuc_freq_thr=0.75):
        seqs_list = list()
        summ_l = 0
        for seq in self.seqs_array:
            seqs_list.append(seq.lower().replace("-", ""))
            summ_l += len(seqs_list[-1])
        let_counts = Counter("".join(seqs_list))
        if float(let_counts.get("a", 0) + let_counts.get("c", 0) + let_counts.get("g", 0) +
                 let_counts.get("t", 0)) / float(summ_l) >= nuc_freq_thr:
            return "nucl"
        else:
            return "prot"


def seq_from_fasta(fasta_path, seq_id, ori=+1, start=1, end=-1):
    fasta_dict = load_fasta_to_dict(fasta_path)
    if end == -1:
        end = len(fasta_dict[seq_id])
    if start < 0:
        start = len(fasta_dict[seq_id]) + start + 1
    if end >= start:
        if ori > 0:
            return fasta_dict[seq_id][start-1: end]
        else:
            return str(Seq(fasta_dict[seq_id][start-1: end]).reverse_complement())
    else:
        if ori > 0:
            return fasta_dict[seq_id][end-1: start]
        else:
            return str(Seq(fasta_dict[seq_id][end-1: start]).reverse_complement())


def shred_seqs(fasta_dict, part_l=50000, parts_ov=5000):
    shredded_seqs = defaultdict(list)
    for seq_id in fasta_dict:
        i = 0
        l_seq = len(fasta_dict[seq_id])
        while i < l_seq:
            if i+part_l < l_seq:
                shredded_seqs[seq_id].append(fasta_dict[seq_id][i: i+part_l])
                last_ov_c = i + part_l + int(parts_ov/2)
                if last_ov_c < l_seq:
                    shredded_seqs[seq_id].append(fasta_dict[seq_id][i+part_l-int(parts_ov/2): last_ov_c])
                else:
                    shredded_seqs[seq_id].append(fasta_dict[seq_id][i+part_l-int(parts_ov/2):])
            else:
                shredded_seqs[seq_id].append(fasta_dict[seq_id][i:])
            i += part_l
    return shredded_seqs


def load_fasta_to_dict(fasta_path, low_memory="auto", **kwargs):
    return SeqsDict.load_from_file(seqs_path=fasta_path, seqs_format="fasta", low_memory=low_memory, **kwargs)


def dump_fasta_dict(fasta_dict, fasta_path, overwrite=True, **kwargs):
    if isinstance(fasta_dict, dict):
        fasta_dict = SeqsDict.load_from_dict(in_dict=fasta_dict)
    fasta_dict.dump(seqs_path=fasta_path, seqs_format="fasta", overwrite=overwrite, **kwargs)


def reduce_seq_names(fasta_dict, num_letters=10, num_words=4):
    if num_letters < 6:
        print("Number of letters must be at least 6")
        return 1
    if num_words < 2:
        print("Number of words must be at least 2")
        return 1
    all_less = True
    for seq_id in fasta_dict:
        if len(seq_id) > num_letters:
            all_less = False
            break
    if all_less:
        return fasta_dict, dict((seq_id, seq_id) for seq_id in fasta_dict)
    splitters_repl = {"_": " ",
                      "\t": " ",
                      ",": " ",
                      ";": " ",
                      ".": " ",
                      ":": " ",
                      "|": " ",
                      "/": " ",
                      "\\": " ",
                      "(": " ",
                      ")": " ",
                      "[": " ",
                      "]": " ",
                      "{": " ",
                      "}": " "}
    parts_size_list = _get_part_size_list(num_letters, num_words)
    reduced_fasta_dict = dict()
    seq_names_dict = dict()
    for seq_name in fasta_dict.keys():
        reduced_seq_name = None
        seq_name_list = filter_list("".join([splitters_repl.get(s, s) for s in seq_name]).split())
        parts = list()
        for i in range(num_words):
            try:
                parts.append(seq_name_list[i][:parts_size_list[i]])
            except IndexError:
                break
        reduced_seq_name = "".join(parts)
        res_len = num_letters - len(reduced_seq_name)
        un_num = 0
        un_fix = get_un_fix(un_num, res_len)
        while seq_names_dict.get((reduced_seq_name+un_fix).upper(), None):
            un_fix = None
            un_num += 1
            un_fix = get_un_fix(un_num, res_len)
        reduced_fasta_dict[(reduced_seq_name+un_fix).upper()] = fasta_dict[seq_name]
        seq_names_dict[(reduced_seq_name+un_fix).upper()] = seq_name
    return reduced_fasta_dict, seq_names_dict


def _get_part_size_list(num_letters, num_words):
    if num_letters == 6:
        return [2, 3]
    if num_letters == 7:
        if num_words >= 3:
            return [2, 3, 1]
        else:
            return [2, 3]
    if num_letters == 8:
        if num_words >= 4:
            return [2, 3, 1, 1]
        elif num_words == 3:
            return [3, 3, 1]
        else:
            return [3, 3]
    if num_letters == 9:
        if num_words >= 4:
            return [3, 3, 1, 1]
        elif num_words == 3:
            return [3, 3, 1]
        else:
            return [3, 4]
    if num_letters == 10:
        if num_words >= 4:
            return [3, 4, 1, 1]
        elif num_words == 3:
            return [3, 4, 1]
        else:
            return [4, 4]
    if num_letters >= 11:
        if num_words >= 4:
            return [3, 4, 1, 1]
        elif num_words == 3:
            return [4, 4, 1]
        else:
            return [4, 5]


def get_orfs(in_fasta_path, out_fasta_path, minsize=180, emboss_inst_dir=conf_constants.emboss_inst_dir):
    orfs_info = dict()
    subprocess.call(os.path.join(emboss_inst_dir, "getorf") + " " + in_fasta_path + " " + out_fasta_path + " -minsize "
                    + str(minsize), shell=True)
    orfs_fasta_dict = load_fasta_to_dict(out_fasta_path)
    corr_orfs_ids = dict()
    orfs_ids = orfs_fasta_dict.keys()
    for orf_id in orfs_ids:
        corr_orf_id = None
        ori = None
        corr_orf_id, c_start, c_end, ori = _get_orf_info(orf_id)
        corr_orfs_ids[orf_id] = corr_orf_id
        orfs_info[corr_orf_id] = {
            "seqid": corr_orf_id,
            "source": "EAGLE",
            "type": "ORF",
            "start": c_start,
            "end": c_end,
            "score": "-",
            "strand": ori,
            "frame": ".",
            "attribute": dict(),
        }
    orfs_fasta_dict.rename_seqs(corr_orfs_ids)
    dump_fasta_dict(orfs_fasta_dict, out_fasta_path)
    return orfs_info


def _get_orf_info(orf_title):
    orf_title_list = orf_title.split()
    c1 = int(orf_title_list[1].strip("[]"))
    c2 = int(orf_title_list[3].strip("[]"))
    if c2 >= c1:
        ori = "+"
        c_start = c1
        c_end = c2
        orf_id = "_".join(orf_title_list[0].split("_")[:-1]) + "|:" + str(c1) + "-" + str(c2)
    else:
        ori = "-"
        c_start = c2
        c_end = c1
        orf_id = "_".join(orf_title_list[0].split("_")[:-1]) + "|:c" + str(c1) + "-" + str(c2)
    return orf_id, c_start, c_end, ori


def parse_orf_id(orf_id):
    orf_id_list = orf_id.split("|:")
    orf_chr = orf_id_list[0]
    if orf_id_list[1][0] == "c":
        ori = "-"
        c_list = orf_id_list[1][1:].split("-")
        c2 = c_list[0]
        c1 = c_list[1]
    else:
        ori = "+"
        c_list = orf_id_list[1].split("-")
        c1 = c_list[0]
        c2 = c_list[1]
    return orf_chr, int(c1), int(c2), ori


def read_blast_out(blast_out_path, ev_thr=1.0e-06, aln_l_thr=180, ident_thr=0.35):
    """
    Reads blast outfmt 6 or 7
    :param blast_out_path:
    :param ev_thr:
    :param aln_l_thr:
    :param ident_thr:
    :return:
    """
    blast_res_dict = defaultdict(list)
    blast_out_f = open(blast_out_path)
    for line_ in blast_out_f:
        orf_id = None
        line = None
        line = line_.strip()
        if not line:
            continue
        if line[0] == "#":
            continue
        line_list = line.split("\t")
        ev = float(line_list[10].strip())
        aln_l = abs(int(line_list[9].strip())-int(line_list[8].strip())) + 1
        ident = float(line_list[2].strip())/100.0
        if ev <= ev_thr and aln_l >= aln_l_thr and ident >= ident_thr:
            blast_res_dict[line_list[0].strip()].append({
                "subj_id": line_list[1].strip(),
                "identity": ident,
                "aln_l": aln_l,
                "evalue": ev,
                "subj_start": int(line_list[8].strip()),
                "subj_end": int(line_list[9].strip()),
            })
    return blast_res_dict
