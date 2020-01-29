import os
import shutil
from copy import deepcopy
from collections import defaultdict, Counter, OrderedDict
import subprocess

import psutil
import numpy as np
from Bio.Seq import Seq

from eagle.constants import conf_constants
from eagle.lib.general import generate_random_string, np_memmap_astype, filter_list, get_un_fix


class SeqsDict(object):

    prot_type = "prot"
    nucl_type = "nucl"
    prot_types = ("protein", prot_type, "p")
    nucl_types = ("nucleotide", nucl_type, "n")
    _chunk_size0 = 100

    def __init__(self, seqs_order, seqs_array, seq_info_dict=None, seqs_type=None, low_memory=False,
                 chunk_size=_chunk_size0):
        if seq_info_dict is None:
            seq_info_dict = defaultdict(lambda: defaultdict(lambda: None))

        self.seqs_order = seqs_order
        self.seqs_array = seqs_array
        self.seq_info_dict = seq_info_dict
        self.seqs_type = seqs_type
        self.low_memory = low_memory
        self._chunk_size = chunk_size

        self._empty_rows = list()  # TODO: implement special object

    def __getitem__(self, item):
        return "".join(self.seqs_array[i].decode() for i in self.seqs_order[item])

    def get_sample(self, seqs, low_memory='auto', **kwargs):
        if low_memory == "auto":
            if self.low_memory:
                low_memory = True
            else:
                low_memory = False
        seqs_array_l = sum(len(self.seqs_order[seq]) for seq in seqs)
        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=self.seqs_array.dtype, mode='w+', shape=(seqs_array_l,))
        else:
            seqs_array = np.zeros(seqs_array_l, dtype=self.seqs_array.dtype)
        seqs_order = dict()
        i = 0
        for seq in seqs:
            r = self.seqs_order[seq]
            seqs_order[seq] = range(i, i+len(r))
            for n, j in enumerate(seqs_order[seq]):
                seqs_array[j] = self.seqs_array[r.start+n]
            i += len(r)
        return SeqsDict(seqs_order=seqs_order, seqs_array=seqs_array, low_memory=low_memory,
                        chunk_size=self._chunk_size)

    def __setitem__(self, key, value):
        # TODO: implement self._empty_rows
        """
        if self._empty_rows:
            self.seqs_order[key] = self._empty_rows.pop(0)
            self.seqs_array[self.seqs_order[key]] = value
        """
        # elif key in self.seqs_order:
        if key in self.seqs_order:
            del self[key]
            self[key] = value
        else:
            self.seqs_order[key] = range(self.seqs_array.shape[0],
                                         self.seqs_array.shape[0]+len(value)//self._chunk_size+1)
            if self.low_memory:
                self.seqs_array = np.memmap(self.seqs_array.filename,
                                            dtype=self.seqs_array.dtype,
                                            mode='r+',
                                            shape=(self.seqs_order[key].stop,),
                                            order='C')
                for n, i in enumerate(self.seqs_order[key]):
                    self.seqs_array[i] = value[n*self._chunk_size: (n+1)*self._chunk_size]
            else:
                self.seqs_array = np.concatenate(
                    (self.seqs_array,
                     np.array([value[n*self._chunk_size: (n+1)*self._chunk_size] for n, i in enumerate(self.seqs_order[key])],
                              dtype=self.seqs_array.dtype))
                )

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
            self._empty_rows.append(self.seqs_order.pop(key))
        else:
            self.seqs_array = np.concatenate((self.seqs_array[:self.seqs_order[key].start],
                                              self.seqs_array[self.seqs_order[key].stop:]))
            i = 0
            for seq_id in self.keys():
                if seq_id == key:
                    del self.seqs_order[seq_id]
                else:
                    self.seqs_order[seq_id] = i
                    i += 1

    def __del__(self):
        if self.low_memory:
            dat_path = self.seqs_array.filename
            del self.seqs_array
            os.remove(dat_path)

    def __contains__(self, item):
        return item is self.seqs_order

    def keys(self):
        return map(lambda si: si[0], sorted(self.seqs_order.items(), key=lambda si: si[1].start))

    def values(self):
        for seq_name in self:
            yield self[seq_name]

    def items(self):
        return map(lambda seq_id: (seq_id, self[seq_id]), self.keys())

    def rename_seqs(self, old_to_new_dict):
        for old_seq in old_to_new_dict:
            if old_seq in self.seqs_order:
                self.seqs_order[old_to_new_dict[old_seq]] = self.seqs_order.pop(old_seq)

    def copy(self):
        # TODO: implement _sub_seqs_dict method
        return deepcopy(self)

    @classmethod
    def load_from_file(cls, seqs_path, seqs_format="fasta", low_memory='auto', **kwargs):
        chunk_size = kwargs.get("chunk_size", cls._chunk_size0)
        if low_memory == 'auto':
            low_memory = cls._check_low_memory(seqs_path=seqs_path)

        n_chunks = 0
        if seqs_format.lower() == "fasta":
            with open(seqs_path) as fasta_f:
                seq_l = 0
                for line_ in fasta_f:
                    line = None
                    line = line_.strip()
                    if not line:
                        continue
                    if line[0] == ">" and seq_l > 0:
                        n_chunks += (seq_l-1)//chunk_size + 1
                        seq_l = 0
                    else:
                        seq_l += len(line)

        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=np.dtype("S%s" % chunk_size), mode='w+', shape=(n_chunks,))
        else:
            seqs_array = np.zeros(n_chunks, dtype=np.dtype("S%s" % chunk_size))

        i = 0
        seqs_order = dict()
        if seqs_format.lower() == "fasta":
            with open(seqs_path) as fasta_f:
                seq_list = list()
                title = None
                for line_ in fasta_f:
                    line = None
                    line = line_.strip()
                    if not line:
                        continue
                    if line[0] == ">":
                        if title is not None:
                            new_seq = None
                            if kwargs.get("restore_stops", False):
                                new_seq = "".join(seq_list).replace("X", "*")
                            else:
                                new_seq = "".join(seq_list)
                            seqs_order[title] = range(i, i + (len(new_seq)-1)//chunk_size + 1)
                            for n, j in enumerate(seqs_order[title]):
                                seqs_array[j] = new_seq[n*chunk_size: (n+1)*chunk_size]
                            i = seqs_order[title].stop
                            title = None
                        title = line[1:]
                        seq_list = list()
                    else:
                        seq_list.append(line)
                if title is not None:
                    new_seq = None
                    if kwargs.get("restore_stops", False):
                        new_seq = "".join(seq_list).replace("X", "*")
                    else:
                        new_seq = "".join(seq_list)
                    seqs_order[title] = range(i, i + (len(new_seq)-1) // chunk_size + 1)
                    for n, j in enumerate(seqs_order[title]):
                        seqs_array[j] = new_seq[n * chunk_size: (n + 1) * chunk_size]
                    i = seqs_order[title].stop
                    seq_list = list()
                    title = None
        return cls(seqs_order=seqs_order, seqs_array=seqs_array, low_memory=low_memory, chunk_size=chunk_size)

    @classmethod
    def load_from_dict(cls, in_dict, low_memory=False, **kwargs):
        chunk_size = kwargs.get("chunk_size", cls._chunk_size0)
        n_chunks = 0
        for seq in in_dict.values():
            if seq:
                n_chunks += (len(seq)-1)//chunk_size + 1

        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=np.dtype("S%s" % chunk_size), mode='w+', shape=(n_chunks,))
        else:
            seqs_array = np.zeros(n_chunks, dtype=np.dtype("S%s" % chunk_size))

        i = 0
        seqs_order = dict()
        for seq_id, seq in in_dict.items():
            seqs_order[seq_id] = range(i, i + (len(seq)-1) // chunk_size + 1)
            for n, j in enumerate(seqs_order[seq_id]):
                seqs_array[j] = seq[n * chunk_size: (n + 1) * chunk_size]
            i = seqs_order[seq_id].stop
        return cls(seqs_order=seqs_order, seqs_array=seqs_array, low_memory=low_memory, chunk_size=chunk_size)

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
        return seqs_path

    def detect_seqs_type(self, nuc_freq_thr=0.75):
        summ_l = 0
        let_counts = Counter()
        for seq_ in self.values():
            seq = None
            seq = seq_.lower().replace("-", "")
            summ_l += len(seq)
            let_counts += Counter(seq)
        if float(let_counts.get("a", 0) + let_counts.get("c", 0) + let_counts.get("g", 0) +
                 let_counts.get("t", 0)) / float(summ_l) >= nuc_freq_thr:
            return self.nucl_type
        else:
            return self.prot_type


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


def shred_fasta(in_fasta, shredded_fasta_path, part_l=50000, parts_ov=5000):
    if isinstance(in_fasta, SeqsDict):
        in_seqs_dict = in_fasta
    else:
        in_seqs_dict = load_fasta_to_dict(fasta_path=in_fasta)

    shredded_in_seqs = shred_seqs(seqs_dict=in_seqs_dict, part_l=part_l, parts_ov=parts_ov)
    seqs_to_scan_dict = OrderedDict()
    for seq_id in shredded_in_seqs:
        i = 0
        for seq in shredded_in_seqs[seq_id]:
            start = i*(part_l-parts_ov)
            seqs_to_scan_dict[seq_id + "|:" + str(start+1) + "-" + str(start+part_l)] = seq
            i += 1
    dump_fasta_dict(fasta_dict=seqs_to_scan_dict, fasta_path=shredded_fasta_path)
    return shredded_fasta_path


def shred_seqs(seqs_dict, part_l=50000, parts_ov=5000):
    shredded_seqs = defaultdict(list)
    for seq_id in seqs_dict:
        i = 0
        l_seq = len(seqs_dict[seq_id])
        while i < l_seq:
            if i+part_l < l_seq:
                shredded_seqs[seq_id].append(seqs_dict[seq_id][i: i + part_l])
                last_ov_c = i + part_l + int(parts_ov/2)
                if last_ov_c < l_seq:
                    shredded_seqs[seq_id].append(seqs_dict[seq_id][i + part_l - int(parts_ov / 2): last_ov_c])
                else:
                    shredded_seqs[seq_id].append(seqs_dict[seq_id][i + part_l - int(parts_ov / 2):])
            else:
                shredded_seqs[seq_id].append(seqs_dict[seq_id][i:])
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
        return
    if num_words < 1:
        print("Number of words must be at least 1")
        return
    splitters_repl = {"_": " ",
                      "\t": " ",
                      ",": " ",
                      ";": " ",
                      ".": " ",
                      ":": " ",
                      "|": " ",
                      "=": " ",
                      "/": " ",
                      "'": " ",
                      '"': " ",
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
    for seq_name in fasta_dict:
        if len(seq_name) <= num_letters:
            reduced_fasta_dict[seq_name] = fasta_dict[seq_name]
            seq_names_dict[seq_name] = seq_name
        else:
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
                if un_fix[-1] == "N":
                    break
            reduced_fasta_dict[(reduced_seq_name+un_fix).upper()] = fasta_dict[seq_name]
            seq_names_dict[(reduced_seq_name+un_fix).upper()] = seq_name
    return reduced_fasta_dict, seq_names_dict


def _get_part_size_list(num_letters, num_words):
    if num_letters == 6:
        if num_words >= 2:
            return [1, 3]
        if num_words == 1:
            return [4]
    if num_letters == 7:
        if num_words >= 3:
            return [1, 3, 1]
        if num_words == 2:
            return [1, 3]
        if num_words == 1:
            return [4]
    if num_letters == 8:
        if num_words >= 4:
            return [1, 3, 1, 1]
        if num_words == 3:
            return [1, 3, 1]
        if num_words == 2:
            return [1, 4]
        if num_words == 1:
            return [5]
    if num_letters == 9:
        if num_words >= 4:
            return [1, 3, 1, 1]
        if num_words == 3:
            return [1, 4, 1]
        if num_words == 2:
            return [2, 4]
        if num_words == 1:
            return [6]
    if num_letters == 10:
        if num_words >= 4:
            return [1, 4, 1, 1]
        if num_words == 3:
            return [2, 4, 1]
        if num_words == 2:
            return [2, 5]
        if num_words == 1:
            return [7]
    if num_letters >= 11:
        if num_words >= 4:
            return [2, 4, 1, 1]
        if num_words == 3:
            return [2, 5, 1]
        if num_words == 2:
            return [3, 5]
        if num_words == 1:
            return [8]


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
            "seqid": corr_orf_id.split("|:")[0],
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
    # TODO: consider to write this function in parallel mode
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


def detect_seqs_type(seqs_dict=None, nuc_freq_thr=0.75, **kwargs):
    # For backward compatibility
    if kwargs.get("fasta_dict", None) is not None:
        seqs_dict = kwargs["fasta_dict"]
    fasta_path = kwargs.get("fasta_path", None)
    if seqs_dict is None:
        if fasta_path is not None:
            seqs_dict = load_fasta_to_dict(fasta_path)
        else:
            return

    summ_l = 0
    let_counts = Counter()
    for seq_key in seqs_dict:
        seq = None
        seq = seqs_dict[seq_key].lower().replace("-", "")
        let_counts += Counter(seq)
        summ_l += len(seq)
    if float(let_counts.get("a", 0)+let_counts.get("c", 0)+let_counts.get("g", 0)+
                        let_counts.get("t", 0))/float(summ_l) >= nuc_freq_thr:
        return SeqsDict.nucl_type
    else:
        return SeqsDict.prot_type
