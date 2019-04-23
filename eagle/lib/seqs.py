import os
from collections import defaultdict
import subprocess

from Bio.Seq import Seq

from eagle.constants import conf_constants
from eagle.lib.general import filter_list, get_un_fix


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


def load_fasta_to_dict(fasta_path):
    fasta_dict = dict()
    seq_list = list()
    title = None
    fasta_f = open(fasta_path)
    for line_ in fasta_f:
        line = None
        line = line_.strip()
        if not line:
            continue
        if line[0] == ">":
            if title:
                fasta_dict[title] = "".join(seq_list)
                seq_list = list()
                title = None
            title = line[1:]
        else:
            seq_list.append(line)
    if title:
        fasta_dict[title] = "".join(seq_list)
        seq_list = list()
        title = None
    return fasta_dict


def dump_fasta_dict(fasta_dict, fasta_path, overwrite=True):
    if overwrite:
        fasta_f = open(fasta_path, 'w')
    else:
        fasta_f = open(fasta_path, 'a')
    for seq_id in fasta_dict.keys():
        fasta_f.write(">"+seq_id+"\n")
        fasta_f.write(fasta_dict[seq_id]+"\n")
    fasta_f.close()


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
    orfs_ids = orfs_fasta_dict.keys()
    for orf_id in orfs_ids:
        corr_orf_id = None
        ori = None
        corr_orf_id, c_start, c_end, ori = _get_orf_info(orf_id)
        orfs_fasta_dict[corr_orf_id] = orfs_fasta_dict.pop(orf_id)
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
