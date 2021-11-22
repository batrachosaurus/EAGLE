import os
import re
import shutil
from copy import deepcopy
from collections import defaultdict, Counter, OrderedDict
import subprocess

import psutil
import numpy as np
from Bio.Seq import Seq

from jsondler import JsonEntry

from eagle.constants import conf_constants
from eaglib._utils.strings import get_un_fix, generate_random_string
from eaglib._utils.logging import send_log_message


class SeqsDict(object):

    prot_type = "prot"
    nucl_type = "nucl"
    rna_type = "rna"
    prot_types = ("protein", prot_type, "p")
    nucl_types = ("nucleotide", nucl_type, "n", rna_type)

    low_memory0 = False
    _chunk_size0 = 100

    def __init__(self, seqs_order, seqs_array,
                 seq_info_dict=None,
                 seqs_type=None,
                 low_memory=None,
                 **kwargs):

        self._seqs_type = None

        if seq_info_dict is None:
            seq_info_dict = defaultdict(lambda: defaultdict(lambda: None))
        if low_memory is None:
            low_memory = self.low_memory0

        self.seqs_order = seqs_order
        self.seqs_array = seqs_array
        self.seq_info_dict = seq_info_dict
        if seqs_type.lower() in self.prot_types + self.nucl_types:
            self.seqs_type = seqs_type.lower()
        self.low_memory = low_memory
        self._chunk_size = int(str(seqs_array.dtype).strip("<U"))
        self.logger = kwargs.get("logger", None)

        self._nucl_seqs_dict = dict()
        self._full_to_short_seq_names = dict()
        self._short_to_full_seq_names = dict()
        self._seq_ids_to_orgs = dict()

        self._empty_rows = list()  # TODO: implement special object

    def _init_subset(self, seqs_order, seqs_array,
                     seq_info_dict=None,
                     seqs_type=None,
                     low_memory=None,
                     **kwargs):

        if seq_info_dict is None:
            seq_info_dict = defaultdict(lambda: defaultdict(lambda: None))
            for seq_name in seqs_order:
                seq_info_dict[seq_name] = deepcopy(self.seq_info_dict[seq_name])
        if seqs_type is None:
            seqs_type = self.seqs_type
        if low_memory is None:
            low_memory = self.low_memory

        seqs_dict = self.__class__(seqs_order, seqs_array,
                                   seq_info_dict=seq_info_dict,
                                   seqs_type=seqs_type,
                                   low_memory=low_memory,
                                   logger=kwargs.pop("logger", self.logger),
                                   **kwargs)

        prev_end = 0
        for rows_range in seqs_order.values():
            d = rows_range.start - prev_end
            if d > 0:
                seqs_dict._empty_rows.append(range(prev_end, prev_end+d))
            prev_end = rows_range.stop

        seqs_dict.full_to_short_seq_names = self.filter_seq_ids(target_dict=self.full_to_short_seq_names,
                                                                seq_ids_to_remain=seqs_order)

        seqs_dict.seq_ids_to_orgs = self.filter_seq_ids(target_dict=self.seq_ids_to_orgs, seq_ids_to_remain=seqs_order)
        seqs_dict.nucl_seqs_dict = self.filter_seq_ids(target_dict=self.nucl_seqs_dict, seq_ids_to_remain=seqs_order)
        return seqs_dict

    @staticmethod
    def filter_seq_ids(target_dict, seq_ids_to_remain):
        return dict(filter(lambda kv: kv[0] in seq_ids_to_remain, target_dict.items()))

    @property
    def seqs_type(self):
        if self._seqs_type is not None:
            return self._seqs_type
        else:
            return self.detect_seqs_type()

    @seqs_type.setter
    def seqs_type(self, seqs_type):
        self._seqs_type = seqs_type

    @property
    def nucl_seqs_dict(self):
        return self._nucl_seqs_dict

    @nucl_seqs_dict.setter
    def nucl_seqs_dict(self, nucl_seqs_dict):
        if self.seqs_type.lower() not in self.prot_types:
            send_log_message(message="The alignment seems to be nucleotide - not applicable",
                             mes_type="w", logger=self.logger)
        else:
            if isinstance(nucl_seqs_dict, dict):
                nucl_seqs_dict = SeqsDict.load_from_dict(nucl_seqs_dict)
            assert isinstance(nucl_seqs_dict, SeqsDict), "ERROR: the value for argument 'nucl_seqs_dict' should be " \
                                                         "an instance of eagle.eaglib.seqs.SeqsDict"
            self._nucl_seqs_dict = nucl_seqs_dict

    @property
    def full_to_short_seq_names(self):
        if not self._full_to_short_seq_names:
            mult_aln_dict_short_ids, self.short_to_full_seq_names = reduce_seq_names(fasta_dict=self, num_letters=10)
        return self._full_to_short_seq_names

    @full_to_short_seq_names.setter
    def full_to_short_seq_names(self, full2short_dict):
        if isinstance(full2short_dict, dict):  # not exhaustive condition
            self._full_to_short_seq_names = full2short_dict
            self._short_to_full_seq_names = {v: k for k, v in full2short_dict}
        else:
            send_log_message(message="the value to assign should be a dict", mes_type="e", logger=self.logger)

    @property
    def short_to_full_seq_names(self):
        if not self._short_to_full_seq_names:
            mult_aln_dict_short_ids, self.short_to_full_seq_names = reduce_seq_names(fasta_dict=self, num_letters=10)
        return self._short_to_full_seq_names

    @short_to_full_seq_names.setter
    def short_to_full_seq_names(self, short2full_dict):
        if isinstance(short2full_dict, dict):  # not exhaustive condition
            self._full_to_short_seq_names = {v: k for k, v in short2full_dict}
            self._short_to_full_seq_names = short2full_dict
        else:
            send_log_message(message="the value to assign should be a dict", mes_type="e", logger=self.logger)

    @property
    def seq_ids_to_orgs(self):
        return self._seq_ids_to_orgs

    @seq_ids_to_orgs.setter
    def seq_ids_to_orgs(self, sito_dict):
        if isinstance(sito_dict, dict):  # not exhaustive condition
            self._seq_ids_to_orgs = sito_dict
        else:
            send_log_message(message="the value to assign should be a dict", mes_type="e", logger=self.logger)

    @property
    def short_id(self):
        if self.seqs_order:
            self_short_ids = deepcopy(self)
            self_short_ids.rename_seqs(self.full_to_short_seq_names)
            return self_short_ids
        else:
            return dict()

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
        return self._init_subset(seqs_order, seqs_array, low_memory=low_memory)

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
        # TODO: implement real delete with elements shift up
        self._empty_rows.append(self.seqs_order.pop(key))

    def __del__(self):
        if self.low_memory:
            dat_path = self.seqs_array.filename
            del self.seqs_array
            os.remove(dat_path)

    def __contains__(self, item):
        return item in self.seqs_order

    def keys(self):
        return map(lambda si: si[0], sorted(self.seqs_order.items(), key=lambda si: si[1].start))

    def values(self):
        for seq_name in self:
            yield self[seq_name]

    def items(self):
        return map(lambda seq_id: (seq_id, self[seq_id]), self.keys())

    @property
    def seq_names(self):
        if self.seqs_order:
            return list(self.keys())
        else:
            return list()

    def rename_seqs(self, old_to_new_dict):
        for old_seq in old_to_new_dict:
            if old_seq in self.seqs_order:
                self.seqs_order[old_to_new_dict[old_seq]] = self.seqs_order.pop(old_seq)

        for short_name in list(self.short_to_full_seq_names.keys()):  # TODO: rename seq in all 'dicts'
            self.short_to_full_seq_names[short_name] = old_to_new_dict.get(self.short_to_full_seq_names[short_name],
                                                                           self.short_to_full_seq_names[short_name])

    @property
    def num_seqs(self):
        return len(self)

    def copy(self):
        return self._init_subset(self.seqs_order, self.seqs_array)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self._init_subset(deepcopy(self.seqs_order), deepcopy(self.seqs_array))

    @classmethod
    def load_from_file(cls, fname=None, format="fasta", seqs_type=None, low_memory='auto', **kwargs):
        if "seqs_path" in kwargs:
            fname = None
            fname = kwargs["seqs_path"]
        if "seqs_format" in kwargs:
            format = None
            format = kwargs["seqs_format"]
        assert fname is not None, "ERROR: no value passed for argument 'fname'"

        chunk_size = kwargs.get("chunk_size", cls._chunk_size0)
        if low_memory == 'auto' or low_memory is None:
            low_memory = cls._check_low_memory(seqs_path=fname)

        n_chunks = 0
        if format.lower() in ("fasta", "fas", "fa"):
            with open(fname) as fasta_f:
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
                if seq_l > 0:
                    n_chunks += (seq_l - 1) // chunk_size + 1
                    seq_l = 0
        if format.lower() in ("phylip", "phy", "ph"):
            print("not implemented yet")  # TODO: implement phylip format reading

        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=np.dtype("U%s" % chunk_size), mode='w+', shape=(n_chunks,))
        else:
            seqs_array = np.zeros(n_chunks, dtype=np.dtype("U%s" % chunk_size))

        i = 0
        seqs_order = dict()
        if format.lower() in ("fasta", "fas", "fa"):
            with open(fname) as fasta_f:
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
        return cls(seqs_order, seqs_array,
                   seqs_type=seqs_type,
                   low_memory=low_memory,
                   **kwargs)

    @classmethod
    def load_from_dict(cls, in_dict, seqs_type=None, low_memory=None, **kwargs):
        chunk_size = kwargs.get("chunk_size", cls._chunk_size0)
        n_chunks = 0
        for seq in in_dict.values():
            if seq:
                n_chunks += (len(seq)-1)//chunk_size + 1

        if low_memory is None:
            low_memory = cls.low_memory0
        if low_memory:
            dat_path = kwargs.get("dat_path", "." + generate_random_string(10) + "_seqs_dict.dat")
            seqs_array = np.memmap(dat_path, dtype=np.dtype("U%s" % chunk_size), mode='w+', shape=(n_chunks,))
        else:
            seqs_array = np.zeros(n_chunks, dtype=np.dtype("U%s" % chunk_size))

        i = 0
        seqs_order = dict()
        for seq_id, seq in in_dict.items():
            seqs_order[seq_id] = range(i, i + (len(seq)-1) // chunk_size + 1)
            for n, j in enumerate(seqs_order[seq_id]):
                seqs_array[j] = seq[n * chunk_size: (n + 1) * chunk_size]
            i = seqs_order[seq_id].stop
        return cls(seqs_order, seqs_array,
                   seqs_type=seqs_type,
                   low_memory=low_memory,
                   **kwargs)

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

    def dump(self, fname=None, format="fasta", overwrite=True, **kwargs):
        if "seqs_path" in kwargs:
            fname = None
            fname = kwargs["seqs_path"]
        if "seqs_format" in kwargs:
            format = None
            format = kwargs["seqs_format"]
        assert fname is not None, "ERROR: no value passed for argument 'fname'"

        if format.lower() == "fasta":
            if overwrite:
                fasta_f = open(fname, 'w')
            else:
                fasta_f = open(fname, 'a')
            if kwargs.get("replace_stops", False):
                for seq_id in self.seqs_order:
                    fasta_f.write(">" + seq_id + "\n")
                    fasta_f.write(self[seq_id].replace("*", "X") + "\n")
            else:
                for seq_id in self.seqs_order:
                    fasta_f.write(">" + seq_id + "\n")
                    fasta_f.write(self[seq_id] + "\n")
            fasta_f.close()

        if format.lower() == "phylip":
            self_short_id = self.short_id
            with open(fname, "w") as phylip_f:
                phylip_f.write("    %s    %s\n" %
                               kwargs.get("shape", (len(self), 100 * len(self.seqs_array) / len(self))))
                for seq_name in self_short_id:
                    num_spaces_to_add = 10 - len(seq_name)
                    spaces_to_add = [" " for i in range(num_spaces_to_add)]
                    phylip_f.write("%s %s\n" % (seq_name + "".join(spaces_to_add), self_short_id[seq_name]))

        return fname

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

    def generate_seqs_pairs(self, only_first_seq=False, raref_base=None):
        """
        Generates pars of aligned sequences from multiple alignment (n(n-1)/2 by deafault)
        :param only_first_seq: set it True to get only pairs containig first sequence (default False)
        :param raref_base: randomly rarefies pairs: n(n-1)/2 -> n*raref_base/2 (10.0 makes sense to use)
        :return:
        """
        seqs_pairs = dict()
        seq_names_list = list(self.seq_names)
        for i, seqi_name in enumerate(seq_names_list[:-1]):
            for seqj_name in seq_names_list[i + 1:]:
                seqs_pairs[frozenset({seqi_name, seqj_name})] = [self[seqi_name], self[seqj_name]]
            if only_first_seq:
                break
        if raref_base is not None and int(raref_base) < len(seq_names_list) and not only_first_seq:
            for seqs_pair in np.random.permutation(list(seqs_pairs.keys()))[
                             (len(seq_names_list) - 1) * int(raref_base / 2.0):]:
                del seqs_pairs[seqs_pair]
        return seqs_pairs

    def stop_codons_stats(self, **kwargs):
        stops_per_seq = list()
        if self.seqs_type in self.prot_types:
            for seq_name in self:
                stops_per_seq.append(self[seq_name].count("*"))
        else:
            for seq_name in self:
                stops_per_seq.append(re.sub("(tag|taa|tga)", "*", self[seq_name].lower()).count("*"))
        stops_per_seq.sort()
        return {"stops_per_seq_median": np.median(stops_per_seq),
                "seqs_with_stops_fract": float(len(list(filter(lambda x: x > 0, stops_per_seq)))) / float(len(self))}

    def rarefy(self, seqs_to_remain=100, **kwargs):
        seqs_ids = self.seq_names
        if len(seqs_ids) <= seqs_to_remain:
            return self
        rarefied_aln_dict = dict()
        for i in range(seqs_to_remain):
            seq_id = None
            seq_id = seqs_ids.pop(np.random.randint(len(seqs_ids)))
            rarefied_aln_dict[seq_id] = self[seq_id]
        rarefied_aln_seqs = SeqsDict.load_from_dict(rarefied_aln_dict)
        return self._init_subset(rarefied_aln_seqs.seqs_order, rarefied_aln_seqs.seqs_array, **kwargs)

    def construct_mult_aln(self,
                           method="MUSCLE",
                           muscle_exec_path=None,
                           mafft_exec_path=None,
                           msaprobs_exec_path=None,
                           emboss_inst_dir=None,
                           hmmer_inst_dir=None,
                           infernal_inst_dir=None,
                           fastme_exec_path=None,
                           kaks_calculator_exec_path=None,
                           aln_name="mult_aln",
                           tmp_dir=None,
                           remove_tmp=True,
                           num_threads=None,
                           **kwargs):

        from eaglib.alignment import MultAln

        if muscle_exec_path is None:
            muscle_exec_path = conf_constants.muscle_exec_path
        if mafft_exec_path is None:
            mafft_exec_path = conf_constants.mafft_exec_path
        if msaprobs_exec_path is None:
            msaprobs_exec_path = conf_constants.msaprobs_exec_path
        if emboss_inst_dir is None:
            emboss_inst_dir = conf_constants.emboss_inst_dir
        if hmmer_inst_dir is None:
            hmmer_inst_dir = conf_constants.hmmer_inst_dir
        if infernal_inst_dir is None:
            infernal_inst_dir = conf_constants.infernal_inst_dir
        if fastme_exec_path is None:
            fastme_exec_path = conf_constants.fastme_exec_path
        if kaks_calculator_exec_path is None:
            kaks_calculator_exec_path = conf_constants.fastme_exec_path
        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_mult_aln_tmp"
        if num_threads is None:
            num_threads = conf_constants.num_threads

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        in_fasta_path = os.path.join(tmp_dir, "to_align.fasta")
        out_fasta_path = os.path.join(tmp_dir, aln_name + ".fasta")

        if method.lower() == "muscle":
            send_log_message("MUSCLE is starting", mes_type="info", logger=self.logger)
            muscle_cmd = muscle_exec_path + " -in " + in_fasta_path + " -out " + out_fasta_path
            subprocess.call(muscle_cmd, shell=True)
            send_log_message("MUSCLE finished", mes_type="info", logger=self.logger)
            mult_aln = MultAln.load_from_file(out_fasta_path, format="fasta", restore_stops=True, **kwargs)

        if method.lower() == "mafft":
            send_log_message("MAFFT is starting", mes_type="info", logger=self.logger)
            mafft_cmd = mafft_exec_path + " --auto" \
                        " --op " + str(kwargs.get("op", kwargs.get("gap_open_penalty", 1.53))) + \
                        " --ep " + str(kwargs.get("ep", kwargs.get("gap_ext_penalty", 0.123))) + \
                        " --thread " + str(num_threads) + \
                        " " + in_fasta_path + " > " + out_fasta_path
            subprocess.call(mafft_cmd, shell=True)
            send_log_message("MAFFT finished", mes_type="info", logger=self.logger)
            mult_aln = MultAln.load_from_file(out_fasta_path, format="fasta", restore_stops=True, **kwargs)

        if method.lower() == "msaprobs":
            send_log_message("MSAProbs is starting", mes_type="info", logger=self.logger)
            msaprobs_cmd = msaprobs_exec_path + " -num_threads " + str(num_threads) + " -v " + \
                           in_fasta_path + " > " + out_fasta_path
            subprocess.call(msaprobs_cmd, shell=True)
            send_log_message("MSAProbs finished", mes_type="info", logger=self.logger)
            mult_aln = MultAln.load_from_file(out_fasta_path, format="fasta", restore_stops=True, **kwargs)

        mult_aln.seqs_type = self.seqs_type
        mult_aln.emboss_inst_dir = emboss_inst_dir
        mult_aln.hmmer_inst_dir = hmmer_inst_dir
        mult_aln.infernal_inst_dir = infernal_inst_dir
        mult_aln.fastme_exec_path = fastme_exec_path
        mult_aln.kaks_calculator_exec_path = kaks_calculator_exec_path
        mult_aln.name = aln_name
        mult_aln.storage_dir = tmp_dir

        if remove_tmp:
            shutil.rmtree(tmp_dir)
        return mult_aln

    def to_blastdb(self, dbtype=None, db_name=None, blast_inst_dir=None, **kwargs):
        from eaglib.alignment import BlastDB

        return BlastDB.make_blastdb(db_name=db_name,
                                    dbtype=dbtype,
                                    blast_inst_dir=blast_inst_dir,
                                    tmp_dir=kwargs.get("blast_tmp_dir", None),
                                    logger=self.logger)

    def to_blast_search(self, blast_type, blast_db, out=None, num_threads=1, outfmt=7, max_hsps=100,
                        **kwargs):
        from eaglib.alignment import BlastDB

        assert isinstance(blast_db, BlastDB), "ERROR: the value for argument 'blast_db' should be " \
                                              "an instance of class eagle.eaglib.alignment.BlastDB"
        return blast_db.run_blast_search(blast_type=blast_type,
                                         query=self,
                                         out=out,
                                         num_threads=num_threads,
                                         outfmt=outfmt,
                                         max_hsps=max_hsps,
                                         **kwargs)


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


def shred_seqs(seqs, part_l=50000, parts_ov=5000, shredded_seqs_fasta=None):
    if isinstance(seqs, str) and os.path.exists(seqs):
        seqs_dict = SeqsDict.load_from_file(seqs, format="fasta")
    else:
        seqs_dict = seqs

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

    shredded_seqs_dict = OrderedDict()
    if shredded_seqs_fasta is not None:
        for seq_id in shredded_seqs:
            i = 0
            for seq in shredded_seqs[seq_id]:
                start = i * (part_l - parts_ov)
                shredded_seqs_dict[seq_id + "|:" + str(start + 1) + "-" + str(start + part_l)] = seq
                i += 1
        SeqsDict.load_from_dict(shredded_seqs_dict).dump(shredded_seqs_fasta, format="fasta")
        return shredded_seqs_fasta
    else:
        return shredded_seqs


def load_fasta_to_dict(fasta_path, low_memory="auto", **kwargs):
    return SeqsDict.load_from_file(fname=fasta_path, format="fasta", low_memory=low_memory, **kwargs)


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
            seq_name_list = list(filter(lambda li: li.strip(),
                                        "".join([splitters_repl.get(s, s) for s in seq_name]).split()))
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


def nucl_accord_prot(prot_seq, nucl_seq):
    nucl_seq_list = list()
    i = 0
    for aa in prot_seq:
        if aa == "-":
            nucl_seq_list.append("---")
        elif aa == ".":
            nucl_seq_list.append("...")
        else:
            nucl_seq_list.append(nucl_seq[i*3:(i+1)*3])
            i += 1
    return "".join(nucl_seq_list)


class GenomeInfo(JsonEntry):

    # json keys
    genome_id_key = "genome_id"
    org_name_key = "org_name"
    taxonomy_key = "taxonomy"
    fna_seqs_key = "fna_seqs"
    fna_seq_fasta_key = "fasta"
    fna_seq_ids_key = "ids"
    repr_seqs_key = "repr_seqs"  # repr - for example, btc or btr
    repr_seq_fasta_key = "fasta"
    repr_seq_id2profile_key = "id2profile"
    extra_info_key = "extra_info"

    # default values
    genome_id_0 = None
    org_name_0 = None
    taxonomy_0 = list()
    fna_seq_fasta_0 = list()
    fna_seq_ids_0 = list()
    repr_seq_fasta_0 = list()
    repr_seq_id2profile_0 = defaultdict(str)  # {repr_seq_id: repr_seq_profile_name}
    extra_info_0 = dict()

    def __init__(self,
                 genome_id=genome_id_0,
                 org_name=org_name_0,
                 taxonomy=None,
                 fna_seq_fasta=None,
                 fna_seq_ids=None,
                 repr_seq_fasta=None,
                 repr_seq_id2profile=None,
                 extra_info=None):

        # attribute names must match keys form GenomeInfo.attr_scheme()
        if taxonomy is None:
            taxonomy = self.taxonomy_0
        if fna_seq_fasta is None:
            fna_seq_fasta = self.fna_seq_fasta_0
        if fna_seq_ids is None:
            fna_seq_ids = self.fna_seq_ids_0
        if repr_seq_fasta is None:
            repr_seq_fasta = self.repr_seq_fasta_0
        if repr_seq_id2profile is None:
            repr_seq_id2profile = self.repr_seq_id2profile_0
        if extra_info is None:
            extra_info = self.extra_info_0
        
        self.genome_id = genome_id
        self.org_name = org_name
        self.taxonomy = taxonomy
        self.fna_seq_fasta = fna_seq_fasta  # list of paths or links
        self.fna_seq_ids = fna_seq_ids
        self.repr_seq_fasta = repr_seq_fasta  # list of paths
        self.repr_seq_id2profile = repr_seq_id2profile
        self.extra_info = extra_info

    @classmethod
    def attr_scheme(cls):
        # json scheme (the keys must match attribute names defined in __init__)
        # CAUTION: if attribute is not defined in __init__ method (sublevel key in a json)
        # don't include it to the result dict
        return {
            "genome_id": (cls.genome_id_key,),
            "org_name": (cls.org_name_key,),
            "taxonomy": (cls.taxonomy_key,),
            "fna_seq_fasta": (cls.fna_seqs_key, cls.fna_seq_fasta_key),
            "fna_seq_ids": (cls.fna_seqs_key, cls.fna_seq_ids_key),
            "repr_seq_fasta": (cls.repr_seqs_key, cls.repr_seq_fasta_key,),
            "repr_seq_id2profile": (cls.repr_seqs_key, cls.repr_seq_id2profile_key,),
            "extra_info": (cls.extra_info_key,),
        }

    @classmethod
    def org_name_from_dict(cls, in_dict):
        return in_dict["org_name"]
