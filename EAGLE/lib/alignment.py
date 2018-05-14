import collections
import os
import shutil
import subprocess
from collections import defaultdict
import pandas

from EAGLE.lib.general import load_fasta_to_dict, dump_fasta_dict, load_phylip_dist_matrix, reduce_seq_names, ConfBase


class MultAln(ConfBase):

    def __init__(self,
                 mult_aln_dict=None,
                 aln_type=None,
                 states_seq=None,
                 tmp_dir="tmp",
                 emboss_inst_dir="",
                 hmmer_inst_dir="",
                 config_path=None,
                 logger=None):

        if mult_aln_dict:
            self.mult_aln_dict = mult_aln_dict
            self.mult_aln_dict_short_id, self.short_to_full_seq_names = reduce_seq_names(fasta_dict=mult_aln_dict,
                                                                                num_letters=10)
            self.full_to_short_seq_names = dict(map(lambda x: (x[1], x[0]), self.short_to_full_seq_names.items()))
        else:
            self.mult_aln_dict = None
            self.mult_aln_dict_short_id = None
            self.short_to_full_seq_names = None
            self.full_to_short_seq_names = None
        self.states_seq = states_seq
        self.aln_type = aln_type
        if not self.aln_type:
            self.aln_type = detect_seqs_type(self.mult_aln_dict)
        self.distance_matrix = None
        self.tmp_dir = tmp_dir
        self.emboss_inst_dir = emboss_inst_dir
        self.hmmer_inst_dir = hmmer_inst_dir
        self.logger = logger
        super(MultAln, self).__init__(config_path=config_path)

    def __getitem__(self, seq_id):
        return self.mult_aln_dict[seq_id]

    def __setitem__(self, seq_id, seq):
        self.mult_aln_dict[seq_id] = seq

    def seqs(self):
        if self.mult_aln_dict:
            return [seq_id for seq_id in self.mult_aln_dict.keys()]
        else:
            return list()

    def dump_alignment(self, aln_f_path):
        pass

    def improve_aln(self,
                    max_gap_fract=0.95,  # maximal fraction of gaps in a column to keep it in alignment
                    max_mismatch_fract=1.0,  # maximal fraction of mismatches in a column to keep it in alignment
                    remove_seq=False,  # if True can remove sequences for cleaning gaps between aln blocks
                    dist_filt=False,  # if True returns a list of alignments (or blocks coordinates) with different strictness of sequences removing (different thresholds for clustering)
                    output='alignment',
                    inplace=False):

        if dist_filt:
            # Constructs a list of rarefied by different distance thresholds alignments and runs improve_aln on it with dist_filt=False
            pass
        coords = list()
        seq_id_list = self.mult_aln_dict.keys()
        for i in range(len(self.mult_aln_dict[seq_id_list[0]])):
            if self.states_seq:
                if self.states_seq[i] in ("match", "specific"):
                    coords = self._update_coords(i+1, coords)
                    continue
            let_counts = defaultdict(int)
            for seq_id in seq_id_list:
                let_counts[self.mult_aln_dict[seq_id][i]] += 1
            if float(let_counts.pop("-"))/float(len(seq_id_list)) <= max_gap_fract:
                if sorted(let_counts.values(), reverse=True)[0] >= 1.0-max_mismatch_fract:
                    coords = self._update_coords(i+1, coords)
            elif remove_seq:
                # TODO: write detection of seqs to remove
                pass
        if len(coords[-1]) == 1:
            coords.pop()
        if not remove_seq and len(coords) >= 1:
            coords = (coords[0][0], coords[-1][1])
        if output.lower() in ("coordinates", "coords", "coord", "c"):
            seq_coord_list = list()
            for seq_id in self.mult_aln_dict.iterkeys():
                seq_c_dict = {"seq_id": seq_id}
                coords_for_seq = self._get_coords_for_seq(coords, self.mult_aln_dict[seq_id])
                for k in range(len(coords_for_seq)):
                    seq_c_dict["c%s" % ((k+1)*2-1)] = coords_for_seq[k][0]
                    seq_c_dict["c%s" % (k+1)*2] = coords_for_seq[k][1]
                seq_coord_list.append(seq_c_dict)
            return pandas.DataFrame(seq_coord_list)
        else:
            pass  # TODO: write output as alignment

    @staticmethod
    def _update_coords(i, coords):
        if not coords:
            coords.append([i])
        elif len(coords[-1]) == 1:
            coords[-1].append(i)
        elif coords[-1][1] == i - 1:
            coords[-1][1] = i
        else:
            coords.append([i])
        return coords

    @staticmethod
    def _get_coords_for_seq(coords, gapped_seq):
        no_gaps_coords = list()
        for coord_pair in coords:
            c1_shift = gapped_seq[:coord_pair[0]].count("-")
            c2_shift = gapped_seq[:coord_pair[1]].count("-")
            no_gaps_coords.append((coord_pair[0]-c1_shift, coord_pair[1]-c2_shift))
        return no_gaps_coords

    def get_distance_matrix(self, method="phylip"):
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if method.lower() == "phylip":
            aln_fasta_path = os.path.join(self.tmp_dir, "alignment.fasta")
            phylip_matrix_path = os.path.join(self.tmp_dir, "dist_matrix.phylip")
            if self.mult_aln_dict:
                if not self.mult_aln_dict_short_id:
                    self.mult_aln_dict_short_id, self.short_to_full_seq_names = \
                        reduce_seq_names(fasta_dict=self.mult_aln_dict, num_letters=10)
                    self.full_to_short_seq_names = dict(map(lambda x: (x[1], x[0]),
                                                            self.short_to_full_seq_names.items()))
            else:
                if self.logger:
                    self.logger.warning("No sequences in alignment")
                else:
                    print("No sequences in alignment")
                return 1
            dump_fasta_dict(fasta_dict=self.mult_aln_dict_short_id, fasta_path=aln_fasta_path)
            if self.aln_type.lower() in ("protein", "prot", "p"):
                phylip_cmd = os.path.join(self.emboss_inst_dir, "fprotdist") + " -sequence " + aln_fasta_path + \
                             " -outfile " + phylip_matrix_path
            else:
                phylip_cmd = os.path.join(self.emboss_inst_dir, "fdnadist") + " -sequence " + aln_fasta_path + \
                             " -method f -outfile " + phylip_matrix_path
            subprocess.call(phylip_cmd, shell=True)
            self.distance_matrix = load_phylip_dist_matrix(matrix_path=phylip_matrix_path)
        shutil.rmtree(self.tmp_dir)
        return self.distance_matrix

    def remove_paralogs(self, seq_ids_to_orgs, method="min_dist", inplace=False):  ###
        short_seq_ids_to_org = self._check_seq_ids_to_org(seq_ids_to_orgs)
        if not self.distance_matrix:
            self.get_distance_matrix()
        org_dist_dict = defaultdict(dict)
        seq_ids = self.distance_matrix.index
        for seq_id in seq_ids:
            org_dist_dict[seq_ids_to_orgs[seq_id]][seq_id] = sum(map(float, list(self.distance_matrix.loc(seq_id))))
        for org in org_dist_dict.keys():
            org_dist_dict[org] = sorted(org_dist_dict[org].items(), key=lambda x: x[1])
        if method.lower() in ("minimal_distance", "min_dist", "md"):
            full_seq_names_filt = dict(filter(lambda x: x[0] == org_dist_dict.get(short_seq_ids_to_org.get(x[0], None),
                                                                                  ((None, None), None))[0][0],
                                              self.short_to_full_seq_names.items()))
            mult_aln_dict_filt = dict(filter(lambda x: x[0] in full_seq_names_filt.keys(), self.mult_aln_dict.items()))
        else:
            # TODO: write spec_pos method
            full_seq_names_filt = None
            mult_aln_dict_filt = None
        if inplace:
            self.mult_aln_dict = mult_aln_dict_filt
            self.short_to_full_seq_names = full_seq_names_filt
        else:
            filtered_aln = self
            filtered_aln.mult_aln_dict = mult_aln_dict_filt
            filtered_aln.short_to_full_seq_names = full_seq_names_filt
            return filtered_aln

    def _check_seq_ids_to_org(self, seq_ids_to_orgs):
        to_short_dict = dict()
        long_suc_num = 0
        short_suc_num = 0
        for key in self.short_to_full_seq_names.keys():
            if seq_ids_to_orgs.get(key, None):
                short_suc_num += 1
            if seq_ids_to_orgs.get(self.short_to_full_seq_names[key], None):
                to_short_dict[self.short_to_full_seq_names[key]] = key
                long_suc_num += 1
        if short_suc_num >= long_suc_num:
            return seq_ids_to_orgs
        else:
            return dict(map(lambda x: (to_short_dict.get(x[0], None),  x[1]), seq_ids_to_orgs.items()))

    def get_blocks_tsv(self, tsv_path, fasta_path, split_into_blocks=False, meta_dict=None):
        # Three block types: CB - conservetive block, SB - specific block, MB - major block. UB - unclassified block (not any of three types => MB=NA)
        cut_ends_df = self.improve_aln(output='coords')
        if split_into_blocks:
            blocks_df = self.define_aln_blocks(cut_ends_coord_dict=cut_ends_df)
        else:
            blocks_df = cut_ends_df
            blocks_df["block_type"] = pandas.Series(["UB"]*cut_ends_df.shape[0])
            blocks_df["block_descr"] = pandas.Series(['major_block_id "NA"; block_id "UB1"; block_name "NA"']*
                                                     cut_ends_df.shape[0])
            blocks_df = blocks_df[["seq_id", "block_type", "start", "end", "block_descr"]]
        tsv_f = open(tsv_path, 'w')
        seqs_dict = dict(list(blocks_df.apply(self._write_blocks_tsv, axis=1, args=(tsv_f, meta_dict))))
        tsv_f.close()
        dump_fasta_dict(fasta_dict=seqs_dict, fasta_path=fasta_path)

    def _write_blocks_tsv(self, block_info, tsv_f, meta_dict=None):
        tsv_f.write("\t".join(list(block_info))+"; %s\n" %
                    "; ".join('%s "%s"' % (meta[0], meta[1]) for meta in meta_dict[block_info["seq_id"]].items()))
        return block_info["seq_id"], self.mult_aln_dict[block_info["seq_id"]].replace("-", "")

    def remove_seqs(self, seqs_list):
        pass

    def define_aln_blocks(self, cut_ends_coord_dict):
        blocks_info = list()  # list of dicts
        pass
        return pandas.DataFrame(blocks_info)

    def get_hmm_profile(self, profile_path, method="hmmer"):
        if method.lower() == "hmmer":
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)
            aln_fasta_path = os.path.join(self.tmp_dir, "aln.fasta")
            dump_fasta_dict(self.mult_aln_dict, aln_fasta_path)
            hmmer_handler = HmmerHandler(inst_dir=self.hmmer_inst_dir, config_path=self.config_path, logger=self.logger)
            hmmer_handler.build_hmm_profile()
            shutil.rmtree(self.tmp_dir)

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


class BlastHandler(ConfBase):

    def __init__(self,
                 inst_dir="",
                 config_path=None,
                 logger=None):

        self.inst_dir = inst_dir
        self.logger = logger

        super(BlastHandler, self).__init__(config_path=config_path)

    def make_blastdb(self, in_fasta, dbtype, db_name=None):
        if db_name:
            makeblastdb_cmd = os.path.join(self.inst_dir, "makeblastdb") + " -in " + in_fasta + " -dbtype " + dbtype + \
                              " -out " + db_name
        else:
            makeblastdb_cmd = os.path.join(self.inst_dir, "makeblastdb") + " -in " + in_fasta + " -dbtype " + dbtype
        subprocess.call(makeblastdb_cmd, shell=True)

    def run_blast_search(self, blast_type, query, db, out, outfmt=7, max_hsps=100):
        pass


class HmmerHandler(ConfBase):

    def __init__(self,
                 inst_dir="",
                 config_path=None,
                 logger=None):

        self.inst_dir = inst_dir
        self.logger = logger

        super(HmmerHandler, self).__init__(config_path=config_path)

    def build_hmm_profile(self, profile_path, in_fasta_path):
        hmmbuild_cmd = os.path.join(self.inst_dir, "hmmbuild") + " " + profile_path + " " + in_fasta_path
        subprocess.call(hmmbuild_cmd, shell=True)

    def run_hmmscan(self):
        pass


def construct_mult_aln(seq_dict=None,
                       fasta_path=None,
                       method="MUSCLE",
                       aln_type=None,
                       muscle_exec_path="",
                       emboss_inst_dir="",
                       hmmer_inst_dir="",
                       tmp_dir="tmp",
                       remove_tmp=True,
                       config_path=None,
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
        muscle_cmd = muscle_exec_path + " -in " + fasta_path + " -out " + out_fasta_path
        subprocess.call(muscle_cmd, shell=True)

    mult_aln_dict = load_fasta_to_dict(out_fasta_path)
    if remove_tmp:
        shutil.rmtree(tmp_dir)
    return MultAln(mult_aln_dict=mult_aln_dict,
                   aln_type=aln_type,
                   emboss_inst_dir=emboss_inst_dir,
                   hmmer_inst_dir=hmmer_inst_dir,
                   tmp_dir=tmp_dir,
                   config_path=config_path,
                   logger=logger)


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
