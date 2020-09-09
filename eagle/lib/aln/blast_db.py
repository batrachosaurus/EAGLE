import os
import shutil
import subprocess
import multiprocessing as mp

from Bio.Blast import NCBIXML

from eagle.constants import conf_constants
from eagle.lib.general import send_log_message, generate_random_string, join_files
from eagle.lib.seqs import SeqsDict


class BlastDB(object):

    def __init__(self,
                 dbtype=None,
                 db_name=None,
                 blast_inst_dir=None,
                 tmp_dir=None,
                 logger=None,
                 **kwargs):

        if tmp_dir is None:
            tmp_dir = generate_random_string(10) + "_blast_tmp"

        self.dbtype=dbtype
        self.db_name = db_name
        self.blast_inst_dir = blast_inst_dir
        if self.blast_inst_dir is None:
            self.blast_inst_dir = conf_constants.blast_inst_dir
        self.tmp_dir = tmp_dir
        self.logger = logger

    @classmethod
    def make_blastdb(cls, in_seqs=None, dbtype=None, db_name=None, blast_inst_dir=None, **kwargs):
        if "in_fasta" in kwargs:
            in_seqs = kwargs["in_fasta"]
        assert in_seqs is not None, "ERROR: no value passed for argument 'in_seqs'"

        blast_db = cls(dbtype=dbtype, db_name=db_name, blast_inst_dir=blast_inst_dir, **kwargs)

        if isinstance(in_seqs, SeqsDict):
            assert blast_db.db_name is not None, "ERROR: no value for argument 'db_name'"
            in_fasta = in_seqs.dump(os.path.splitext(blast_db.db_name)[0]+".fasta", format="fasta")
        else:
            in_fasta = in_seqs
            if blast_db.db_name is None and type(in_fasta) is str and os.path.exists(in_fasta):
                blast_db.db_name = os.path.splitext(in_fasta)[0]

        makeblastdb_cmd = os.path.join(blast_db.blast_inst_dir, "makeblastdb") + " -in " + in_fasta + \
                          " -dbtype " + blast_db.dbtype + \
                          " -out " + blast_db.db_name
        subprocess.call(makeblastdb_cmd, shell=True)
        return blast_db

    def run_blast_search(self, blast_type, query, db=None, out=None, num_threads=1, outfmt=7, max_hsps=100, **kwargs):
        if db is None:
            db = self.db_name
        read_output = False
        if out is None:
            out = os.path.splitext(self.db_name)[0] + "_%s.xml" % generate_random_string(10)
            outfmt = 5
            read_output = True
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        if num_threads > 1 and kwargs.get("split_input", True):
            send_log_message("splitting '%s' into %s parts" % (query, num_threads), mes_type="info", logger=self.logger)
            if self.logger is not None:
                self.logger.info()
            else:
                print("INFO: splitting '%s' into %s parts" % (query, num_threads))
            if isinstance(query, SeqsDict):
                query_dict = query
                query_path = os.path.join(self.tmp_dir, "query.fasta")
            else:
                query_dict = SeqsDict.load_from_file(query, format="fasta",
                                                     dat_path=os.path.join(self.tmp_dir, ".%s.dat" % query))
                query_path = query
            query_chunk_size = len(query_dict) // num_threads + 1
            p_list = list()
            query_seqs = query_dict.seq_names
            chunk_kwargs = kwargs.copy()
            chunk_kwargs["split_input"] = False
            chunk_kwargs["remove_tmp"] = False
            query_chunks_list = list()
            for i in range(num_threads):
                query_chunk_path = None
                query_chunk_path = os.path.join(self.tmp_dir,
                                                ("_%s" % i).join(os.path.splitext(os.path.basename(query_path))))
                query_chunks_list.append([query_chunk_path, query_chunk_path + ".bl"])
                if len(query_dict) == i*query_chunk_size:
                    del query_chunks_list[-1]
                    continue
                elif (i+1) * query_chunk_size > len(query_dict):
                    query_dict.get_sample(query_seqs[i * query_chunk_size:]).dump(seqs_path=query_chunks_list[-1][0])
                else:
                    query_dict.get_sample(query_seqs[i*query_chunk_size:
                                                     (i+1)*query_chunk_size]).dump(seqs_path=query_chunks_list[-1][0])
                p = mp.Process(target=self.run_blast_search,
                               args=(blast_type, query_chunks_list[-1][0], db,
                                     query_chunks_list[-1][1], 1, outfmt, max_hsps),
                               kwargs=chunk_kwargs)
                p.start()
                p_list.append(p)
            for p in p_list:
                p.join()
            join_files(in_files_list=list(map(lambda p: p[1], query_chunks_list)), out_file_path=out)

        else:
            if isinstance(query, SeqsDict):
                query_path = os.path.join(self.tmp_dir, "query.fasta")
                query.dump(query_path, format="fasta")
            else:
                query_path = query
            blast_search_cmd = os.path.join(self.blast_inst_dir, blast_type) + \
                               " -query " + query_path + \
                               " -db " + db + \
                               " -out " + out + \
                               " -word_size " + kwargs.get("word_size", str(3)) + \
                               " -num_threads " + str(num_threads) + \
                               " -outfmt " + str(outfmt) + \
                               " -max_hsps " + str(max_hsps)
            send_log_message("run '%s' command" % blast_search_cmd, mes_type="info", logger=self.logger)
            if self.logger is not None:
                self.logger.info()
            else:
                print("INFO: run '%s' command" % blast_search_cmd)
            subprocess.call(blast_search_cmd, shell=True)

        if kwargs.get("remove_tmp", True):
            shutil.rmtree(self.tmp_dir, ignore_errors=True)
        if read_output:
            with open(out) as out_f:
                return NCBIXML.parse(out_f)
        else:
            return out
