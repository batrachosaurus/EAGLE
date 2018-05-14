# This code can have only standard Python imports
import ConfigParser
import multiprocessing as mp
import sys
import time
import shutil
from collections import OrderedDict

import pandas


class ConfBase(object):

    def __init__(self, config_path):
        self.config_path = None
        self.config = None
        if config_path:
            self.update_by_config(config_path=config_path)

    def update_by_config(self, config_path):
        self.config_path = config_path
        self.config = _config_parser(config_path=config_path)
        config_sections = self.config.sections()
        for param in self.__dict__.keys():
            new_value = None
            for config_section in config_sections:
                if self._get_from_config(self.config, config_section, param, fallback=None):
                    new_value = self.config.get(config_section, param)
                    break
            if new_value:
                if type(self.__dict__[param]) is list:
                    self.__dict__[param] = [elm.strip() for elm in new_value.split(",")]
                elif type(self.__dict__[param]) is bool:
                    self.__dict__[param] = bool_from_str(new_value)
                elif new_value.lower() in ("0", "false", "f", "no", "n", "na", "none", "null"):
                    self.__dict__[param] = None
                else:
                    self.__dict__[param] = type(self.__dict__[param])(new_value)

    @staticmethod
    def _get_from_config(config_obj, section, option, fallback):
        try:
            return config_obj.get(section, option)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            return fallback


def _config_parser(config_path):
    """ Function parses config file and puts the result into an object of ConfigParser class
      :param config_path: path to config file
      :return: a ConfigParser object
      """
    config = ConfigParser.ConfigParser()
    config.read(config_path)
    return config


def bool_from_str(string):
    if str(string).lower() in ("0", "false", "f", "no", "n", "na", "none", "null"):
        return False
    else:
        return True


def run_proc_pool(num_threads, queue, constant_params=None, end_message="done"):
    #  queue is multiprocessing Manager list
    proc_list = list()
    for i in range(num_threads):
        p = mp.Process(target=_queue_reader,
                       args=(queue, constant_params, end_message,))
        p.start()
        proc_list.append(p)
    time.sleep(10)
    for proc in proc_list:
        proc.join()
    queue = None
    time.sleep(10)
    proc_list = None


def _queue_reader(queue, constant_params=None, end_message="done", timeout=1):
    while True:
        if queue:
            q_mess = queue.pop(0)
            if q_mess == end_message:
                queue.append(q_mess)
                break
            else:
                if constant_params:
                    q_mess.update(constant_params)
                worker(q_mess)
        else:
            time.sleep(timeout)


def worker(kwargs, use_try=False):
    func = kwargs.pop('function', None)
    if 'try_err_message' in kwargs.keys():
        use_try = True
    logger = kwargs.get('logger', None)
    if func:
        if use_try:
            try:
                func(**kwargs)
            except:
                if logger:
                    logger.warning("%s %s" % (kwargs['try_err_message'], sys.exc_info()))
                else:
                    print "%s %s" % (kwargs['try_err_message'], sys.exc_info())
        else:
            func(**kwargs)
    else:
        if logger:
            logger.warning("No function to run")
        else:
            print "No function to run"


def filter_list(in_list):
    filtered_list = list()
    for elm_ in in_list:
        elm = None
        elm = elm_.strip()
        if elm:
            filtered_list.append(elm)
    return filtered_list


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
        fasta_f.write(seq_id+"\n")
        fasta_f.write(fasta_dict[seq_id]+"\n")
    fasta_f.close()


def load_phylip_dist_matrix(matrix_path):
    matr_f = open(matrix_path)
    lines_dict = OrderedDict()
    seqs_list = list()
    for line_ in matr_f:
        line = None
        line = line_.strip()
        if not line:
            continue
        line_list = filter_list(line.split())
        if len(line_list) == 1:
            continue
        seqs_list.append(line_list[0])
        lines_dict[line_list[0]] = line_list[1:]
    dist_matrix = pandas.DataFrame.from_dict(data=lines_dict, orient='index')
    dist_matrix.columns = seqs_list
    return dist_matrix


def dump_phylip_dist_matrix(dist_matrix, matrix_path):
    matr_f = open(matrix_path, 'w')
    matr_f.write("    %s\n" % len(dist_matrix.columns))
    for seq in dist_matrix.index:
        num_spaces_to_add = 10 - len(seq)
        spaces_to_add = [" " for i in range(num_spaces_to_add)]
        matr_f.write("%s %s\n" % (seq+"".join(spaces_to_add), " ".join(dist_matrix.loc[seq].tolist())))
    matr_f.close()


def reduce_seq_names(fasta_dict, num_letters=10, num_words=4):
    if num_letters < 6:
        print("Number of letters must be at least 6")
        return 1
    if num_words < 2:
        print("Number of words must be at least 2")
        return 1
    splitters_repl = {"\t": " ",
                      ",": " ",
                      ";": " ",
                      ".": " ",
                      ":": " ",
                      "|": " ",
                      "/": " ",
                      "\\": " "}
    parts_size_list = _get_part_size_list(num_letters, num_words)
    reduced_fasta_dict = dict()
    seq_names_dict = dict()
    for seq_name in fasta_dict.keys():
        if len(seq_name) <= num_letters:
            prepared_seq_name = None
            prepared_seq_name = seq_name+"".join("_" for i in range(num_letters-len(seq_name)))
            seq_names_dict[prepared_seq_name] = seq_name
            reduced_fasta_dict[prepared_seq_name] = fasta_dict[seq_name]
            continue
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
        un_num = 1
        un_fix = get_un_fix(un_num, res_len)
        while seq_names_dict.get(reduced_seq_name+un_fix, None):
            un_fix = None
            un_num += 1
            un_fix = get_un_fix(un_num, res_len)
        reduced_fasta_dict[reduced_seq_name+un_fix] = fasta_dict[seq_name]
        seq_names_dict[reduced_seq_name+un_fix] = seq_name
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


def get_un_fix(un_num, fix_len):
    un_codes = ["_", '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E']
    # 'N' - undefined (num duplicates is bigger than len(un_codes))
    if fix_len == 1:
        try:
            return un_codes[un_num-1]
        except IndexError:
            return 'N'
    elif fix_len == 0:
        return ""
    elif un_num <= len(un_codes):
        return un_codes[0] + get_un_fix(un_num, fix_len - 1)
    else:
        filled_rank = len(un_codes)**(fix_len-1)
        return un_codes[un_num//filled_rank - 1] + get_un_fix(un_num % filled_rank, fix_len - 1)


def load_newick(newick_f_path):
    newick_f = open(newick_f_path)
    tree_list = list()
    for line_ in newick_f:
        line = None
        line = line_.strip()
        tree_list.append(line)
        if line[-1] == ";":
            break
    return "".join(tree_list)


def dump_tree_newick(tree_newick, newick_f_path):
    newick_f = open(newick_f_path, "w")
    newick_f.write(tree_newick)
    newick_f.close()


def join_files(in_files_list, out_file_path):
    with open(out_file_path, 'wb') as out_file:
        for f_path in in_files_list:
            f = open(f_path, 'rb')
            shutil.copyfileobj(f, out_file)
            f.close()
        out_file.close()
