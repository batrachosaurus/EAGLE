import sys

from EAGLE.constants import EAGLE_logger


def worker(kwargs, use_try=False):
    func = kwargs.pop('function', None)
    if 'try_err_message' in kwargs.keys():
        use_try = True
    if func:
        if use_try:
            try:
                func(**kwargs)
            except:
                EAGLE_logger.warning("%s %s" % (kwargs['try_err_message'], sys.exc_info()))
        else:
            func(**kwargs)
    else:
        EAGLE_logger.warning("No function to run")


def read_fasta_to_dict(fasta_path):
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
