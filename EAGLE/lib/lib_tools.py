import os
import pickle
import sys


def worker(kwargs, use_try=False):
    func = kwargs.pop('function', None)
    if 'try_err_message' in kwargs.keys():
        use_try = True
    if func:
        if use_try:
            try:
                func(**kwargs)
            except:
                print kwargs['try_err_message'], sys.exc_info()
        else:
            func(**kwargs)
    else:
        print "No function to run"


def load_list_from_file(f_path, remove_list_f=False):
    loaded_list = []
    f = open(f_path, 'rb')
    for line_ in f:
        line = None
        line = line_.strip()
        loaded_list.append(pickle.loads(line))
    if remove_list_f:
        os.remove(f_path)
    return loaded_list
