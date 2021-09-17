# This code can have only standard Python imports
import gzip
import io
import re
from collections import OrderedDict
import os
import shutil
import subprocess
import string
import random
import platform

import wget
import numpy as np


def filter_list(in_list):  # This can be reduced with 'list(filter(lambda li: li.strip(), in_list))'
    filtered_list = list()
    for elm_ in in_list:
        elm = None
        elm = elm_.strip()
        if elm:
            filtered_list.append(elm)
    return filtered_list


def revert_dict(in_dict):  # This can be reduced with '{v: k for k, v in in_dict}'
    out_dict = OrderedDict()
    for key in in_dict:
        out_dict[in_dict[key]] = key
    return out_dict


def get_un_fix(un_num, fix_len):
    un_codes = ["_", '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E']
    # 'N' - undefined (num duplicates is bigger than len(un_codes))
    if fix_len == 1:
        try:
            return un_codes[un_num]
        except IndexError:
            return 'N'
    elif fix_len == 0:
        return ""
    elif un_num < len(un_codes):
        return un_codes[0] + get_un_fix(un_num, fix_len - 1)
    else:
        filled_rank = len(un_codes)**(fix_len-1)
        return un_codes[un_num//filled_rank] + get_un_fix(un_num % filled_rank, fix_len - 1)


def join_files(in_files_list, out_file_path, files_transform=None, **kwargs):  # TODO: move to eaglib._utils.files
    if type(in_files_list) not in (list, tuple):
        if in_files_list == out_file_path:
            return 1
        in_files_LIST = [in_files_list]
        in_files_list = None
        in_files_list = in_files_LIST
        in_files_LIST = None
    with open(out_file_path, 'wb') as out_file:
        for f_path in in_files_list:
            f = open(f_path, 'rb')
            if callable(files_transform):
                shutil.copyfileobj(files_transform(f, **kwargs), out_file)
            else:
                shutil.copyfileobj(f, out_file)
            f.close()
        out_file.close()
    if kwargs.get("remove_infiles", False):
        for f_path in in_files_list:
            os.remove(f_path)
    return out_file_path


def gunzip(in_path, out_path, remove_input=True):  # TODO: move to eaglib._utils.files
    with gzip.open(in_path, 'rb') as input_f_gz, \
            io.open(out_path, 'wb') as output_f:
        shutil.copyfileobj(input_f_gz, output_f)
        output_f.close()
    if remove_input:
        os.remove(in_path)


def compare_files(f1_path, f2_path):  # TODO: move to eaglib._utils.files
    # returns True if files are equal else returns False
    f1 = open(f1_path, 'rb')
    f2 = open(f2_path, 'rb')
    f1_lines = f1.readlines()
    f2_lines = f2.readlines()
    f1.close()
    f2.close()
    if len(f1_lines) != len(f2_lines):
        return False
    for i in range(len(f1_lines)):
        if f1_lines[i].strip() != f2_lines[i].strip():
            return False
    return True


def generate_random_string(l=10):  # TODO: move to eaglib._utils.strings
    return "".join(random.choice(string.ascii_letters + string.digits) for i in range(l))


def np_memmap_astype(dat_path, old_dtype, new_dtype, shape):  # TODO: move to eaglib._utils.types
    old_dat_path = dat_path+".old"
    shutil.move(dat_path, old_dat_path)
    memmap_array = np.memmap(dat_path, dtype=new_dtype, mode="w+", shape=shape)
    for i, x in enumerate(np.memmap(old_dat_path, dtype=old_dtype, mode='r', shape=shape)):
        memmap_array[i] = x
    os.remove(old_dat_path)
    return memmap_array


def fullmatch_regexp_list(pattern, target_list):  # TODO: move to eaglib._utils.strings
    return list(map(lambda x: re.fullmatch(pattern, x), target_list))


def download_file(file_link, download_dir="./", logger=None): # TODO: find more pythonic way; move to eaglib._utils.files

    if platform.system() == 'Windows':
        try:
            wget.download(file_link, out=download_dir)
        except IOError:
            if logger is not None:
                logger.warning("'%s' file has not been found" % file_link)
            else:
                print("'%s' file has not been found" % file_link)
    else:
        subprocess.call("wget " + file_link + " -P " + download_dir + "/", shell=True)
    return os.path.join(download_dir, os.path.basename(file_link))


def bool_from_str(string):  # TODO: move to eaglib._utils.types
    if str(string).lower() in ("0", "false", "f", "no", "n", "na", "none", "null"):
        return False
    else:
        return True
