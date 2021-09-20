import gzip
import io
import os
import platform
import shutil
import subprocess

import wget


def join_files(in_files_list, out_file_path, files_transform=None, **kwargs):
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


def gunzip(in_path, out_path, remove_input=True):
    with gzip.open(in_path, 'rb') as input_f_gz, \
            io.open(out_path, 'wb') as output_f:
        shutil.copyfileobj(input_f_gz, output_f)
        output_f.close()
    if remove_input:
        os.remove(in_path)


def compare_files(f1_path, f2_path):
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


def download_file_(file_link, download_dir="./", logger=None):  # This can be reduced with 'urllib.request.urlretrieve()'

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
