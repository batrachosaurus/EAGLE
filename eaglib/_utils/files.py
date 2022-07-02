import gzip
import os
import shutil
import urllib.request

from deprecated import deprecated


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
    with gzip.open(in_path, 'rb') as in_f:
        with open(out_path, 'wb') as out_f:
            shutil.copyfileobj(in_f, out_f)
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


@deprecated(reason="too simple code to be a separate function")
def download_file_(file_link, download_dir="./", logger=None):
    local_f_path = os.path.join(download_dir, file_link.split("/")[-1])
    with open(local_f_path, 'wb') as local_f:
        local_f.write(urllib.request.urlopen(file_link).read())
