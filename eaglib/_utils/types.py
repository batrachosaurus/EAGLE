import os
import shutil

import numpy as np


def np_memmap_astype(dat_path, old_dtype, new_dtype, shape):
    old_dat_path = dat_path+".old"
    shutil.move(dat_path, old_dat_path)
    memmap_array = np.memmap(dat_path, dtype=new_dtype, mode="w+", shape=shape)
    for i, x in enumerate(np.memmap(old_dat_path, dtype=old_dtype, mode='r', shape=shape)):
        memmap_array[i] = x
    os.remove(old_dat_path)
    return memmap_array


def bool_from_str(string):
    if str(string).lower() in ("0", "false", "f", "no", "n", "na", "none", "null"):
        return False
    else:
        return True