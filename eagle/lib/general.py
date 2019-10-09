# This code can have only standard Python imports
import configparser
import gzip
import io
import gc
from collections import OrderedDict
import logging
import logging.config
import multiprocessing as mp
import os
import pickle
import shutil
import numbers
import subprocess
import sys
import time
import string
import random

import yaml
import numpy as np


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
                elif not isinstance(self.__dict__[param], numbers.Number) and \
                        new_value.lower() in ("0", "0.0", "false", "f", "no", "n", "na", "none", "null"):
                    self.__dict__[param] = None
                else:
                    self.__dict__[param] = type(self.__dict__[param])(new_value)

    @staticmethod
    def _get_from_config(config_obj, section, option, fallback):
        try:
            return config_obj.get(section, option)
        except (configparser.NoSectionError, configparser.NoOptionError):
            return fallback


def _config_parser(config_path):
    """ Function parses config file and puts the result into an object of ConfigParser class
      :param config_path: path to config file
      :return: a ConfigParser object
      """
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def setup_logging(default_path,
                  default_level=logging.INFO,
                  env_key='LOG_CFG'):
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):

        with open(path, 'rt') as f:
            string = f.read()
            config = yaml.load(string)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)


class RedisQueue:

    def __init__(self, queue_name, redis_connection, max_size=0):
        self.name = queue_name
        self.conn = redis_connection
        self.max_size = max_size

    def qsize(self):
        return self.conn.llen(self.name)

    def empty(self):
        return self.qsize() == 0

    def full(self):
        return self.qsize() != 0

    def put(self, message):
        if self.max_size > 0:
            while self.qsize() >= self.max_size:
                time.sleep(1)
        self.conn.rpush(self.name, pickle.dumps(message))

    def get(self, block=True, timeout=None):
        if block:
            message = self.conn.blpop(self.name, timeout=timeout)
        else:
            message = self.conn.lpop(self.name)

        if message:
            return pickle.loads(message[1])


def bool_from_str(string):
    if str(string).lower() in ("0", "false", "f", "no", "n", "na", "none", "null"):
        return False
    else:
        return True


def run_proc_pool(num_threads, queue, constant_params=None, end_message="done"):
    #  queue is RedisQueue object
    proc_dict = dict()
    proc_states = mp.Manager().dict()
    for i in range(num_threads):
        proc_states[i] = 0
        p = mp.Process(target=_queue_reader,
                       args=(queue, i, proc_states, constant_params, end_message,))
        p.start()
        proc_dict[i] = p
    time.sleep(10)
    proc_dict = _maintain_procs(proc_dict, proc_states,
                                p_args=[queue, 0, proc_states, constant_params, end_message,],
                                end_message="done")
    for proc_id in proc_dict.keys():
        proc_dict[proc_id].join()
    proc_dict = None


def _queue_reader(queue, proc_num, proc_states, constant_params=None, end_message="done", timeout=1):
    while True:
        proc_states[proc_num] = 0
        if queue.full():
            q_mess = queue.get()
            if q_mess == end_message:
                queue.put(q_mess)
                break
            else:
                if constant_params:
                    q_mess.update(constant_params)
                worker(q_mess)
        else:
            time.sleep(timeout)


def _maintain_procs(p_dict, p_states, p_args, end_message="done", no_response_max=20):
    while len(p_states.keys()) > 0:
        for proc_id in p_states.keys():
            if p_states[proc_id] == end_message:
                p_states.pop(proc_id)
            elif p_states[proc_id] > no_response_max:
                p_args[1] = proc_id
                p = mp.Process(target=_queue_reader,
                               args=tuple(p_args))
                p.start()
                p_dict[proc_id] = p
                p_states[proc_id] = 0
            else:
                p_states[proc_id] += 1
        time.sleep(20)
    return p_dict


def worker(kwargs, use_try=False):
    func = kwargs.pop('function', None)
    if 'try_err_message' in kwargs.keys():
        use_try = True
    logger_name = kwargs.get('logger_name', None)
    if logger_name is not None:
        logger = logging.getLogger(logger_name)
    else:
        logger = None
    if callable(func):
        if use_try:
            try:
                func(**kwargs)
            except:
                if logger:
                    logger.warning("%s %s" % (kwargs['try_err_message'], sys.exc_info()))
                else:
                    print("%s %s" % (kwargs['try_err_message'], sys.exc_info()))
        else:
            func(**kwargs)
    else:
        if logger:
            logger.warning("No function to run")
        else:
            print("No function to run")
    gc.collect()


def filter_list(in_list):
    filtered_list = list()
    for elm_ in in_list:
        elm = None
        elm = elm_.strip()
        if elm:
            filtered_list.append(elm)
    return filtered_list


def reverse_dict(in_dict):
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
    return out_file_path


def get_redis_server(host='localhost', port=6379, restart=True):
    if host == 'localhost' or host == '127.0.0.1':
        if restart:
            subprocess.call("echo 1 > /proc/sys/net/ipv4/tcp_tw_reuse", shell=True)
            subprocess.call("echo 1 > /proc/sys/net/ipv4/tcp_tw_recycle", shell=True)
            subprocess.call("redis-cli -p " + str(port) + " shutdown", shell=True)
        subprocess.Popen("redis-server --port " + str(port), shell=True)
        time.sleep(10)
        return "connected to Redis server"


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


def generate_random_string(l=10):
    return "".join(random.choice(string.ascii_letters + string.digits) for i in range(l))


def np_memmap_astype(dat_path, old_dtype, new_dtype, shape):
    old_dat_path = dat_path+".old"
    shutil.move(dat_path, old_dat_path)
    memmap_array = np.memmap(dat_path, dtype=new_dtype, mode="w+", shape=shape)
    for i, x in enumerate(np.memmap(old_dat_path, dtype=old_dtype, mode='r', shape=shape)):
        memmap_array[i] = x
    os.remove(old_dat_path)
    return memmap_array


def send_log_message(message, mes_type="info", logger=None):
    if isinstance(logger, logging.Logger):
        if mes_type.lower() in ("info", "i"):
            logger.info(message)
        if mes_type.lower() in ("warning", "warn", "w"):
            logger.warning(message)
        if mes_type.lower() in ("error", "err", "e"):
            logger.error(message)
    else:
        if mes_type.lower() in ("info", "i"):
            print("INFO: %s" % message)
        if mes_type.lower() in ("warning", "warn", "w"):
            print("WARNING: %s" % message)
        if mes_type.lower() in ("error", "err", "e"):
            print("ERROR: %s" % message)
