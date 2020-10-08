import asyncio
import gc
import logging
import multiprocessing as mp
import queue
import random
import string
import sys
import time
from threading import Thread


class AsyncWorker(object):

    def __init__(self):
        self.tasks = dict()
        self.loop = asyncio.new_event_loop()
        self._t = Thread(target=start_background_loop, args=(self.loop,), daemon=True)
        self._t.start()
        self._is_closed = False

    def __getitem__(self, item):
        return self.tasks[item]

    def submit(self, task):
        if self._is_closed:
            raise RuntimeError("A task cannot be submitted into the closed worker")
        if not self.loop.is_running():
            self._t = Thread(target=start_background_loop, args=(self.loop,), daemon=True)
            self._t.start()

        task_id = ''.join([random.choice(string.ascii_letters + string.digits) for n in range(32)])
        self.tasks[task_id] = asyncio.run_coroutine_threadsafe(task.pop("function")(**task), loop=self.loop)
        return task_id

    def join(self):
        self.loop.call_soon_threadsafe(self.loop.stop)
        self._t.join()
        self.loop.close()
        del self  # TODO: is the result available?

    def close(self):
        self._is_closed = True

    def get_result(self, task_id):
        while True:
            task = self.tasks[task_id]
            if task._state == "FINISHED":
                break
            time.sleep(5)
        return task.result()


class ProcessPool(object):

    def __init__(self, workers=1):
        self.workers = workers
        self.queue = mp.Queue()
        self._processes = list()
        self._process_states = mp.Manager().list()
        for i in range(self.workers):
            self._process_states.append(0)
            p = mp.Process(target=self._queue_reader, args=(self.queue, i, self._process_states))
            p.start()
            self._processes.append(p)
        self._is_closed = False

    def submit(self, task):
        if self._is_closed:
            raise RuntimeError("A task cannot be submitted into the closed pool of workers")
        self.queue.put(task)

    @staticmethod
    def _queue_reader(q, p_index, p_states, timeout=1):
        while True:
            try:
                task = q.get_nowait()
                p_states[p_index] = 1
                process_worker(task)
            except queue.Empty:
                time.sleep(timeout)
            p_states[p_index] = 0

    def join(self, timeout=1):
        is_running = True
        while is_running:
            is_running = False
            for p_state in self._process_states:
                if p_state > 0:
                    is_running = True
                    time.sleep(timeout)
        for p in self._processes:
            p.terminate()
        self.queue.close()
        del self

    def close(self):
        self._is_closed = True


def start_background_loop(loop: asyncio.AbstractEventLoop) -> None:
    asyncio.set_event_loop(loop)
    loop.run_forever()
    loop.run_until_complete(loop.shutdown_asyncgens())


def process_worker(kwargs, use_try=False):
    func = kwargs.pop('function', None)
    res = None
    if 'try_err_message' in kwargs.keys():
        use_try = True
    logger_name = kwargs.get('logger_name', None)
    if isinstance(logger_name, str):
        logger = logging.getLogger(logger_name)
    else:
        logger = None
    if callable(func):
        if use_try:
            try:
                res = func(**kwargs)
            except:
                if isinstance(logger, logging.Logger):
                    logger.warning("%s %s" % (kwargs['try_err_message'], sys.exc_info()))
                else:
                    print("%s %s" % (kwargs['try_err_message'], sys.exc_info()))
        else:
            res = func(**kwargs)
    else:
        if logger:
            logger.warning("No function to run")
        else:
            print("No function to run")
    gc.collect()
    return res


worker = process_worker
