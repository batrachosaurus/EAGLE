import time
import queue
import random
import string
import asyncio
from threading import Thread
import multiprocessing as mp


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
            if task.done():
                break
            time.sleep(5)
        return task.result()


class ProcessPool(object):

    def __init__(self, workers=1, autostart=True):
        self.workers = workers
        self.queue = mp.Queue()
        self._processes = list()
        self._process_states = mp.Manager().list()
        if autostart:
            self.start()
        self._is_closed = False

    def start(self):
        for i in range(self.workers):
            self._process_states.append(0)
            p = mp.Process(target=self._queue_reader, args=(self.queue, i, self._process_states))
            p.start()
            self._processes.append(p)

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
        time.sleep(timeout)
        is_running = True
        while is_running:
            is_running = False
            for i, p in enumerate(self._processes):
                if not p.is_alive():
                    self._process_states[i] = 0
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


class TaskResultsStorage(object):

    def __init__(self, storage_time=3600):
        self.tasks = mp.Manager().dict()
        self.named_tasks = mp.Manager().dict()
        self.storage_time = storage_time
        self._p = mp.Process(target=self._del_obsolete_tasks, args=(self.tasks, self.storage_time), daemon=True)
        self._p.start()

    def __getitem__(self, item):
        return self.get_task_result(item)

    @staticmethod
    def _del_obsolete_tasks(tasks, storage_time):
        time.sleep(storage_time)
        for task_id in list(tasks.keys()):
            if time.time()-tasks[task_id]["time"] > storage_time:
                del tasks[task_id]

    def add_task_result(self, task_result, task_name=None):
        task_id = ''.join([random.choice(string.ascii_letters + string.digits) for n in range(32)])
        self.tasks[task_id] = {"result": task_result, "time": time.time()}
        if task_name is not None:
            self.named_tasks = task_id
        return task_id

    def get_task_result(self, task_id_or_name):
        if task_id_or_name in self.tasks:
            return self.tasks[task_id]["result"]
        elif task_id_or_name in self.tasks:
            return self.tasks[self.named_tasks[task_name]]["result"]


def start_background_loop(loop: asyncio.AbstractEventLoop) -> None:
    asyncio.set_event_loop(loop)
    loop.run_forever()
    loop.run_until_complete(loop.shutdown_asyncgens())


def process_worker(kwargs):
    return kwargs.pop("function")(**kwargs)


def process_coroutine_worker(kwargs):
    loop = asyncio.new_event_loop()
    loop.run_until_complete(kwargs.pop("function")(**kwargs))
    loop.close()


worker = process_worker
