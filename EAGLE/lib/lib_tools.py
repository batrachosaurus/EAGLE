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
