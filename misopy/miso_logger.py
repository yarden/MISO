import os
import sys
import time

import misopy
import logging
import logging.handlers

loggers = {}

def get_logger(logger_name, log_outdir,
               level=logging.WARNING,
               include_stdout=True):
    """
    Return a logging object.
    """
    global loggers
    # Avoid race-conditions
    try:
        os.makedirs(log_outdir)
    except OSError:
        pass
    if loggers.get(logger_name):
        return loggers.get(logger_name)
    logger = logging.getLogger(logger_name)
    formatter = \
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                          datefmt='%m/%d/%Y %I:%M:%S %p')
    # Do not log to file
    #if log_outdir is not None:
    #    log_filename = os.path.join(log_outdir, "%s.log" %(logger_name))
    #    fh = logging.FileHandler(log_filename)
    #    fh.setLevel(level)
    #    fh.setFormatter(formatter)
    #    logger.addHandler(fh)
    logging.root.setLevel(level)
    # Optionally add handler that streams all logs
    # to stdout
    if include_stdout:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    logger.info("Created logger %s" %(logger_name))
    loggers.update({logger_name: logger})
    return logger
