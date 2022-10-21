"""Test case for the multiprocessing logger."""
import multiprocessing as mp
from ase.clease.mp_logger import MultiprocessHandler
import logging as lg
import os
from ase.test import NotAvailable

if __name__ != "__main__":
    msg = "Pickle error when file is not top level module"
    raise NotAvailable(msg)

log_file = "mp_logger_test.txt"
logger = lg.getLogger(__name__)
logger.setLevel(lg.INFO)
handler = MultiprocessHandler(log_file)
logger.addHandler(handler)


def function_logging(msg):
    """Log message via logger."""
    logger.info(msg)


args = ["Indx: {}".format(i) for i in range(32)]

workers = mp.Pool(2)
workers.map(function_logging, args)

# Try the read the logging file back and ensure that
# it has the expected format
record_found = [False for _ in range(len(args))]
with open(log_file, 'r') as infile:
    lines = infile.readlines()

for line in lines:
    # This will raise a ValueError if it does not find the line
    # in args. In that case something went wrong with the logger
    indx = args.index(line[:-1])  # Skip the trailing '\n'

os.remove(log_file)
