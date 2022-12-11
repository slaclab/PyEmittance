import logging
import sys

from . import _version

__version__ = _version.get_versions()["version"]

from pyemittance.pyemittance import PyEmittance

def print_logging(level=logging.INFO):
    """
    Redirects log messages to stdout (simple).
    """
    logging.basicConfig(format='%(asctime)s | %(levelname)s : %(message)s',
                     level=logging.INFO, stream=sys.stdout)
