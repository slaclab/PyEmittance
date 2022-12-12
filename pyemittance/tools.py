import numpy as np
import json
import datetime
import os
import errno


def isotime():
    return (
        datetime.datetime.utcnow()
        .replace(tzinfo=datetime.timezone.utc)
        .astimezone()
        .replace(microsecond=0)
        .isoformat()
    )


def mkdir_p(path):
    """Set up dirs for results in working dir"""
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


class NpEncoder(json.JSONEncoder):
    """
    See: https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable/50916741
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)
