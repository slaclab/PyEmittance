import numpy as np
import json
import datetime


def isotime():
    return (
        datetime.datetime.utcnow()
        .replace(tzinfo=datetime.timezone.utc)
        .astimezone()
        .replace(microsecond=0)
        .isoformat()
    )


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
