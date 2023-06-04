from pyemittance.tools import NpEncoder, isotime
import numpy as np
import json
from epics import caget, caget_many

import logging
logger = logging.getLogger(__name__)

def save_image(im, nrow, ncol,  timestamp, impath="", avg_img=True):
    """Saves images with col,row info and corresp. settings"""

    if avg_img:

        np.save(str(impath) + f"img_avg_{timestamp}.npy", im)
        np.save(str(impath) + f"ncol_avg_{timestamp}.npy", ncol)
        np.save(str(impath) + f"nrow_avg_{timestamp}.npy", nrow)

    else:

        np.save(str(impath) + f"img_{timestamp}.npy", im)
        np.save(str(impath) + f"ncol_{timestamp}.npy", ncol)
        np.save(str(impath) + f"nrow_{timestamp}.npy", nrow)


def save_config(
    xrms,
    yrms,
    xrms_err,
    yrms_err,
    timestamp,
    meas_read_pv,
    im=None,
    configpath="",
    impath="",
):
    if timestamp is None:
        f = open(configpath + "bax_beamsize_config_info.csv", "a+")
        timestamp = isotime()
    else:
        f = open(configpath + "beamsize_config_info.csv", "a+")

    # todo make more general, pandas etc

    bact_cur = caget(meas_read_pv)
    f.write(
        f"{timestamp},{varx_cur},{vary_cur},{varz_cur},"
        f"{bact_cur},{xrms},{yrms},{xrms_err},{yrms_err}\n"
    )

    f.close()

    if im:
        np.save((str(impath) + f"img_config_{timestamp}.npy", im.proc_image))


def save_emit_run(out_dict, path=""):
    timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
    with open(path + f"pyemittance_data_{timestamp}.json", "w") as outfile:
        json.dump(out_dict, outfile, cls=NpEncoder)
