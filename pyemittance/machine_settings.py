# This file contains functions to retrieve  settings for the machine
# (r-matrices, Twiss params, etc)
# These should be specified here or in json config files
import numpy as np

import logging
logger = logging.getLogger(__name__)


def get_twiss0(beamline_info_config_dict):
    """Import Twiss0 from config file"""

    beamline_info = beamline_info_config_dict
    twiss0 = beamline_info["Twiss0"]
    # emit, beta, alpha
    twiss0_by_dim = {
        "x": [twiss0[0], twiss0[2], twiss0[4]],
        "y": [twiss0[1], twiss0[3], twiss0[5]],
    }

    return twiss0_by_dim


def get_rmat(beamline_info_config_dict):
    """Import r-matrix from config file"""
    rMatx = np.array(beamline_info_config_dict["rMatx"]).reshape(2, 2)
    rMaty = np.array(beamline_info_config_dict["rMaty"]).reshape(2, 2)
    return rMatx, rMaty

def get_energy(beamline_info_config_dict):
    """Import beam energy from config file [GeV]"""

    beamline_info = beamline_info_config_dict
    energy = beamline_info["energy"]
    return energy


def get_quad_len(beamline_info_config_dict):
    """Import quad len from config file [m]"""

    beamline_info = beamline_info_config_dict
    L = beamline_info["Lquad"]
    return L
