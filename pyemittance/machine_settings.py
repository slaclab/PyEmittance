# This file contains functions to retrieve  settings for the machine
# (r-matrices, Twiss params, etc)
# These should be specified here or in json config files
import numpy as np


def which_machine(beamline_info_config_dict):
    """Print which machine settings are being used"""
    beamline_info = beamline_info_config_dict
    name = beamline_info['name']
    print(f"Using {name} beamline info.")

def get_twiss0(beamline_info_config_dict):
    """Import Twiss0 from config file"""

    beamline_info = beamline_info_config_dict
    twiss0 = beamline_info['Twiss0']
    # emit, beta, alpha
    twiss0_by_dim = {'x': [twiss0[0], twiss0[2], twiss0[4]],
                     'y': [twiss0[1], twiss0[3], twiss0[5]]
                     }

    return twiss0_by_dim

def get_rmat(beamline_info_config_dict):
    """Import r-matrix from config file"""

    beamline_info = beamline_info_config_dict
    # if only separated by a drift:
    # rMat will be [1, L, 0, 1]

    # try:
    rmatx = beamline_info['rMatx']
    rmaty = beamline_info['rMaty']
    # m11, m12, m21, m22
    rmatx_array = np.array([[rmatx[0], rmatx[1]],
                            [rmatx[2], rmatx[3]]
                           ])
    rmaty_array = np.array([[rmaty[0], rmaty[1]],
                        [rmaty[2], rmaty[3]]
                       ])

    return rmatx_array, rmaty_array


def get_energy(beamline_info_config_dict):
    """Import beam energy from config file [GeV]"""

    beamline_info = beamline_info_config_dict
    energy = beamline_info['energy']
    return energy

def get_quad_len(beamline_info_config_dict):
    """Import quad len from config file [m]"""

    beamline_info = beamline_info_config_dict
    l = beamline_info['l']
    return l