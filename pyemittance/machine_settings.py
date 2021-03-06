# This file contains functions to retrieve  settings for the machine
# (r-matrices, Twiss params, etc)
# These should be specified here or in json config files
import os
import numpy as np
import json

this_dir, this_filename = os.path.split(__file__)
CONFIG_PATH = os.path.join(this_dir, "configs")

def which_machine(filepath=CONFIG_PATH+'/beamline_info.json'):
    """Print which machine settings are being used"""
    beamline_info = json.load(open(filepath))
    name = beamline_info['name']
    print(f"Using {name} beamline info.")

def get_twiss0(filepath= CONFIG_PATH+'/beamline_info.json'):
    """Import Twiss0 from config file"""

    beamline_info = json.load(open(filepath))
    twiss0 = beamline_info['Twiss0']
    # emit, beta, alpha
    twiss0_by_dim = {'x': [twiss0[0], twiss0[2], twiss0[4]],
                     'y': [twiss0[1], twiss0[3], twiss0[5]]
                     }

    return twiss0_by_dim

def get_rmat(filepath=CONFIG_PATH+'/beamline_info.json'):
    """Import r-matrix from config file"""

    beamline_info = json.load(open(filepath))
    # if only separated by a drift:
    # rMat will be [1, L, 0, 1]
    rmat = beamline_info['rMat']
    # m11, m12, m21, m22
    rmat_array = np.array([[rmat[0], rmat[1]],
                            [rmat[2], rmat[3]]
                           ])

    return rmat_array

def get_energy(filepath=CONFIG_PATH+'/beamline_info.json'):
    """Import beam energy from config file [GeV]"""

    beamline_info = json.load(open(filepath))
    energy = beamline_info['energy']
    return energy