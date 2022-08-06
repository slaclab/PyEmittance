# This file contains functions to retrieve  settings for the machine
# (r-matrices, Twiss params, etc)
# These should be specified here or in json config files
import os
import numpy as np
import json

### '''''CHANGE HERE ''' #TODO: update config settings
meas_type = 'OTRS'
#################

if meas_type == 'WIRE':
    add_path = '/LCLS_WS02'
elif meas_type == 'OTRS':
    add_path = '/LCLS2_OTR3'
else:
    add_path = ''

this_dir, this_filename = os.path.split(__file__)
CONFIG_PATH = os.path.join(this_dir, 'configs' + add_path)

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

    try:
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
    except KeyError:
        # TODO: dumb hack for old json 
        rmat = beamline_info['rMat']
        rmat = [np.array([[rmat[0], rmat[1]],[rmat[2], rmat[3]]]), np.array([[rmat[0], rmat[1]],[rmat[2], rmat[3]]])]
        return rmat

def get_energy(filepath=CONFIG_PATH+'/beamline_info.json'):
    """Import beam energy from config file [GeV]"""

    beamline_info = json.load(open(filepath))
    energy = beamline_info['energy']
    return energy

def get_quad_len(filepath=CONFIG_PATH+'/beamline_info.json'):
    """Import quad len from config file [m]"""

    beamline_info = json.load(open(filepath))
    l = beamline_info['l']
    return l