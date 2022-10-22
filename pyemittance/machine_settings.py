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
    """
    Import Twiss0 from config file
    """

    beamline_info = beamline_info_config_dict
    twiss0 = beamline_info['Twiss0']
    # emit, beta, alpha
    twiss0_by_dim = {'x': [twiss0[0], twiss0[2], twiss0[4]],
                     'y': [twiss0[1], twiss0[3], twiss0[5]]
                     }

    return twiss0_by_dim

def get_rmat(beamline_info_config_dict):
    """
    Import r-matrix from config file
    This function is used in the quad scan meas.
    """

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


def get_rmat_wires(beamline_info_config_dict):
    """
    Import r-matrix from config file
    This function is used in the multiwire meas.
    """
    beamline_info = beamline_info_config_dict

    # setup for 3 or 4 wires
    # if there is only 3, fourth is ignored in calculations
    rmat_wires = {'wire1': None,
                  'wire2': None,
                  'wire3': None,
                  'wire4': None
                  }
    for w in rmat_wires.keys():
        if w in beamline_info.keys():
            rmat_wires[w]['rMatx'] = beamline_info[w]['rMatx']
            rmat_wires[w]['rMaty'] = beamline_info[w]['rMaty']
        else:
            if w == "wire1":
                raise Exception("Wire configuration in beamline_info.json is set incorrectly.")
            del rmat_wires[w]
        # m11, m12, m21, m22
        # in the form below in the Matlab model:
        # R{n} = Rab(1:4, 1: 4);
        # Rx(n,:)=[Rab(1, 1), Rab(1, 2)];
        # Ry(n,:)=[Rab(3, 3), Rab(3, 4)];
        rmat_wires[w]['rMatx'] = np.array([[rmat_wires[w]['rMatx'][0], rmat_wires[w]['rMatx'][1]],
                                           [rmat_wires[w]['rMatx'][2], rmat_wires[w]['rMatx'][3]]
                                           ])
        rmat_wires[w]['rMaty'] = np.array([[rmat_wires[w]['rMaty'][0], rmat_wires[w]['rMaty'][1]],
                                           [rmat_wires[w]['rMaty'][2], rmat_wires[w]['rMaty'][3]]
                                           ])
    return rmat_wires


def get_loc_wires(beamline_info_config_dict):
    """
    Import wire locations from config file
    This function is used in the multiwire meas.
    """
    beamline_info = beamline_info_config_dict

    # setup for 3 or 4 wires
    # if there is only 3, fourth is ignored in calculations
    # set it to None in json file if possible
    loc_wires = {'wire1': None,
                 'wire2': None,
                 'wire3': None,
                 'wire4': None
                 }
    for w in loc_wires.keys():
        if w in beamline_info.keys():
            loc_wires[w] = beamline_info[w]['location']
        else:
            if w == "wire1":
                raise Exception("Wire configuration in beamline_info.json is set incorrectly.")
            del loc_wires[w]

    return loc_wires

