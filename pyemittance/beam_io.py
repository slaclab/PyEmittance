import numpy as np
from os import path
import json
import time

from epics import PV
from pyemittance.otrs_io import get_beamsizes_otrs
from pyemittance.wire_io import get_beamsizes_wire

class MachineIO():
    """Class for handling all machine I/O"""
    def __init__(self, name='LCLS', meas_type='OTRS'):
        # machine name: only LCLS or FACET are currently supported
        self.name = name
        # specify OTRS or WIRE scans
        self.meas_type = meas_type
        self.online = False
        self.use_profmon = False
        self.settle_time = 3 # sleep time in seconds

        # load info about PVs used in measurements (e.g. quad scan PV, image PV)
        self.dir, self.filename = path.split(__file__)
        self.CONFIG_PATH = path.join(self.dir, "configs")
        self.meas_pv_info = json.load(open(self.CONFIG_PATH + '/meas_pv_info.json'))

        self.meas_read_pv = PV(self.meas_pv_info['meas_device']['pv']['read'])
        if self.online:
            # load info about settings to optimize
            self.opt_pv_info = json.load(open(self.CONFIG_PATH + '/opt_pv_info.json'))
            self.opt_pvs = self.opt_pv_info['opt_vars']
            self.meas_cntrl_pv = PV(self.meas_pv_info['meas_device']['pv']['cntrl'])
            self.sol_cntrl_pv = PV(self.opt_pvs[0])
            self.cq_cntrl_pv = PV(self.opt_pvs[1])
            self.sq_cntrl_pv = PV(self.opt_pvs[2])

    def get_beamsizes_machine(self, config, quad_val):
        """Fn that pyemittance.observer calls
        Takes quad value as input,
        Returns xrms, yrms, xrms_err, yrms_err
        """
        if self.online and quad_val is not None:
            self.setquad(quad_val)
            self.setinjector(config)
            time.sleep(self.settle_time)
            
        else:
            print("Not setting online values.")

        if self.meas_type == 'OTRS':
            return get_beamsizes_otrs(self.use_profmon)

        elif self.meas_type == 'WIRE':
            if self.online:
                return get_beamsizes_wire(self.online)
            else:
                print("Not running wire scans.")

        else:
            raise NotImplementedError('No valid measurement type defined.')

    def setinjector(self, set_list):
        if self.online:
            self.sol_cntrl_pv.put(set_list[0])
            self.cq_cntrl_pv.put(set_list[1])
            self.sq_cntrl_pv.put(set_list[2])
        elif not set_list:
            pass
        else:
            print("Not setting inj online values.")
            pass

    def setquad(self, value):
        """Sets Q525 to new scan value"""
        if self.online:
            self.meas_cntrl_pv.put(value)
        else:
            print("Not setting quad online values.")
            pass

    def get_beamsize_inj(self, set_list, quad):
        """Get beamsize fn that changes upstream cu injector
        Returns xrms and yrms in [m]
        """            
        beamsize = self.get_beamsizes_machine(set_list, quad)
        return np.array([beamsize[0], beamsize[1]])
