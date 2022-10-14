import numpy as np
import time

from epics import PV
from pyemittance.saving_io import save_config
from pyemittance.load_json_configs import load_configs


class MachineIO():
    """Class for handling all machine I/O"""
    def __init__(self, config_name='LCLS_OTR2', config_dict=None, meas_type='OTRS'):
        # specify OTRS or WIRE scans
        self.meas_type = meas_type
        self.online = False
        self.use_profmon = False
        self.settle_time = 3  # sleep time in seconds

        # Set configs for measurement
        # if config is not provided, use LCLS OTR2 as default
        if config_dict is None and config_name is None:
            print("No configuration specified. Taking default LCLS-OTR2 configs.")
            self.config_name = "LCLS_OTR2"
            self.config_dict = self.load_config()
        else:
            self.config_name = config_name
            self.config_dict = config_dict if config_dict else self.load_config()

        self.meas_pv_info = self.config_dict['meas_pv_info']
        self.meas_read_pv = PV(self.meas_pv_info['meas_device']['pv']['read'])

        # load info about settings to optimize
        self.opt_pv_info = self.config_dict['opt_pv_info']
        self.opt_pvs = self.opt_pv_info['opt_vars']
        self.meas_cntrl_pv = PV(self.meas_pv_info['meas_device']['pv']['cntrl'])
        self.sol_cntrl_pv = PV(self.opt_pvs[0])
        self.cq_cntrl_pv = PV(self.opt_pvs[1])
        self.sq_cntrl_pv = PV(self.opt_pvs[2])

    def load_config(self):
        # if path exists, load from path
        if self.config_name is not None:
            self.config_dict = load_configs(self.config_name)
        return self.config_dict

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
            print("Running offline.")

        if self.meas_type == 'OTRS' and self.online:
            from pyemittance.otrs_io import get_beamsizes_otrs
            return get_beamsizes_otrs(self.config_dict, self.use_profmon)
        elif self.meas_type == 'WIRE' and self.online:
            from pyemittance.wire_io import get_beamsizes_wire
            print("Running wire scanner")
            return get_beamsizes_wire(self.online, self.config_dict)
        elif not self.online:
            return np.random.uniform(0.5e-4,5e-4), np.random.uniform(1e-4,6e-4), 0, 0
        else:
            raise NotImplementedError('No valid measurement type defined.')

    def setinjector(self, set_list):
        if self.online and set_list is not None:
            self.sol_cntrl_pv.put(set_list[0])
            self.cq_cntrl_pv.put(set_list[1])
            self.sq_cntrl_pv.put(set_list[2])
        elif not set_list:
            pass
        else:
            print("Not setting injector online values.")
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
        # save BAX beam size data
        save_config(beamsize[0], beamsize[1], beamsize[2], beamsize[3], None)
        return np.array([beamsize[0], beamsize[1]])
