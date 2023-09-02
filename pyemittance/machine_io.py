import numpy as np
import time

from epics import PV
from pyemittance.load_json_configs import load_configs
from pyemittance.otrs_io import get_beamsizes_otrs
from pyemittance.wire_io import get_beamsizes_wire

import logging
logger = logging.getLogger(__name__)

class MachineIO:
    """Class for handling all machine I/O"""

    def __init__(self, config_dict=None, meas_type="OTRS", online=False):
        # specify OTRS or WIRE scans
        self.meas_type = meas_type
        self.online = online
        
        # Allow str to load
        if isinstance(config_dict, str):
            config_dict = load_configs(self.config_name)        
        self.config_dict = config_dict
        
        self.meas_pv_info = self.config_dict["meas_pv_info"]

        meas_device = self.meas_pv_info["meas_device"]
        if "settle_time" in meas_device:
            self.settle_time = meas_device["settle_time"]
        else:
            self.settle_time = 0
            logger.warning("No settle_time found in in meas_device, setting to zero")

        if "bounds" in meas_device:
            self.bounds = meas_device["bounds"]
        else:
            self.bounds = None

        # Connect to PVs 
        self.meas_read_pv = PV(self.meas_pv_info["meas_device"]["pv"]["read"])
        self.meas_cntrl_pv = PV(self.meas_pv_info["meas_device"]["pv"]["cntrl"])


    def get_beamsizes_machine(self, quad_val):
        """Fn that pyemittance.observer calls
        Takes quad value as input,
        Returns xrms, yrms, xrms_err, yrms_err
        """
        if self.online and quad_val is not None:
            self.setquad(quad_val)
            if self.settle_time > 0:  
                logger.info(f"Settling for {self.settle_time} s...")
                time.sleep(self.settle_time)
        
        if self.meas_type == "OTRS" and self.online:
            return get_beamsizes_otrs(self.config_dict)
        elif self.meas_type == "WIRE" and self.online:
            logger.info("Running wire scanner")
            return get_beamsizes_wire(self.online, self.config_dict)
        elif not self.online:
            # Some dummy
            return {'xrms': np.random.uniform(0.5e-4, 5e-4),
                    'yrms': np.random.uniform(1e-4, 6e-4),
                    'xrms_err': 0,
                    'yrms_err': 0,
                   }
        else:
            raise NotImplementedError("No valid measurement type defined.")
    def setquad(self, value):
        """Sets quad to new scan value"""

        if self.bounds is not None:
            lower, upper = self.bounds
            if value < lower:
                raise ValueError(f'Quad set point value {value} < {lower} lower limit')
            if value > upper:
                raise ValueError(f'Quad set point value {value} > {upper} upper limit')                
        
        if self.online:
            logger.info(f'EPICS put {self.meas_cntrl_pv.pvname} = {value}')
            self.meas_cntrl_pv.put(value)
        else:
            logger.info("Not setting quad online values.")
            pass
