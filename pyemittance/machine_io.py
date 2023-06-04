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

    def __init__(self, config_name="LCLS_OTR2", config_dict=None, meas_type="OTRS"):
        # specify OTRS or WIRE scans
        self.meas_type = meas_type
        self.online = False


        # Set configs for measurement
        # if config is not provided, use LCLS OTR2 as default
        if config_dict is None and config_name is None:
            logger.info("No configuration specified. Taking default LCLS-OTR2 configs.")
            self.config_name = "LCLS_OTR2"
            self.config_dict = self.load_config()
        else:
            self.config_name = config_name
            self.config_dict = config_dict if config_dict else self.load_config()

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

    def load_config(self):
        # if path exists, load from path
        if self.config_name is not None:
            self.config_dict = load_configs(self.config_name)
        return self.config_dict

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
            return np.random.uniform(0.5e-4, 5e-4), np.random.uniform(1e-4, 6e-4), 0, 0
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
