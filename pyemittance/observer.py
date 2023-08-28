import bisect
import numpy as np
from pyemittance.machine_io import MachineIO
from time import sleep

import logging
logger = logging.getLogger(__name__)

class Observer:
    """
    Observer reads beamsizes and sets measurement quad
    Observer stores values for beamsizes and quad settings
    """

    def __init__(self, config_dict=None,
                 quad_meas=None,
                 beam_meas=None,
                 beam_meas_err=None,
                 use_prev_meas = False,
                 tolerance=0.1,
                 online=False,
                 meas_type = None
                ):
        
        self.config_dict = config_dict        
        self.quad_meas = quad_meas
        self.beam_meas = beam_meas
        self.beam_meas_err = beam_meas_err
        self.use_prev_meas = use_prev_meas
        self.tolerance = tolerance
        self.online = online
        self.meas_type = meas_type
        
        # Init PVs
        self.io = MachineIO(self.config_dict, meas_type=self.meas_type, online=self.online)

    def measure_beam(self, quad_list):
        xrms = []
        yrms = []
        xrms_err = []
        yrms_err = []

        if not self.quad_meas or self.use_prev_meas is False:
            # if no measurements exist yet, measure all
            for val in quad_list:
                # measure bs at this value
                bdat = self.get_beamsizes(val)
                xrms.append(    bdat['xrms'])
                yrms.append(    bdat['yrms'])
                xrms_err.append(bdat['xrms_err'])
                yrms_err.append(bdat['yrms_err'])

                # update saved values
                self.quad_meas.append(val)
                self.beam_meas["x"].append(xrms[-1])
                self.beam_meas["y"].append(yrms[-1])
                self.beam_meas_err["x"].append(xrms_err[-1])
                self.beam_meas_err["y"].append(yrms_err[-1])

        else:
            for val in quad_list:
                # find loc within sorted list
                loc = bisect.bisect_left(self.quad_meas, val)

                if (
                    loc != 0
                    and loc != len(self.quad_meas) - 1
                    and loc < len(self.quad_meas)
                ):
                    # compare to values before and after
                    diff_prev = abs(val - self.quad_meas[loc - 1])
                    diff_next = abs(self.quad_meas[loc] - val)
                elif loc == 0:
                    diff_prev = np.inf
                    diff_next = abs(self.quad_meas[loc] - val)
                elif loc == len(self.quad_meas) - 1:
                    diff_prev = abs(val - self.quad_meas[loc - 1])
                    diff_next = abs(self.quad_meas[loc] - val)
                elif loc >= len(self.quad_meas):
                    diff_prev = abs(val - self.quad_meas[-1])
                    diff_next = np.inf

                if (
                    diff_prev > self.tolerance
                    and diff_next > self.tolerance
                    or loc >= len(self.quad_meas)
                ):
                    # if no neighboring value is within tol
                    # or if value has not been measured

                    # add in list and measure value
                    self.quad_meas.insert(loc, val)

                    # measure bs at this value
                    # returns xrms, yrms, xrms_err, yrms_err
                    bdat = self.get_beamsizes(val)

                    # add new quad value in same location
                    self.beam_meas["x"].insert(loc, bdat['xrms'])
                    self.beam_meas["y"].insert(loc, bdat['yrms'])
                    self.beam_meas_err["x"].insert(loc, bdat['xrms_err'])
                    self.beam_meas_err["y"].insert(loc, bdat['yrms_err'])

                    xrms.append(self.beam_meas["x"][loc])
                    yrms.append(self.beam_meas["y"][loc])
                    xrms_err.append(self.beam_meas_err["x"][loc])
                    yrms_err.append(self.beam_meas_err["y"][loc])

                else:  # if either is <= tolerance
                    if diff_prev <= diff_next:
                        use_loc = loc - 1
                    else:
                        use_loc = loc

                    # return already measured value (closest)
                    xrms.append(self.beam_meas["x"][use_loc])
                    yrms.append(self.beam_meas["y"][use_loc])
                    xrms_err.append(self.beam_meas_err["x"][use_loc])
                    yrms_err.append(self.beam_meas_err["y"][use_loc])

        return xrms, yrms, xrms_err, yrms_err

    def get_beamsizes(self, val):
        #self.io = MachineIO(self.config_name, self.config_dict, self.meas_type)
        #sleep(.1) 
        return self.io.get_beamsizes_machine(val)
