import bisect
import numpy as np

class Observer:
    """
    Observer reads beamsizes and sets measurement quad
    Observer stores values for beamsizes and quad settings
    """

    def __init__(self, quad_meas, beam_meas, beam_meas_err):
        self.quad_meas = quad_meas
        self.beam_meas = beam_meas
        self.beam_meas_err = beam_meas_err
        self.use_prev_meas = False
        self.tolerance = 0.1

        self.inj_config = None  # injector configuration settings for SOL, CQ, SQ
        self.use_model = True
        # if using the surrogate model
        self.get_beamsizes_model = None
        self.add_noise = False
        self.noise_red = 100000
        
        # if using machine
        self.use_profmon = False
        self.online = False
        self.config_name = 'sim'
        self.config_dict = None
        self.meas_type = 'OTRS'

        self.verbose = True

    def measure_beam(self, quad_list):
        xrms = []
        yrms = []
        xrms_err = []
        yrms_err = []

        if not self.quad_meas or self.use_prev_meas is False:
            # if no measurements exist yet, measure all
            for val in quad_list:
                # measure bs at this value
                beamsizes = self.get_beamsizes(val)
                xrms.append(beamsizes[0])
                yrms.append(beamsizes[1])
                xrms_err.append(beamsizes[2])
                yrms_err.append(beamsizes[3])

                # update saved values
                self.quad_meas.append(val)
                self.beam_meas['x'].append(xrms[-1])
                self.beam_meas['y'].append(yrms[-1])
                self.beam_meas_err['x'].append(xrms_err[-1])
                self.beam_meas_err['y'].append(yrms_err[-1])
                
        else:
            for val in quad_list:
                # find loc within sorted list
                loc = bisect.bisect_left(self.quad_meas, val)

                if loc != 0 and loc != len(self.quad_meas)-1 and loc < len(self.quad_meas):
                    # compare to values before and after
                    diff_prev = abs(val - self.quad_meas[loc-1])
                    diff_next = abs(self.quad_meas[loc] - val)
                elif loc == 0:
                    diff_prev = np.inf
                    diff_next = abs(self.quad_meas[loc] - val)
                elif loc == len(self.quad_meas)-1:
                    diff_prev = abs(val - self.quad_meas[loc-1])
                    diff_next = abs(self.quad_meas[loc] - val)
                elif loc >= len(self.quad_meas):
                    diff_prev = abs(val - self.quad_meas[-1])
                    diff_next = np.inf

                if diff_prev > self.tolerance and diff_next > self.tolerance or loc >= len(self.quad_meas):
                    # if no neighboring value is within tol
                    # or if value has not been measured
                    
                    # add in list and measure value
                    self.quad_meas.insert(loc, val)

                    # measure bs at this value
                    # returns xrms, yrms, xrms_err, yrms_err
                    beamsizes = self.get_beamsizes(val)

                    # add new quad value in same location
                    self.beam_meas['x'].insert(loc, beamsizes[0])
                    self.beam_meas['y'].insert(loc, beamsizes[1])
                    self.beam_meas_err['x'].insert(loc, beamsizes[2])
                    self.beam_meas_err['y'].insert(loc, beamsizes[3])

                    xrms.append(self.beam_meas['x'][loc])
                    yrms.append(self.beam_meas['y'][loc])
                    xrms_err.append(self.beam_meas_err['x'][loc])
                    yrms_err.append(self.beam_meas_err['y'][loc])

                else: # if either is <= tolerance
                    if diff_prev <= diff_next:
                        use_loc = loc-1
                    else:
                        use_loc = loc

                    # return already measured value (closest)
                    xrms.append(self.beam_meas['x'][use_loc])
                    yrms.append(self.beam_meas['y'][use_loc])
                    xrms_err.append(self.beam_meas_err['x'][use_loc])
                    yrms_err.append(self.beam_meas_err['y'][use_loc])

        return xrms, yrms, xrms_err, yrms_err

    def get_beamsizes(self, val):
        """Define where the beamsizes are acquired from"""
        if self.use_model:
            if self.add_noise:
                # this is outdated/not used right now
                beamsizes = self.get_beamsizes_model(self.inj_config, val)
                xrms = beamsizes[0] + np.random.rand(1) / self.noise_red
                yrms = beamsizes[1] + np.random.rand(1) / self.noise_red
                xrms_err = beamsizes[2] + np.random.rand(1) / self.noise_red / 10
                yrms_err = beamsizes[3] + np.random.rand(1) / self.noise_red / 10
                return xrms[0], yrms[0], xrms_err[0], yrms_err[0]

            if len(np.asarray(self.inj_config).shape) == 2:
                # temp workaround for testing
                out = self.get_beamsizes_model(self.inj_config[0]+[val])
                return np.array([out["sigma_x"][0], out["sigma_y"][0], 0., 0.])

            return self.get_beamsizes_model(self.inj_config, val)

        else:
            from pyemittance.beam_io import MachineIO
            io = MachineIO(self.config_name, self.config_dict, self.meas_type)
            io.online = self.online
            io.use_profmon = self.use_profmon
            # note we are not setting the injector on the machine here
            return io.get_beamsizes_machine(self.inj_config, val)
        
