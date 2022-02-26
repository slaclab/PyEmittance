import numpy as np
import bisect
from pyemittance.beam_io import get_beamsizes_machine

class Observer:
    '''
    Observer reads beamsizes and sets measurement quad
    Observer stores values for beamsizes and quad settings
    '''

    def __init__(self, quad_meas, beam_meas, beam_meas_err):
        self.quad_meas = quad_meas
        self.beam_meas = beam_meas
        self.beam_meas_err = beam_meas_err
        self.use_prev_meas = True 
        self.tolerance = 0.1
        
        # if using the surrogate model
        self.use_model = True
        self.get_beamsizes_model = None
        self.config = None

    def measure_beam(self, quad_list):
        '''ADD ERRORS TOO!!!'''
        xrms = []
        yrms = []
        xrms_err = []
        yrms_err = []
                
        if not self.quad_meas:
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
                    self.beam_meas['x'] = xrms
                    self.beam_meas['y'] = yrms
                    self.beam_meas_err['x'] = xrms_err
                    self.beam_meas_err['y']= yrms_err
                
        else:
            for val in quad_list:              
                # find loc within sorted list
                loc = bisect.bisect_left(self.quad_meas, val)

                if (loc != 0 and 
                    loc != len(self.quad_meas)-1 and 
                    loc != len(self.quad_meas)
                   ):
                    
                    # compare to values before and after
                    diff_prev = abs(val - self.quad_meas[loc-1])
                    diff_next = abs(val - self.quad_meas[loc+1]) 

                if (loc == 0 or 
                    loc == len(self.quad_meas)-1 or 
                    loc == len(self.quad_meas) or
                    (diff_prev > self.tolerance and diff_next > self.tolerance)
                   ):
                    
                    # add in list and measure value
                    self.quad_meas[loc:loc] = [val]

                    # measure bs at this value
                    # returns xrms, yrms, xrms_err, yrms_err
                    beamsizes = self.get_beamsizes(val)

                    # add new quad value in same location
                    self.beam_meas['x'][loc:loc] = [beamsizes[0]]
                    self.beam_meas['y'][loc:loc] = [beamsizes[1]]
                    self.beam_meas_err['x'][loc:loc] = [beamsizes[2]]
                    self.beam_meas_err['y'][loc:loc] = [beamsizes[3]]

                    xrms.append(self.beam_meas['x'][loc])
                    yrms.append(self.beam_meas['y'][loc])
                    xrms_err.append(self.beam_meas_err['x'][loc])
                    yrms_err.append(self.beam_meas_err['y'][loc])

                else: # if either is <= tolerance 
                    if diff_prev <= diff_next:
                        use_loc = loc-1
                    else:
                        use_loc = loc+1

                    # return already measured value (closest)
                    xrms.append(self.beam_meas['x'][use_loc])
                    yrms.append(self.beam_meas['y'][use_loc])
                    xrms_err.append(self.beam_meas_err['x'][use_loc])
                    yrms_err.append(self.beam_meas_err['y'][use_loc])
                
        return xrms, yrms, xrms_err, yrms_err
                                    
        
    def get_beamsizes(self, val):
        '''Define where the beamsizes are acquired from'''
        if self.use_model == True:
            return self.get_beamsizes_model(self.config, val)

        if self.use_model == False:
            return get_beamsizes_machine
        
