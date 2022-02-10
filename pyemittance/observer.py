import numpy as np

class Observer:
    '''
    Observer reads beamsizes and sets measurement quad
    Observer stores values for beamsizes and quad settings
    '''

    def __init__(self, quad_meas, beam_meas):
        self.quad_meas = quad_meas
        self.beam_meas = beam_meas
        self.use_model = True

    def measure_beam(self, ):

        if self.use_model == True:
            xrms = np.array(model1(self.quad_meas[-1]))
            yrms = np.array(model1(self.quad_meas[-1]))

        if self.use_model == False:
            # todo: add code for machine io
            pass

        self.beam_meas['x'] = np.concatenate((self.beam_meas['x'], np.array([xrms])))
        self.beam_meas['y'] = np.concatenate((self.beam_meas['y'], np.array([yrms])))

        return xrms, yrms

    def set_quad(self, quad_val):

        if self.use_model == True:
            pass

        if self.use_model == False:
            # todo: add code for machine io
            # set quad
            # read quad back as quad_val
            pass

        self.quad_meas = np.concatenate((self.quad_meas, np.array([quad_val])))