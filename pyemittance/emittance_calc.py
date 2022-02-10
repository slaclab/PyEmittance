import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
from pyemittance.optics import quad_drift_mat2, propagate_sigma, estimate_sigma_mat_thick_quad_drift

class EmitCalc:
    '''
    Uses info recorded in Observer to do an emittance fit
    '''

    def __init__(self, quad_vals=None, beam_vals=None):
        self.quad_vals = np.empty(0, ) if quad_vals is None else quad_vals
        self.beam_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals is None else beam_vals
        self.beam_vals_err = {dim: beam_vals[dim]*0.1 for dim in beam_vals}
        self.x_use = np.arange(0, len(beam_vals['x']), 1)
        self.y_use = np.arange(0, len(beam_vals['y']), 1)

        self.sig_11 = []
        self.sig_12 = []
        self.sig_22 = []
        self.sig_quad = []
        self.covariance_matrix = None
        
        self.test_mode = False
        self.noise_red = 50000

    def check_conditions(self, ):

        self.x_use = np.arange(0, len(beam_vals['x']), 1)
        self.y_use = np.arange(0, len(beam_vals['y']), 1)

        minx = np.min(self.beam_vals['x'])
        miny = np.min(self.beam_vals['y'])

        self.x_use = np.argwhere(self.beam_vals['x'] < 2.0 * minx)
        self.y_use = np.argwhere(self.beam_vals['y'] < 2.0 * miny)

    def weighting_func(self, beamsizes, beamsizes_err):
        '''
        Weigh the fit with Var(sizes) and the sizes themselves
        :param beamsizes: RMS beamsizes measured on screen
        :param err_beamsizes: error on RMS estimate
        :return: weights for fitting
        '''
        var_bs = ( 2 * beamsizes * beamsizes_err )**2
        weights = 1 / beamsizes + 1 / var_bs
        return weights

    def error_propagation(self, gradient):
        '''
        Propagate error from var(y) to emittance
        :param gradient: gradient of emittance
        :return: error on emittance from fit
        '''
        return np.sqrt( (gradient.T @ self.covariance_matrix) @ gradient)

    def get_emit_gradient(self, sig_ele, emit):
        '''
        Gradient of emittance to calculate cov of parameters estimated
        :param s11: ME11 of r matrix at screen
        :param s12: ME12 of r matrix at screen
        :param s22: ME12 of r matrix at screen
        @param emit: emittance from fit
        :return: gradient of emittance
        '''
        s11 = sig_ele[0]
        s12 = sig_ele[1]
        s22 = sig_ele[2]
        emit_gradient = 1/(2*emit) * np.array([[s11, -2*s12, s22]]).T
        return emit_gradient

    def do_emit_fit(self, dim='x'):
        '''
        Get emittance at quad from beamsizes and quad scan
        :param dim: 'x' or 'y'
        :return: emittance and error
        '''

        # todo update based on x_use, y_use for throwing away fit points
        q = self.quad_vals

        bs = self.beam_vals[dim]
        bs_err = self.beam_vals_err[dim]

        weights = self.weighting_func(bs, bs_err)

        if self.test_mode == False:
            emit, self.sig_11, self.sig_12, self._sig_22 = self.estimate_sigma_mat_thick_quad_drift(bs, q, weights)
            plt.show()

        if self.test_mode == True:
            bs = bs + np.random.rand(len(bs)) / self.noise_red
            print("NOISE")
            emit, self.sig_11, self.sig_12, self._sig_22 = self.estimate_sigma_mat_thick_quad_drift(bs, q, weights)
            plt.show()

        err = np.std(np.absolute(self.sig_11 - bs))
        
        return emit, err

    self.sig_quad = [s11, s12, s22]


