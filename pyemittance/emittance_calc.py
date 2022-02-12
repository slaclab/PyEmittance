import numpy as np
import matplotlib.pyplot as plt
from pyemittance.optics import estimate_sigma_mat_thick_quad, twiss_and_bmag
from pyemittance.beam_io import get_twiss0

class EmitCalc:
    '''
    Uses info recorded in Observer to do an emittance fit
    '''

    def __init__(self, quad_vals=None, beam_vals=None):
        self.quad_vals = np.empty(0, ) if quad_vals is None else quad_vals
        self.beam_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals is None else beam_vals
        self.beam_vals_err = {dim: self.beam_vals[dim]*0.05 for dim in self.beam_vals}
        self.x_use = np.arange(0, len(self.beam_vals['x']), 1)
        self.y_use = np.arange(0, len(self.beam_vals['y']), 1)

        self.sig_mat_screen = {'x': [], 'y': []}
        self.twiss0 = get_twiss0() # emit, beta, alpha
        self.twiss_screen = {'x': [], 'y': []} # emit, beta, alpha
        self.beta_err = None
        self.alpha_err = None

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
        sig_bs = 2 * beamsizes * beamsizes_err
        weights = 1 / beamsizes + 1 / sig_bs
        return weights

    def error_propagation(self, gradient):
        '''
        Propagate error from var(y) to emittance
        :param gradient: gradient of emittance
        :return: error on emittance from fit
        '''
        return np.sqrt( (gradient.T @ self.covariance_matrix) @ gradient)

    def get_emit(self, dim='x'):
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
            emit, emit_err, beta_rel_err, alpha_rel_err, sig_11, sig_12, sig_22 = estimate_sigma_mat_thick_quad(bs, q, weights)
            plt.show()

        if self.test_mode == True:
            bs = bs + np.random.rand(len(bs)) / self.noise_red
            print("NOISE")
            emit, emit_err, beta_rel_err, alpha_rel_err, sig_11, sig_12, sig_22 = estimate_sigma_mat_thick_quad(bs, q, weights)
            plt.show()

        err = np.std(np.absolute(np.sqrt(sig_11) - bs))

        self.sig_mat_screen[dim] = [sig_11, sig_12, sig_22]
        self.beta_err = beta_rel_err
        self.alpha_err = alpha_rel_err

        print(f"emit: {emit/1e-6:.3}, emit err: {emit_err/1e-6:.3}, bs er: {err/1e-6:.3}")
        return emit, err

    def get_twiss_bmag(self, dim='x'):

        sig_11, sig_12, sig_22 = self.sig_mat_screen[dim][0], self.sig_mat_screen[dim][1], self.sig_mat_screen[dim][2]

        # twiss0 in x or y AT THE SCREEN
        beta0, alpha0 = self.twiss0[dim][1], self.twiss0[dim][2]

        # return dict of emit, beta, alpha, bmag
        twiss = twiss_and_bmag(sig_11, sig_12, sig_22,
                               self.beta_err, self.alpha_err,
                               beta0=beta0, alpha0=alpha0)
        # Save twiss at screen
        self.twiss_screen[dim] = twiss['emit'], twiss['beta'], twiss['alpha']

        return twiss['bmag'], twiss['bmag_err']



