import numpy as np
import matplotlib.pyplot as plt
from pyemittance.optics import estimate_sigma_mat_thick_quad, twiss_and_bmag, get_kL, normalize_emit
from pyemittance.beam_io import get_twiss0

class EmitCalc:
    """
    Uses info recorded in Observer to do an emittance fit
    """

    def __init__(self, quad_vals=None, beam_vals=None, beam_vals_err=None):
        self.quad_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if quad_vals is None else quad_vals # in kG
        self.beam_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals is None else beam_vals
        self.beam_vals_err = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals_err is None else beam_vals_err
        self.x_use = np.arange(0, len(self.beam_vals['x']), 1)
        self.y_use = np.arange(0, len(self.beam_vals['y']), 1)

        self.sig_mat_screen = {'x': [], 'y': []}
        self.twiss0 = get_twiss0() # emit, beta, alpha
        self.twiss_screen = {'x': [], 'y': []} # emit, beta, alpha
        self.beta_err = None
        self.alpha_err = None

        self.test_mode = False
        self.noise_red = 50000
        self.plot = True
        self.calc_bmag = False

    def check_conditions(self, ):

        self.x_use = np.arange(0, len(beam_vals['x']), 1)
        self.y_use = np.arange(0, len(beam_vals['y']), 1)

        minx = np.min(self.beam_vals['x'])
        miny = np.min(self.beam_vals['y'])

        self.x_use = np.argwhere(self.beam_vals['x'] < 2.0 * minx)
        self.y_use = np.argwhere(self.beam_vals['y'] < 2.0 * miny)

    def weighting_func(self, beamsizes, beamsizes_err):
        """
        Weigh the fit with Var(sizes) and the sizes themselves
        :param beamsizes: RMS beamsizes measured on screen
        :param err_beamsizes: error on RMS estimate
        :return: weights for fitting
        """
        beamsizes = np.array(beamsizes)
        beamsizes_err = np.array(beamsizes_err)

        sig_bs = 2 * beamsizes * beamsizes_err
        # Here the weight is 1/sigma
        weights = 1 / beamsizes + 1 / sig_bs
        return weights

    def error_propagation(self, gradient):
        """
        Propagate error from var(y) to emittance
        :param gradient: gradient of emittance
        :return: error on emittance from fit
        """
        return np.sqrt( (gradient.T @ self.covariance_matrix) @ gradient)

    def get_emit(self, dim='x'):
        """
        Get emittance at quad from beamsizes and quad scan
        :param dim: 'x' or 'y'
        :return: normalized emittance and error
        """

        q = self.quad_vals[dim]
        # quad vals are passed in machine units
        kL = get_kL(q)

        bs = self.beam_vals[dim]
        bs_err = self.beam_vals_err[dim]

        weights = self.weighting_func(bs, bs_err) # 1/sigma

        if self.test_mode == False:
            res = estimate_sigma_mat_thick_quad(bs, kL, bs_err, weights, 
                                                calc_bmag=self.calc_bmag, plot=self.plot)
            if np.isnan(res[0]):
                return np.nan, np.nan
            else:
                emit, emit_err, beta_rel_err, alpha_rel_err = res[0:4]
                if self.calc_bmag:
                    sig_11, sig_12, sig_22 = res[4:]

        if self.test_mode == True:
            bs = bs + np.random.rand(len(bs)) / self.noise_red
            print("NOISE")
            res = estimate_sigma_mat_thick_quad(bs, kL, bs_err, weights, 
                                                calc_bmag=self.calc_bmag, plot=self.plot)
            if np.isnan(res[0]):
                return np.nan, np.nan
            else:
                emit, emit_err, beta_rel_err, alpha_rel_err = res[0:4]
                if self.calc_bmag:
                    sig_11, sig_12, sig_22 = res[4:]

        emit, emit_err = normalize_emit(emit, emit_err)

        if self.calc_bmag:
            self.sig_mat_screen[dim] = [sig_11, sig_12, sig_22]
            self.beta_err = beta_rel_err
            self.alpha_err = alpha_rel_err

            bmag, bmag_err = self.get_twiss_bmag(dim=dim)
            return emit, emit_err, bmag, bmag_err

        return emit, emit_err

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



