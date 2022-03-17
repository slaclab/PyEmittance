import numpy as np
from pyemittance.optics import estimate_sigma_mat_thick_quad, twiss_and_bmag, get_kL, normalize_emit
from pyemittance.machine_settings import get_twiss0

class EmitCalc:
    """
    Uses info recorded in Observer to do an emittance fit
    """

    def __init__(self, quad_vals=None, beam_vals=None, beam_vals_err=None):
        self.quad_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if quad_vals is None else quad_vals # in kG
        self.beam_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals is None else beam_vals
        self.beam_vals_err = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals_err is None else beam_vals_err

        self.sig_mat_screen = {'x': [], 'y': []}
        self.twiss0 = get_twiss0() # emit, beta, alpha
        self.twiss_screen = {'x': [], 'y': []} # emit, beta, alpha
        self.beta_err = None
        self.alpha_err = None

        self.calc_bmag = False
        self.plot = False
        self.verbose = False

        # Main output of emittance calc
        self.out_dict = {'nemitx': None,
                         'nemity': None,
                         'nemitx_err': None,
                         'nemity_err': None,
                         'bmagx': None,
                         'bmagy': None,
                         'bmagx_err': None,
                         'bmagy_err': None,
                         'opt_q_x': None,
                         'opt_q_y': None}

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

        res = estimate_sigma_mat_thick_quad(bs, kL, bs_err, weights,
                                            calc_bmag=self.calc_bmag, plot=self.plot, verbose=self.verbose)
        if np.isnan(res[0]):
            self.out_dict[f'nemit{dim}'] = np.nan
            self.out_dict[f'nemit{dim}_err'] = np.nan
            self.out_dict[f'bmag{dim}'] = np.nan
            self.out_dict[f'bmag{dim}_err'] = np.nan
            return self.out_dict
        else:
            emit, emit_err, beta_rel_err, alpha_rel_err = res[0:4]
            if self.calc_bmag:
                sig_11, sig_12, sig_22 = res[4:]

        norm_emit_res = normalize_emit(emit, emit_err)
        self.out_dict[f'nemit{dim}'] = normalize_emit(emit, emit_err)[0]
        self.out_dict[f'nemit{dim}_err'] = normalize_emit(emit, emit_err)[1]

        if self.calc_bmag:
            self.sig_mat_screen[dim] = [sig_11, sig_12, sig_22]
            self.beta_err = beta_rel_err
            self.alpha_err = alpha_rel_err

            bmag_calc_res = self.get_twiss_bmag(dim=dim)
            # Get bmag and bmag_err
            self.out_dict[f'bmag{dim}'] = bmag_calc_res[0]
            self.out_dict[f'bmag{dim}_err'] = bmag_calc_res[1]
            # Get best value for scanning quad
            self.out_dict[f'opt_q_{dim}'] = q[bmag_calc_res[2]]

        return self.out_dict

    def get_twiss_bmag(self, dim='x'):

        sig_11 = self.sig_mat_screen[dim][0]
        sig_12 = self.sig_mat_screen[dim][1]
        sig_22 = self.sig_mat_screen[dim][2]

        # twiss0 in x or y AT THE SCREEN
        beta0, alpha0 = self.twiss0[dim][1], self.twiss0[dim][2]

        # return dict of emit, beta, alpha, bmag
        twiss = twiss_and_bmag(sig_11, sig_12, sig_22,
                               self.beta_err, self.alpha_err,
                               beta0=beta0, alpha0=alpha0)
        # Save twiss at screen
        self.twiss_screen[dim] = twiss['emit'], twiss['beta'], twiss['alpha']

        return twiss['bmag'], twiss['bmag_err'], twiss['min_idx']

    def get_gmean_emit(self):

        try:
            nemit = np.sqrt( self.out_dict['nemitx'] * self.out_dict['nemity'] )
            nemit_err = nemit * ( (self.out_dict['nemitx_err']/self.out_dict['nemitx'])**2 +
                                  (self.out_dict['nemity_err']/self.out_dict['nemity'])**2 )**0.5

            self.out_dict['nemit'] = nemit
            self.out_dict['nemit_err'] = nemit_err

        except TypeError:
            self.out_dict['nemit'] = None
            self.out_dict['nemit_err'] = None



