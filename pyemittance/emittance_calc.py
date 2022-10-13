import numpy as np
from pyemittance.optics import estimate_sigma_mat_thick_quad, twiss_and_bmag, get_kL, normalize_emit
from pyemittance.machine_settings import get_twiss0
from pyemittance.saving_io import save_emit_run
from pyemittance.load_json_configs import load_configs


class EmitCalc:
    """
    Uses info recorded in Observer to do an emittance fit
    """
    def __init__(self,
                 quad_vals: list = None,
                 beam_vals: list = None,
                 beam_vals_err: list = None,
                 config: dict = None,
                 config_path: str = None
                 ):

        self.quad_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if quad_vals is None else quad_vals # in kG
        self.beam_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals is None else beam_vals

        # Make sure error is added to beamsizes if none is provided
        if beam_vals_err is None or sum(beam_vals_err['x'])==0 or sum(beam_vals_err['y'])==0:
            self.bs_error = (0.015, 0.015) # Define some error on beamsizes in each dimension
            self.beam_vals_err = {'x': np.asarray(self.beam_vals['x'])*self.bs_error[0],
                                  'y': np.asarray(self.beam_vals['y'])*self.bs_error[1]}
        else:
            self.beam_vals_err = beam_vals_err

        # if config is not provided, use LCLS OTR2 as default
        if config is None and config_path is None:
            print("No configuration specified. Taking default LCLS-OTR2 configs.")
            self.config_path = "/LCLS_OTR2"
            self.config = self.load_config()
        else:
            self.config_path = config_path
            self.config = config if config else self.load_config()

        self.dims = ['x', 'y'] # TODO: make code use provided in self.dims instead of hardcoding 'x' and 'y'
        self.sig_mat_screen = {'x': [], 'y': []}
        self.twiss0 = get_twiss0() # emit, beta, alpha
        self.twiss_screen = {'x': [], 'y': []} # emit, beta, alpha
        self.beta_err = None
        self.alpha_err = None

        self.calc_bmag = False
        self.plot = False
        self.verbose = False
        self.save_runs = False

        # Main output of emittance calc
        self.out_dict = {'nemitx': None,
                         'nemity': None,
                         'nemitx_err': None,
                         'nemity_err': None,
                         'nemit': None,
                         'nemit_err': None,
                         'bmagx': None,
                         'bmagy': None,
                         'bmagx_err': None,
                         'bmagy_err': None,
                         'bmag_emit': None,
                         'bmag_emit_err': None,
                         'opt_q_x': None,
                         'opt_q_y': None,
                         'total_points_measured': None
                         }

    def load_config(self):
        # if path exists, load from path
        if self.config_path is not None:
            self.config = load_configs(self.config_path)
        return self.config

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

    def get_emit(self):
        """
        Get emittance at quad from beamsizes and quad scan
        :param dim: 'x' or 'y'
        :return: normalized emittance and error
        """

        for dim in self.dims:
            # run emit calc for x and y

            q = self.quad_vals[dim]
            # quad vals are passed in machine units
            kL = get_kL(q)

            bs = self.beam_vals[dim]
            bs_err = self.beam_vals_err[dim]

            weights = self.weighting_func(bs, bs_err) # 1/sigma

            # Storing quadvals and beamsizes in self.out_dict for plotting purposes
            self.out_dict[f'quadvals{dim}'] = list(q)
            self.out_dict[f'beamsizes{dim}'] = list(bs)
            self.out_dict[f'beamsizeserr{dim}'] =  list(bs_err)

            res = estimate_sigma_mat_thick_quad(bs, kL, bs_err, weights, dim=dim,
                                                calc_bmag=self.calc_bmag,
                                                plot=self.plot, verbose=self.verbose)
            if np.isnan(res[0]):
                self.out_dict['nemitx'], self.out_dict['nemity'] = np.nan, np.nan
                self.out_dict['nemitx_err'], self.out_dict['nemity_err'] = np.nan, np.nan
                self.out_dict['bmagx'], self.out_dict['bmagy'] = np.nan, np.nan
                self.out_dict['bmagx_err'], self.out_dict['bmagy_err'] = np.nan, np.nan
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

        # get geometric mean (set to None when x and y are not calculated)
        self.get_gmean_emit()

        if self.save_runs:
            self.save_run()

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

            if self.out_dict['bmagx'] is not None and self.out_dict['bmagy'] is not None:
                nbmag = np.sqrt( self.out_dict['bmagx'] * self.out_dict['bmagy'] )
                bmag_emit_err = nemit*nbmag * (
                    (self.out_dict['nemitx_err']/self.out_dict['nemitx'])**2 +
                    (self.out_dict['nemity_err']/self.out_dict['nemity'])**2 +
                    (self.out_dict['bmagx_err']/self.out_dict['bmagx'])**2 +
                    (self.out_dict['bmagy_err']/self.out_dict['bmagy'])**2)**0.5
                self.out_dict['bmag_emit'] = nemit * nbmag
                self.out_dict['bmag_emit_err'] = bmag_emit_err

        except TypeError:
            self.out_dict['nemit'] = np.nan
            self.out_dict['nemit_err'] = np.nan
            self.out_dict['bmag_emit'] = np.nan
            self.out_dict['bmag_emit_err'] = np.nan

    def save_run(self):
        save_emit_run(self.out_dict)