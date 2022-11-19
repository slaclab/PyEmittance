import numpy as np
from pyemittance.machine_settings import get_twiss0, get_rmat, get_energy
from pyemittance.load_json_configs import load_configs
from pyemittance.optics import twiss_and_bmag


class EmitCalcBase:
    """
    Base class for emittance measurement types
    """
    def __init__(self,
                 meas_vals: dict = None,  # quad vals or wire locations/rMats
                 beam_vals: dict = None,
                 beam_vals_err: dict = None,
                 config_dict: dict = None,
                 config_name: str = None
                 ):

        self.meas_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if meas_vals is None else meas_vals
        self.beam_vals = {'x': np.empty(0, ), 'y': np.empty(0, )} if beam_vals is None else beam_vals

        # Make sure error is added to beamsizes if none is provided
        if beam_vals_err is None or sum(beam_vals_err['x']) == 0 or sum(beam_vals_err['y']) == 0:
            self.bs_error = (0.015, 0.015)  # Define some error on beamsizes in each dimension
            self.beam_vals_err = {'x': np.asarray(self.beam_vals['x'])*self.bs_error[0],
                                  'y': np.asarray(self.beam_vals['y'])*self.bs_error[1]}
        else:
            self.beam_vals_err = beam_vals_err

        # if config is not provided, use LCLS OTR2 as default
        if config_dict is None and config_name is None:
            print("No configuration specified.")
            self.define_default_config()
        else:
            self.config_name = config_name
            self.config_dict = config_dict if config_dict else self.load_config()

        self.dims = ['x', 'y'] # TODO: use user passed dimensions in self.dims, and make it extensible
        self.sig_mat_meas = {'x': [], 'y': []}
        # Twiss0 should be for the first wire of multiwire measurements
        self.twiss0 = get_twiss0(self.config_dict['beamline_info'])  # emit0, beta0, alpha0
        self.twiss_screen = {'x': [], 'y': []} # emit, beta, alpha
        self.beta_err = None
        self.alpha_err = None
        self.energy = get_energy(self.config_dict['beamline_info'])
        self.rmat = get_rmat(self.config_dict['beamline_info'])

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
                         }

        # Define class specific attributes
        self.init_class_attr()

    def init_class_attr(self):
        pass

    def define_default_config(self):
        pass

    def load_config(self):
        # if path exists, load from path
        if self.config_name is not None:
            self.config_dict = load_configs(self.config_name)
        return self.config_dict

    def weighting_func(self, beamsizes, beamsizes_err):
        """
        Weigh the fit with Var(sizes) and the sizes themselves
        :param beamsizes: RMS beamsizes measured on screen
        :param beamsizes_err: error on RMS estimate
        :return: weights for fitting
        """
        beamsizes = np.array(beamsizes)
        beamsizes_err = np.array(beamsizes_err)

        sig_bs = 2 * beamsizes * beamsizes_err
        # Here the weight is 1/sigma
        weights = 1 / beamsizes + 1 / sig_bs
        return weights

    def get_emit(self):
        """
        Get emittance at location specified in configs from beamsizes and meas. values
        :return: normalized emittance and error (dict)
        """
        pass

    def get_twiss_bmag(self, dim='x'):

        sig_11 = self.sig_mat_meas[dim][0]
        sig_12 = self.sig_mat_meas[dim][1]
        sig_22 = self.sig_mat_meas[dim][2]

        # twiss0 in x or y at the screen/wire
        beta0, alpha0 = self.twiss0[dim][1], self.twiss0[dim][2]

        # return dict of emit, beta, alpha, bmag
        twiss = twiss_and_bmag(sig_11, sig_12, sig_22,
                               self.beta_err, self.alpha_err,
                               beta0=beta0, alpha0=alpha0)
        # Save twiss at screen/wire
        self.twiss_screen[dim] = [twiss['emit'], twiss['beta'], twiss['alpha']]

        return twiss['bmag'], twiss['bmag_err'], twiss['min_idx']

    def get_emittance(self):
        """Run emittance calc in optics module
        and return sigma matrix
        """
        pass

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

