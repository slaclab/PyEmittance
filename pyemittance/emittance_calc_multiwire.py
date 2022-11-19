import numpy as np
from os import path, makedirs
import errno
import os
from pyemittance.emittance_calc_base import EmitCalcBase
from pyemittance.optics import estimate_sigma_mat_multiwire, twiss_and_bmag
from pyemittance.machine_settings import get_rmat_wires, get_loc_wires
from pyemittance.saving_io import save_emit_run
from pyemittance.load_json_configs import load_configs


class MultiWireCalc(EmitCalcBase):
    """
    Multiwire emittance measurement type
    """
    def init_class_attr(self):
        # here self.meas_vals attribute is constant for x and y
        self.wire_rmat = get_rmat_wires(self.config_dict['beamline_info'])
        self.wire_loc = get_loc_wires(self.config_dict['beamline_info'])
        self.meas_vals = self.wire_rmat
        # Initialize paths and dirs for saving
        self.init_saving()

    def define_default_config(self):
        print("No configuration specified. Taking default LCLS-WS02 configs.")
        self.config_name = "LCLS_WS02"
        self.config_dict = self.load_config()

    def load_config(self):
        # if path exists, load from path
        if self.config_name is not None:
            self.config_dict = load_configs(self.config_name)
        return self.config_dict

    def get_emit(self):
        """
        Get emittance at the FIRST wire from beamsizes and wire locations
        :return: normalized emittance and error
        """

        for dim in self.dims:
            # run emit calc for x and y

            bs = self.beam_vals[dim]
            bs_err = self.beam_vals_err[dim]

            weights = self.weighting_func(bs, bs_err) # 1/sigma

            # Storing quadvals and beamsizes in self.out_dict for plotting purposes
            self.out_dict[f'locations'] = list(self.wire_loc.values())
            self.out_dict[f'beamsizes{dim}'] = list(bs)
            self.out_dict[f'beamsizeserr{dim}'] = list(bs_err)

            res = estimate_sigma_mat_multiwire(self.wire_rmat, bs, bs_err, weights, self.out_dict[f'locations'],
                                               dim=dim, plot=self.plot, verbose=self.verbose)

            # Add all results
            print(res)
            self.output.update(res)

            # Skip further calcs if there was an error
            if res[f'error_{dim}']:
                continue

            if self.calc_bmag:
                # TODO: implement match at first wire
                # if dim == 'x':
                #     sig_11 = res['screen_sigma_11']
                #     sig_12 = res['screen_sigma_12']
                #     sig_22 = res['screen_sigma_22']
                #
                # else:
                #     sig_11 = res['screen_sigma_33']
                #     sig_12 = res['screen_sigma_34']
                #     sig_22 = res['screen_sigma_44']
                # self.sig_mat_screen[dim] = [sig_11, sig_12, sig_22]
                #
                # beta_rel_err = res[f'beta_{dim}_rel_err']
                # alpha_rel_err = res[f'alpha_{dim}_rel_err']
                #
                # self.beta_err = beta_rel_err
                # self.alpha_err = alpha_rel_err
                #
                # bmag_calc_res = self.get_twiss_bmag(dim=dim)
                # # Get bmag and bmag_err
                # self.output[f'screen_bmag{dim}'] = bmag_calc_res[0]
                # self.output[f'screen_bmag{dim}_err'] = bmag_calc_res[1]
                print("Match not implemented for multiwire scan yet.")
                pass

        # get geometric mean if possible
        if (not self.output['error_x']) and (not self.output['error_y']) :
            self.get_gmean_emit()

        if self.save_runs:
            self.save_run()

        return self.out_dict

    def get_twiss_bmag(self, dim='x'):
        '''Not Implemented'''
        # TODO: implement match at first wire
        # sig_11 = self.sig_mat_screen[dim][0]
        # sig_12 = self.sig_mat_screen[dim][1]
        # sig_22 = self.sig_mat_screen[dim][2]
        #
        # # twiss0 in x or y AT THE SCREEN
        # beta0, alpha0 = self.twiss0[dim][1], self.twiss0[dim][2]
        #
        # # return dict of emit, beta, alpha, bmag
        # twiss = twiss_and_bmag(sig_11, sig_12, sig_22,
        #                        self.beta_err, self.alpha_err,
        #                        beta0=beta0, alpha0=alpha0)
        # # Save twiss at screen
        # self.twiss_screen[dim] = [twiss['emit'], twiss['beta'], twiss['alpha']]
        #
        # return twiss['bmag'], twiss['bmag_err'], twiss['min_idx']
        pass

    def get_gmean_emit(self):

        try:
            nemit = np.sqrt( self.output['norm_emit_x'] * self.output['norm_emit_y'] )
            nemit_err = nemit * ( (self.output['norm_emit_x_err']/self.output['norm_emit_x'])**2 +
                                  (self.output['norm_emit_y_err']/self.output['norm_emit_y'])**2 )**0.5

            self.output['sqrt_norm_emit_4d'] = nemit
            self.output['sqrt_norm_emit_4d_err'] = nemit_err

            # TODO: implement match stats at first wire
            # if 'bmag_x' in self.output and 'bmag_y' in self.output:
            #     nbmag = np.sqrt( self.output['bmag_x'] * self.output['bmag_y'] )
            #     bmag_emit_err = nemit*nbmag * (
            #         (self.output['norm_emit_x_err']/self.output['norm_emit_x'])**2 +
            #         (self.output['norm_emit_y_err']/self.output['norm_emit_y'])**2 +
            #         (self.output['bmag_x_err']/self.output['bmag_x'])**2 +
            #         (self.output['bmag_y_err']/self.output['bmag_y'])**2)**0.5
            #     self.output['bmag_emit'] = nemit * nbmag
            #     self.output['bmag_emit_err'] = bmag_emit_err

        except TypeError:
            self.output['sqrt_norm_emit_4d'] = np.nan
            self.output['sqrt_norm_emit_4d_err'] = np.nan
            self.output['bmag_emit'] = np.nan
            self.output['bmag_emit_err'] = np.nan

    def save_run(self):
        save_emit_run(self.out_dict, path=self.config_dict['savepaths']['fits'])

    def init_saving(self):
        """Initialize dirs and files for saving"""

        savepaths = self.config_dict['savepaths']

        def mkdir_p(path_):
            """Set up dirs for results in working dir"""
            try:
                makedirs(path_)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(path_):
                    pass
                else:
                    raise

        # Make directories if needed
        try:
            mkdir_p(savepaths['summaries'])
            mkdir_p(savepaths['fits'])
            mkdir_p(savepaths['raw_saves'])
        except OSError:
            print("Savepaths not set. Please set them in 'configs/savepaths.json'")
            from pathlib import Path
            parent = Path(__file__).resolve().parent
            examples_dir = str(parent)[:-11] + "docs/examples"
            print("Using docs examples directory: ", examples_dir)
            savepaths['summaries'] = examples_dir + "/summaries/"
            savepaths['fits'] = examples_dir + "/saved_fits/"
            savepaths['raw_saves'] = examples_dir + "/raw_saves/"
            mkdir_p(savepaths['summaries'])
            mkdir_p(savepaths['fits'])
            mkdir_p(savepaths['raw_saves'])

        file_exists = path.exists(savepaths['summaries'] + "beamsize_config_info.csv")

        if not file_exists:
            # TODO: check how many wires are saving and save appropriately.
            # This is hardcoded for single measurements right now.
            f = open(savepaths['summaries'] + "beamsize_config_info.csv", "a+")
            f.write(
                f"{'timestamp'},{'varx_cur'},{'vary_cur'},{'varz_cur'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'}\n"
            )
            f.close()


