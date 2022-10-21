import numpy as np
from os import path, makedirs
import errno
import os
from pyemittance.emittance_calc_base import EmitCalcBase
from pyemittance.optics import estimate_sigma_mat_multiwire, normalize_emit
from pyemittance.machine_settings import get_rmat_wires, get_loc_wires
from pyemittance.saving_io import save_emit_run


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
            # res is [emit, emit_err, beta_err/beta, alpha_err/alpha, s11, s12, s22]
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

            norm_emit_res = normalize_emit(emit, emit_err, self.energy)
            self.out_dict[f'nemit{dim}'] = norm_emit_res[0]
            self.out_dict[f'nemit{dim}_err'] = norm_emit_res[1]

            if self.calc_bmag:
                # bmag is calculated at the first wire

                self.sig_mat_meas[dim] = [sig_11, sig_12, sig_22]
                self.beta_err = beta_rel_err
                self.alpha_err = alpha_rel_err

                bmag_calc_res = self.get_twiss_bmag(dim=dim)
                # Get bmag and bmag_err
                self.out_dict[f'bmag{dim}'] = bmag_calc_res[0]
                self.out_dict[f'bmag{dim}_err'] = bmag_calc_res[1]

        # get geometric mean (set to None when x and y are not calculated)
        self.get_gmean_emit()

        if self.save_runs:
            self.save_run()

        return self.out_dict

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
            examples_dir = str(parent)[:-11] + "examples"
            print("Using examples directory: ", examples_dir)
            savepaths['summaries'] = examples_dir + "/summaries/"
            savepaths['fits'] = examples_dir + "/saved_fits/"
            savepaths['raw_saves'] = examples_dir + "/raw_saves/"
            mkdir_p(savepaths['summaries'])
            mkdir_p(savepaths['fits'])
            mkdir_p(savepaths['raw_saves'])

        file_exists = path.exists(savepaths['summaries'] + "beamsize_config_info.csv")

        if not file_exists:
            # TODO: check how many wires are saving and save appr. num of meas.
            f = open(savepaths['summaries'] + "beamsize_config_info.csv", "a+")
            f.write(
                f"{'timestamp'},{'varx_cur'},{'vary_cur'},{'varz_cur'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'},"
                f"{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'}\n"
            )
            f.close()


