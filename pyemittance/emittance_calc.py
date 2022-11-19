import numpy as np
from os import path, makedirs
import errno
import os
from pyemittance.emittance_calc_base import EmitCalcBase
from pyemittance.optics import estimate_sigma_mat_thick_quad, get_kL, normalize_emit
from pyemittance.machine_settings import get_quad_len
from pyemittance.saving_io import save_emit_run


class EmitCalc(EmitCalcBase):
    """ Quad scan emittance measurement type
    Uses info recorded in Observer to do an emittance fit
    """
    def init_class_attr(self):
        self.quad_vals = self.meas_vals
        self.quad_len = get_quad_len(self.config_dict['beamline_info'])
        # Initialize paths and dirs for saving
        self.init_saving()

    def define_default_config(self):
        print("No configuration specified. Taking default LCLS-OTR2 configs.")
        self.config_name = "LCLS_OTR2"
        self.config_dict = self.load_config()

    def get_emit(self):
        """
        Get emittance at quad from beamsizes and quad scan
        :return: normalized emittance and error
        """

        for dim in self.dims:
            # run emit calc for x and y

            q = self.quad_vals[dim]
            # quad vals are passed in machine units
            kL = get_kL(q, self.quad_len, self.energy)

            bs = self.beam_vals[dim]
            bs_err = self.beam_vals_err[dim]

            weights = self.weighting_func(bs, bs_err) # 1/sigma

            # Storing quadvals and beamsizes in self.out_dict for plotting purposes
            self.out_dict[f'quadvals{dim}'] = list(q)
            self.out_dict[f'beamsizes{dim}'] = list(bs)
            self.out_dict[f'beamsizeserr{dim}'] =  list(bs_err)

            res = estimate_sigma_mat_thick_quad(bs, kL, bs_err, weights, dim=dim, Lquad=self.quad_len,
                                                energy=self.energy, rmat=self.rmat, calc_bmag=self.calc_bmag,
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

            norm_emit_res = normalize_emit(emit, emit_err, self.energy)
            self.out_dict[f'nemit{dim}'] = norm_emit_res[0]
            self.out_dict[f'nemit{dim}_err'] = norm_emit_res[1]

            if self.calc_bmag:
                self.sig_mat_meas[dim] = [sig_11, sig_12, sig_22]
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
            mkdir_p(savepaths['images'])
            mkdir_p(savepaths['summaries'])
            mkdir_p(savepaths['fits'])
            mkdir_p(savepaths['raw_saves'])
        except OSError:
            print("Savepaths not set. Please set them in 'configs/savepaths.json'")
            from pathlib import Path
            parent = Path(__file__).resolve().parent
            examples_dir = str(parent)[:-11] + "docs/examples"
            print("Using docs examples directory: ", examples_dir)
            savepaths['images'] = examples_dir + "/saved_images/"
            savepaths['summaries'] = examples_dir + "/summaries/"
            savepaths['fits'] = examples_dir + "/saved_fits/"
            savepaths['raw_saves'] = examples_dir + "/raw_saves/"
            mkdir_p(savepaths['images'])
            mkdir_p(savepaths['summaries'])
            mkdir_p(savepaths['fits'])
            mkdir_p(savepaths['raw_saves'])

        # Start headings
        file_exists = path.exists(savepaths['summaries'] + "image_acq_quad_info.csv")

        if not file_exists:
            f = open(savepaths['summaries'] + "image_acq_quad_info.csv", "a+")
            f.write(
                f"{'timestamp'},{'ncol'},{'nrow'},{'roi_xmin'},{'roi_xmax'}"
                f",{'roi_ymin'},{'roi_ymax'},{'resolution'},{'bact'},"
                f"{'x_size'},{'y_size'},{'xrms'},{'yrms'},"
                f"{'xrms_err'},{'yrms_err]'}\n")
            f.close()

        file_exists = path.exists(savepaths['summaries'] + "beamsize_config_info.csv")

        if not file_exists:
            f = open(savepaths['summaries'] + "beamsize_config_info.csv", "a+")
            f.write(
                f"{'timestamp'},{'varx_cur'},{'vary_cur'},{'varz_cur'},"
                f"{'bact_cur'},{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'}\n")
            f.close()
