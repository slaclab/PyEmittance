import numpy as np
import os
from pathlib import Path

from pyemittance.optics import (
    estimate_sigma_mat_thick_quad,
    twiss_and_bmag,
    kL_from_machine_value,
    normalize_emit,
)
from pyemittance.machine_settings import get_twiss0, get_rmat, get_energy, get_quad_len
from pyemittance.saving_io import save_emit_run
from pyemittance.load_json_configs import load_configs
from pyemittance.tools import mkdir_p

import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)

class EmitCalc:
    """
    Uses info recorded in Observer to do an emittance fit

    quad_vals: dict of quad values in machine units:
        kL = 


    """

    def __init__(
        self,
        quad_vals: dict = None,
        beam_vals: dict = None,
        beam_vals_err: dict = None,
        config_dict: dict = None, 
    ):

        self.quad_vals = (
            {
                "x": np.empty(
                    0,
                ),
                "y": np.empty(
                    0,
                ),
            }
            if quad_vals is None
            else quad_vals
        )  # in kG
        self.beam_vals = (
            {
                "x": np.empty(
                    0,
                ),
                "y": np.empty(
                    0,
                ),
            }
            if beam_vals is None
            else beam_vals
        )

        # Make sure error is added to beamsizes if none is provided
        if (
            beam_vals_err is None
            or sum(beam_vals_err["x"]) == 0
            or sum(beam_vals_err["y"]) == 0
        ):
            self.bs_error = (
                0.015,
                0.015,
            )  # Define some error on beamsizes in each dimension
            self.beam_vals_err = {
                "x": np.asarray(self.beam_vals["x"]) * self.bs_error[0],
                "y": np.asarray(self.beam_vals["y"]) * self.bs_error[1],
            }
        else:
            self.beam_vals_err = beam_vals_err

        if config_dict is None:
            raise ValueError("Must provide config_dict to EmitCalc class")

        self.config_dict = config_dict

        self.dims = [
            "x",
            "y",
        ]  # TODO: make code use provided in self.dims, and make it extensible
        self.sig_mat_screen = {"x": [], "y": []}
        self.twiss0 = get_twiss0(self.config_dict["beamline_info"])  # emit, beta, alpha
        self.twiss_screen = {"x": [], "y": []}  # emit, beta, alpha
        self.beta_err = None
        self.alpha_err = None
        self.energy = get_energy(self.config_dict["beamline_info"])
        self.rmat = get_rmat(self.config_dict["beamline_info"])
        self.quad_len = get_quad_len(self.config_dict["beamline_info"])

        self.calc_bmag = False
        self.plot = False
        self.save_runs = False
        
        # Initialize paths and dirs for saving
        ### self.init_saving()

        # Main output of emittance calc
        self.output = {}

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

    def error_propagation(self, gradient):
        """
        Propagate error from var(y) to emittance
        :param gradient: gradient of emittance
        :return: error on emittance from fit
        """
        return np.sqrt((gradient.T @ self.covariance_matrix) @ gradient)

    def get_emit(self):
        """
        Get emittance at quad from beamsizes and quad scan
        :return: normalized emittance and error
        """

        for dim in self.dims:
            # run emit calc for x and y

            q = np.array(self.quad_vals[dim])
            # quad vals are passed in machine units
            kL = kL_from_machine_value(q, self.energy)

            bs = self.beam_vals[dim]
            bs_err = self.beam_vals_err[dim]

            weights = self.weighting_func(bs, bs_err)  # 1/sigma

            # Storing quadvals and beamsizes in self.output for plotting purposes
            self.output[f"quadvals{dim}"] = np.array(q)
            self.output[f"beamsizes{dim}"] = np.array(bs)
            self.output[f"beamsizeserr{dim}"] = np.array(bs_err)

            res = estimate_sigma_mat_thick_quad(
                bs,
                kL,
                bs_err,
                weights,
                dim=dim,
                Lquad=self.quad_len,
                energy=self.energy,
                rmats=self.rmat,
                plot=self.plot,
            )
            # Add all results
            self.output.update(res)

            # Skip further calcs if there was an error
            if res[f"error_{dim}"]:
                continue

            if self.calc_bmag:
                if dim == "x":
                    sig_11 = res["screen_sigma_11"]
                    sig_12 = res["screen_sigma_12"]
                    sig_22 = res["screen_sigma_22"]

                else:
                    sig_11 = res["screen_sigma_33"]
                    sig_12 = res["screen_sigma_34"]
                    sig_22 = res["screen_sigma_44"]
                self.sig_mat_screen[dim] = [sig_11, sig_12, sig_22]

                beta_rel_err = res[f"beta_{dim}_rel_err"]
                alpha_rel_err = res[f"alpha_{dim}_rel_err"]

                self.beta_err = beta_rel_err
                self.alpha_err = alpha_rel_err

                bmag_calc_res = self.get_twiss_bmag(dim=dim)
                # Get bmag and bmag_err
                self.output[f"screen_bmag{dim}"] = bmag_calc_res[0]
                self.output[f"screen_bmag{dim}_err"] = bmag_calc_res[1]
                # Get best value for scanning quad
                self.output[f"optimal_quadval_{dim}"] = q[bmag_calc_res[2]]

            if self.plot:
                self.plot_output(dim=dim)

        # get geometric mean if possible
        if (not self.output["error_x"]) and (not self.output["error_y"]):
            self.get_gmean_emit()

        if self.save_runs:
            self.save_run()

        return self.output

    def plot_output(self, dim='x'):
        output = self.output
        fig, ax = plt.subplots(figsize=(8,2))
        x = output[f"quadvals{dim}"]
        if dim == 'x':
            suffix = '11'
            color = 'black'
        else:
            suffix = '33'
            color = 'red'
        scale = 1e6 # m -> um 
       # y1err = output[f'beamsizeserr{dim}'] * scale 
        ax.plot(x, np.sqrt(output[f'screen_sigma_{suffix}']) * scale, label=f'Model', color=color)
        ax.scatter(x, output[f'beamsizes{dim}'] * scale , marker='x', label=f'Measured', color=color)
        #ax.errorbar(x, y1, yerr=y1err, fmt='o', label=f'Measured', color=color)
        ax.set_xlabel('Quad Value (kG)')
        ax.set_ylabel(fr'Beam Size {dim} ($\mu$m)')
        plt.legend()


    def summary(self):
        data = self.output
        s = f""""
    Emittance Calculation Summary
    
    Emittance x: {data['norm_emit_x']/1e-6:.3f} +/- {data['norm_emit_x_err']/1e-6:.3f} mm mrad
    Emittance y: {data['norm_emit_y']/1e-6:.3f} +/- {data['norm_emit_y_err']/1e-6:.3f} mm mrad
    
    Before scanning quad:
                    x        y
    norm_emit  {data['norm_emit_x']/1e-6:8.2f}  {data['norm_emit_y']/1e-6:8.2f} (mm-mrad)                    
    beta       {data['beta_x']:8.2f}  {data['beta_y']:8.2f} (m)
    alpha      {data['alpha_x']:8.2f}  {data['alpha_y']:8.2f} (1)

    
    """
        return s





    def get_twiss_bmag(self, dim="x"):

        sig_11 = self.sig_mat_screen[dim][0]
        sig_12 = self.sig_mat_screen[dim][1]
        sig_22 = self.sig_mat_screen[dim][2]

        # twiss0 in x or y AT THE SCREEN
        beta0, alpha0 = self.twiss0[dim][1], self.twiss0[dim][2]

        # return dict of emit, beta, alpha, bmag
        twiss = twiss_and_bmag(
            sig_11,
            sig_12,
            sig_22,
            self.beta_err,
            self.alpha_err,
            beta0=beta0,
            alpha0=alpha0,
        )
        # Save twiss at screen
        self.twiss_screen[dim] = [twiss["emit"], twiss["beta"], twiss["alpha"]]

        return twiss["bmag"], twiss["bmag_err"], twiss["min_idx"]

    def get_gmean_emit(self):

        try:
            nemit = np.sqrt(self.output["norm_emit_x"] * self.output["norm_emit_y"])
            nemit_err = (
                nemit
                * (
                    (self.output["norm_emit_x_err"] / self.output["norm_emit_x"]) ** 2
                    + (self.output["norm_emit_y_err"] / self.output["norm_emit_y"]) ** 2
                )
                ** 0.5
            )

            self.output["sqrt_norm_emit_4d"] = nemit
            self.output["sqrt_norm_emit_4d_err"] = nemit_err

            if "bmag_x" in self.output and "bmag_y" in self.output:
                nbmag = np.sqrt(self.output["bmag_x"] * self.output["bmag_y"])
                bmag_emit_err = (
                    nemit
                    * nbmag
                    * (
                        (self.output["norm_emit_x_err"] / self.output["norm_emit_x"])
                        ** 2
                        + (self.output["norm_emit_y_err"] / self.output["norm_emit_y"])
                        ** 2
                        + (self.output["bmag_x_err"] / self.output["bmag_x"]) ** 2
                        + (self.output["bmag_y_err"] / self.output["bmag_y"]) ** 2
                    )
                    ** 0.5
                )
                self.output["bmag_emit"] = nemit * nbmag
                self.output["bmag_emit_err"] = bmag_emit_err

        except TypeError:
            self.output["sqrt_norm_emit_4d"] = np.nan
            self.output["sqrt_norm_emit_4d_err"] = np.nan
            self.output["bmag_emit"] = np.nan
            self.output["bmag_emit_err"] = np.nan

    def save_run(self):
        save_emit_run(self.output, path=self.config_dict["savepaths"]["fits"])

    def init_saving(self):
        """Initialize dirs and files for saving"""

        savepaths = self.config_dict["savepaths"]

        # Make directories if needed
        try:
            mkdir_p(savepaths["images"])
            mkdir_p(savepaths["summaries"])
            mkdir_p(savepaths["fits"])
            mkdir_p(savepaths["raw_saves"])
        except OSError:
            logger.info("Savepaths not set. Please set them in 'configs/savepaths.json'")

            parent = Path(__file__).resolve().parent
            examples_dir = str(parent)[:-11] + "examples"
            logger.info("Using examples directory: ", examples_dir)
            savepaths["images"] = examples_dir + "/saved_images/"
            savepaths["summaries"] = examples_dir + "/summaries/"
            savepaths["fits"] = examples_dir + "/saved_fits/"
            savepaths["raw_saves"] = examples_dir + "/raw_saves/"
            mkdir_p(savepaths["images"])
            mkdir_p(savepaths["summaries"])
            mkdir_p(savepaths["fits"])
            mkdir_p(savepaths["raw_saves"])

        # Start headings
        file_exists = os.path.exists(savepaths["summaries"] + "image_acq_quad_info.csv")

        if not file_exists:

            # TODO: add others as inputs
            f = open(savepaths["summaries"] + "image_acq_quad_info.csv", "a+")
            f.write(
                f"{'timestamp'},{'ncol'},{'nrow'},{'roi_xmin'},{'roi_xmax'}"
                f",{'roi_ymin'},{'roi_ymax'},{'resolution'},{'bact'},"
                f"{'x_size'},{'y_size'},{'xrms'},{'yrms'},"
                f"{'xrms_err'},{'yrms_err]'}\n"
            )
            f.close()

        file_exists = os.path.exists(savepaths["summaries"] + "beamsize_config_info.csv")

        if not file_exists:
            # todo add others as inputs
            f = open(savepaths["summaries"] + "beamsize_config_info.csv", "a+")
            f.write(
                f"{'timestamp'},{'varx_cur'},{'vary_cur'},{'varz_cur'},"
                f"{'bact_cur'},{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'}\n"
            )
            f.close()
