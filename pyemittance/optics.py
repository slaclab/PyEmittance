# Module containing functions for beam optics calculations
import numpy as np
from numpy import sin, cos, sinh, cosh, sqrt
import scipy.linalg

from scipy.constants import c as c_light 
mec2 = scipy.constants.value('electron mass energy equivalent in MeV')*1e6

import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)


def kL_from_machine_value(quad_val, energy):
    """
    Calculates quadrupole strength from machine value.

    Parameters
    ----------
    quad_val: ndarray
        integrated quad field -[kG] (SLAC convention)

    energy: float
        Beam energy [eV]

    Returns
    -------
    ndarray
    Quad strength  [1/m^2]

    """
    gamma = energy / mec2
    beta = np.sqrt(1 - 1 / gamma**2)

    return c_light*quad_val * 0.1 / (beta * energy) # 1/m^2

def machine_value_from_kL(kL, energy):
    """
    Calculates machine value from quadrupole strength.

    Parameters
    ----------
    kL: ndarray
        integrated quad focusing strength [1/m]

    energy: float
        Beam energy [eV]

    Returns
    -------
    ndarray
    Integrated quad field [kG] (SLAC positron convention)

    """
    gamma = energy / mec2
    beta = np.sqrt(1 - 1 / gamma**2)
    return kL * beta * energy / c_light * 10 # 1/m -> kG







def sigma_from_twiss(emit, beta, alpha):
    """
    Returns the 2x2 covariance matrix from Twiss parameters
    """
    gamma = (1+alpha**2)/beta
    sigma_x = np.sqrt(emit*beta)
    sigma_px = np.sqrt(emit*gamma)
    cov_x__px = -alpha*emit
    sigma_11 = sigma_x**2
    sigma_12 = cov_x__px
    sigma_22 = sigma_px**2
    sigma0 = np.array([[sigma_11, sigma_12], [sigma_12, sigma_22]])
    return sigma0


def normalize_emit(emit, energy, species="electron"):
    assert species == "electron", "Only electron species supported"
    gamma = energy / mec2
    beta = np.sqrt(1 - 1 / gamma**2)
    return emit * gamma * beta

def thin_quad_mat2(kL):
    """
    Quad transport matrix, 2x2, assuming thin quad
    :param kL: quad strength * quad length (1/m)
    :return: thin quad transport matrix
    """
    return np.array([[1, 0], [-kL, 1]])


def drift_mat2(L):
    """
    2x2 Transport matrix for a simple drift with length L
    """
    return np.array([[1, L], [0, 1]])


def quad_mat2(kL, L=0):
    """
    Quadrupole transfer matrix, 2x2, assuming some quad thickness
    L = 0 returns thin quad matrix
    d not None get drift mat
    :param kL: quad strength * quad length (1/m)
    :param L: quad length (m)
    :return: thick quad transport matrix
    """

    if L == 0:
        return thin_quad_mat2(kL)

    k = kL / L

    if k == 0:
        # Take drift or custom mat
        mat2 = drift_mat2(L)
    elif k > 0:
        # Focusing
        rk = sqrt(k)
        phi = rk * L
        mat2 = [[cos(phi), sin(phi) / rk], [-rk * sin(phi), cos(phi)]]
    else:
        # Defocusing
        rk = sqrt(-k)
        phi = rk * L
        mat2 = [[cosh(phi), sinh(phi) / rk], [rk * sinh(phi), cosh(phi)]]

    return mat2


def propagate_sigma(mat2_init, mat2_ele):
    """
    Propagate a transport matrix through beamline from point A to B
    :return: 2x2 matrix at B
    """
    return (mat2_ele @ mat2_init) @ mat2_ele.T


def estimate_sigma_mat_thick_quad(
    sizes,
    kLlist,
    sizes_err=None,
    weights=None,
    dim="x",
    Lquad=None,
    energy=None,
    rmats=None,
    calc_bmag=False,
    plot=True,
):
    """
    Estimates the beam sigma matrix at a screen by scanning an upstream quad.
    This models the system as a thick quad.
    :param sizes: measured beam sizes at the screen
    :param kLlist: kL of the upstream quad
    :param weights: If present, each measurement will be weighted by this factor.
    :param Lquad:  length of the quadrupole magnet (m)
    :param plot: bool to plot ot not
    :return: emittance, sig11, sig12 and sig22 at measurement screen
    """

    if dim not in ("x", "y"):
        raise ValueError(f"Bad dim: {dim}")

    # Measurement vector
    sizes = np.array(sizes)
    if np.isnan(sizes).any():
        idx_not_nan = ~np.isnan(sizes)
        sizes = sizes[idx_not_nan]
        kLlist = np.array(kLlist)[idx_not_nan]
        sizes_err = np.array(sizes_err)[idx_not_nan]
        if weights is not None:
            weights = np.array(weights)
            weights = weights[idx_not_nan]

    if len(sizes) < 3:
        logger.warning("Less than 3 data points were passed. Returning NaN.")
        return [np.nan, np.nan, np.nan, np.nan]

    b = sizes**2
    n = len(b)

    # Fill in defaults, checking.
    if weights is None:
        weights = np.ones(n)
    assert len(weights) == n

    # Get rmat from configs
    r_mat = rmats[0] if dim == "x" else rmats[1]

    # Multiply by weights. This should correspond to the other weight multiplication below
    b = weights * sizes**2

    # form B matrix
    B = []
    # Collect mat2 for later
    mat2s = []
    for kL, weight in zip(kLlist, weights):
        if dim == 'x':
            mat2 = r_mat @ quad_mat2( kL, L=Lquad)
        else:
            mat2 = r_mat @ quad_mat2(-kL, L=Lquad)
        mat2s.append(mat2)
        r11, r12, r21, r22 = mat2.flatten()
        r_mat_factor = np.array([r11**2, 2 * r11 * r12, r12**2])
        B.append(r_mat_factor * weight)  # corresponding weight multiplication

    B = np.array(B)

    # Invert (pseudoinverse)
    s11, s12, s22 = scipy.linalg.pinv(B) @ b

    # Twiss calc just before the quad
    emit2 = s11 * s22 - s12**2

    # return NaN if emit can't be calculated
    if emit2 < 0:
        logger.error("Emittance can't be computed. Returning error")
        return {f"error_{dim}": True}

    emit = np.sqrt(emit2)
    beta = s11 / emit
    alpha = -s12 / emit

    # Get error on emittance from fitted params
    emit_err, beta_err, alpha_err = get_twiss_error(emit, s11, s12, s22, B)

    # Normalized emittance
    norm_emit = normalize_emit(emit, energy)
    norm_emit_err = normalize_emit(emit_err, energy)

    # Propagate to screen
    s11_screen, s12_screen, s22_screen = propagate_to_screen(
        s11, s12, s22, kLlist, mat2s, Lquad, energy, sizes, sizes_err,
    )

    out = {}
    out[f"error_{dim}"] = False
    out[f"emit_{dim}"] = emit
    out[f"norm_emit_{dim}"] = norm_emit
    out[f"beta_{dim}"] = beta
    out[f"alpha_{dim}"] = alpha
    out[f"emit_{dim}_err"] = emit_err
    out[f"norm_emit_{dim}_err"] = norm_emit_err
    out[f"beta_{dim}_rel_err"] = beta_err / beta
    out[f"alpha_{dim}_rel_err"] = alpha_err / alpha  # Note: alpha can be zero!

    # Sigma matrix info
    if dim == "x":
        out["sigma_11"] = s11
        out["sigma_12"] = s12
        out["sigma_22"] = s22
        out["screen_sigma_11"] = s11_screen
        out["screen_sigma_12"] = s12_screen
        out["screen_sigma_22"] = s22_screen
    elif dim == "y":
        out["sigma_33"] = s11
        out["sigma_34"] = s12
        out["sigma_44"] = s22
        out["screen_sigma_33"] = s11_screen
        out["screen_sigma_34"] = s12_screen
        out["screen_sigma_44"] = s22_screen
    else:
        raise ValueError(f"Bad dim: {dim}")

    return out


def propagate_to_screen(
    s11, s12, s22, kLlist, mat2s, Lquad, energy, sizes, sizes_err,
    ):
    # Matrix form for propagation
    sigma0 = np.array([[s11, s12], [s12, s22]])

    # Propagate forward to the screen
    s11_screen = []
    s12_screen = []
    s22_screen = []
    for kL, mat2 in zip(kLlist, mat2s):
        sigma1 = propagate_sigma(sigma0, mat2)
        s11_screen.append(sigma1[0, 0])
        s12_screen.append(sigma1[0, 1])
        s22_screen.append(sigma1[1, 1])
    s11_screen = np.array(s11_screen)
    s12_screen = np.array(s12_screen)
    s22_screen = np.array(s22_screen)

    return s11_screen, s12_screen, s22_screen


def twiss_and_bmag(sig11, sig12, sig22, beta_err, alpha_err, beta0=1, alpha0=0):
    """
    Calculates Twiss ang Bmag from the sigma matrix.
    """

    # Twiss at screen
    emit = np.sqrt(sig11 * sig22 - sig12**2)
    beta = sig11 / emit
    alpha = -sig12 / emit
    gamma = sig22 / emit

    # Form bmag
    gamma0 = (1 + alpha0**2) / beta0
    bmag = (beta * gamma0 - 2 * alpha * alpha0 + gamma * beta0) / 2
    # Add err in quadrature (assuming twiss0 has no uncertainty)
    # Taking relative error as measured at quad
    bmag_err = bmag * np.sqrt((beta_err) ** 2 + (alpha_err) ** 2)

    # TODO: make this take bmag at the nominal Q5 value and not the min to match Matlab
    # what to do if nominal isn't within scanned range?
    bmag_min = min(bmag)
    bmag_min_err = bmag_err[np.argmin(bmag)]

    d = {}
    d["emit"] = emit
    d["beta"] = beta
    d["alpha"] = alpha
    d["bmag"] = bmag_min
    d["bmag_err"] = bmag_min_err
    d["min_idx"] = np.argmin(bmag)  # best scanning quad val

    return d


def gradient_mat3(emit, a1, a2, a3):
    """
    Gradient of f = { emittance, beta, alpha }
    where f is obtained at the scanning location (quad)
    :param emit: emittance parameter estimate
    :param a1: matrix element s11
    :param a2: matrix element s12
    :param a3: matrix element s22
    :return: gradient of f
    """

    emit_gradient = 1.0 / (2 * emit) * np.array([a3, -2 * a2, a1])
    beta_gradient = (
        1.0
        / (2 * emit**3)
        * np.array([2 * emit**4 - a1 * a3, 2 * a2 * a1, -(a1**2)])
    )
    alpha_gradient = (
        -1.0 / (2 * emit) * np.array([a2 * a3, 2 * emit**2 - 2 * a2**2, a1 * a2])
    )

    f_gradient = np.array([emit_gradient, beta_gradient, alpha_gradient]).T

    return f_gradient


def get_fit_param_error(f_gradient, B):
    """
    Error estimation of the fitted params (s11, s12, s22)
    and the Twiss params. See 10.3204/DESY-THESIS-2005-014 p.10
    :param f_gradient: gradient of the 3-vector Twiss params
    :param B: B matrix incl. all scanned quad values
    :return: sqrt of the diagonal of the error matrix (emit_err, beta_err, alpha_err)
    """

    C = scipy.linalg.pinv(B.T @ B)

    error_matrix = f_gradient.T @ C @ f_gradient
    twiss_error = np.sqrt(np.diag(error_matrix))

    return twiss_error


def get_twiss_error(emit, a1, a2, a3, B):
    """
    Get error on the twiss params from fitted params
    :param emit: emittance parameter estimate
    :param a1: matrix element s11
    :param a2: matrix element s12
    :param a3: matrix element s22
    :param B:  B matrix incl. all scanned quad values
    :return: emit_err, beta_err, alpha_err
    """

    # get gradient of twiss params
    f_gradient = gradient_mat3(emit, a1, a2, a3)
    # calculate errors on twiss from var and covar
    twiss_error = get_fit_param_error(f_gradient, B)

    return twiss_error
