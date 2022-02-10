# Module containing functions for beam optics calculations
import numpy as np
from numpy import sin, cos, sinh, cosh, sqrt
import scipy.linalg
import matplotlib.pyplot as plt

def thin_quad_mat2(kL):
    '''
    Quad transport matrix, 2x2, assuming thin quad
    :param kL: quad strength * quad length (1/m)
    :return: thin quad transport matrix
    '''
    return np.array([[1, 0], [-kL, 1]])

def drift_mat2(L):
    '''
    Drift transport matrix, 2x2
    :param L: drift length (m)
    :return: drift transport matrix
    '''
    return np.array([[1, L], [0, 1]])

def quad_mat2(kL, L=0):
    '''
    Quadrupole transfer matrix, 2x2, assuming some quad thickness
    L = 0 returns thin quad matrix
    :param kL: quad strength * quad length (1/m)
    :param L: quad length (m)
    :return: thick quad transport matrix
    '''

    if L == 0:
        return thin_quad_mat2(kL)

    k = kL / L

    if k == 0:
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

def quad_drift_mat2(kL, *, Ltot=2.0, Lquad=0):
    '''
    Composite [quad, drift] 2x2 transfer matrix
    :param kL: quad strength * quad length (1/m)
    :param Ltot: length between quad and drift (m)
    :param Lquad: quad length (m)
    :return:
    '''

    Ldrift = Ltot - Lquad

    return drift_mat2(Ldrift) @ quad_mat2(kL, Lquad)

def propagate_sigma(mat2_init, mat2_ele):
    '''
    Propagate a transport matrix through beamline from point A to B
    :param sigma_mat2: 2x2 matrix at A
    :param mat2: total 2x2 trasport matrix of elements between point A and B
    :return: 2x2 matrix at B
    '''
    return (mat2_ele @ mat2_init) @ mat2_ele.T

def estimate_sigma_mat_thick_quad_drift(sizes, kLlist, weights=None, Ltot=2.26, Lquad=0.108, plot=True):
        '''
        Estimates the beam sigma matrix at a screen by scanning an upstream quad.
        This models the system as a thick quad.
        :param sizes: measured beam sizes at the screen
        :param kLlist: kL of the upstream quad
        :param weights: If present, each measurement will be weighted by this factor.
        :param Ltot: total length (m)
        :param Lquad:  length of the quadrupole magnet (m)
        :param plot: bool to plot ot not
        :return: emittance, sig11, sig12 and sig22 at measurement screen
        '''

        # measuerement vector
        b = sizes ** 2
        n = len(b)

        # Fill in defaults, checking.
        if weights is None:
            weights = np.ones(n)
        assert len(weights) == n

        # Multiply by weights. This should correspond to the other weight multiplication below
        b = weights * sizes ** 2

        # form B matrix
        B = []
        # Collect mat2 for later
        mat2s = []
        for kL, weight in zip(kLlist, weights):
            mat2 = quad_drift_mat2(kL, Ltot=Ltot, Lquad=Lquad)
            mat2s.append(mat2)
            r11, r12, r21, r22 = mat2.flatten()
            r_mat_factor = np.array([r11 ** 2, 2 * r11 * r12, r12 ** 2])
            B.append(r_mat_factor * weight)  # corresponding weight multiplication

        B = np.array(B)

        # Invert (pseudoinverse)
        s11, s12, s22 = scipy.linalg.pinv(B) @ b

        # Twiss calc just before the quad
        emit2 = s11 * s22 - s12 ** 2
        emit = np.sqrt(emit2)
        beta = s11 / emit
        alpha = -s12 / emit

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

        if plot:
            # Plot the data
            plt.scatter(kLlist, sizes ** 2, marker='x', label=f'Measurements')

            # Model prediction
            plt.scatter(kLlist, s11_screen, marker='.', label=f'Model')

            plt.title(f'beta={beta:.1f} m, alpha = {alpha:0.2f}, emit = {emit * 1e9:0.2f} nm')
            plt.xlabel('kL (1/m)')
            plt.ylabel(f'sizes^2 (m$^2$)')
            plt.ylim(0, None)
            plt.legend()

        return emit, s11_screen, s12_screen, s22_screen
