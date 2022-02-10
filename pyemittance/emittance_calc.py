import numpy as np
from numpy import sin, cos, sinh, cosh, sqrt
import matplotlib.pyplot as plt
import scipy.linalg

class EmitCalc:
    '''
    Uses info recorded in Observer to do an emittance fit

    '''

    def __init__(self, quad_vals=np.empty(0, ), beam_vals={'x': np.empty(0, ), 'y': np.empty(0, )}):
        self.quad_vals = quad_vals
        self.beam_vals = beam_vals
        self.beam_vals_err = {dim: beam_vals[dim]*0.1 for dim in beam_vals}
        self.x_use = np.arange(0, len(beam_vals['x']), 1)
        self.y_use = np.arange(0, len(beam_vals['y']), 1)

        self.sig_11 = []
        self.sig_12 = []
        self.sig_22 = []
        self.sig_quad = []
        self.covariance_matrix = None
        
        self.test_mode = False
        self.noise_red = 50000

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
        var_bs = ( 2 * beamsizes * beamsizes_err )**2
        weights = 1 / beamsizes + 1 / var_bs
        return weights

    def error_propagation(self, gradient):
        """
        Propagate error from var(y) to emittance
        :param gradient: gradient of emittance
        :return: error on emittance from fit
        """
        return np.sqrt( (gradient.T @ self.covariance_matrix) @ gradient)

    def get_emit_gradient(self, sig_ele, emit):
        """
        Gradient of emittance to calculate cov of parameters estimated
        :param s11: ME11 of r matrix at screen
        :param s12: ME12 of r matrix at screen
        :param s22: ME12 of r matrix at screen
        @param emit: emittance from fit
        :return: gradient of emittance
        """
        s11 = sig_ele[0]
        s12 = sig_ele[1]
        s22 = sig_ele[2]
        emit_gradient = 1/(2*emit) * np.array([[s11, -2*s12, s22]]).T
        return emit_gradient

    def do_emit_fit(self, dim='x'):

        # todo update based on x_use, y_use for throwing away fit points
        q = self.quad_vals

        bs = self.beam_vals[dim]
        bs_err = self.beam_vals_err[dim]

        weights = self.weighting_func(bs, bs_err)

        if self.test_mode == False:
            emit, self.sig_11, self.sig_12, self._sig_22 = self.estimate_sigma_mat_thick_quad_drift(bs, q, weights)
            plt.show()

        if self.test_mode == True:
            bs = bs + np.random.rand(len(bs)) / self.noise_red
            print("NOISE")
            emit, self.sig_11, self.sig_12, self._sig_22 = self.estimate_sigma_mat_thick_quad_drift(bs, q, weights)
            plt.show()

        err = np.std(np.absolute(self.sig_11 - bs))
        
#         print(emit, err)
#         err_cov = self.error_propagation(self.get_emit_gradient(self.sig_quad, emit))

#         print(emit, err, err_cov)
        
        return emit, err

    def thin_quad_mat2(self, kL):
        return np.array([[1, 0], [-kL, 1]])

    def drift_mat2(self, L):
        return np.array([[1, L], [0, 1]])

    def quad_mat2(sel, kL, L=0):
        """
        Quadrupole transfer matrix, 2x2. Note that

        """

        if L == 0:
            return self.thin_quad_mat2(kL)

        k = kL / L

        if k == 0:
            mat2 = self.drift_mat2(L)
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

    def quad_drift_mat2(self, kL, *, Ltot=2.0, Lquad=0):
        """
        Composite [quad, drift] 2x2 transfer matrix.
        """

        Ldrift = Ltot - Lquad

        return self.drift_mat2(Ldrift) @ self.quad_mat2(kL, Lquad)

    def propagate_sigma(self, sigma_mat2, mat2):
        return (mat2 @ sigma_mat2) @ mat2.T

    def estimate_sigma_mat_thick_quad_drift(self, sizes, kLlist, weights=None, Ltot=2.26, Lquad=0.108, plot=True):
        """
        Estimates the beam sigma matrix at a screen by scanning an upstream quad.

        This models the system as a thick quad.

        Parameters
        ----------
        sizes : array of float
            measured beam sizes at the screen

        kLlist : array of float
            kL of the upstream


        weights : array of float or None
            If present, each measurement will be weighted by this factor.

        Ltot : float
            total length in meters

        Lquad: float
            Length of the quadrupole magnet in meters

        plot : bool

        Returns
        -------

        s11_screen : array of float

        s12_screen : array of float

        s22_screen : tuple of float

        """

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
        self.covariance_matrix = []
        
        # Collect mat2 for later
        mat2s = []
        for kL, weight in zip(kLlist, weights):
            mat2 = self.quad_drift_mat2(kL, Ltot=Ltot, Lquad=Lquad)
            mat2s.append(mat2)
            r11, r12, r21, r22 = mat2.flatten()
            r_mat_factor = np.array([r11 ** 2, 2 * r11 * r12, r12 ** 2]) 
            B.append(r_mat_factor * weight)  # corresponding weight multiplication
            
#         w_diag = np.eye(len(weights))*weights
#         self.covariance_matrix = scipy.linalg.pinv( 
#             r_mat_factor @ w_diag @ r_mat_factor.T
#         )
                                     
        B = np.array(B)

        # Invert (pseudoinverse)
        s11, s12, s22 = scipy.linalg.pinv(B) @ b

        # Twiss calc just before the quad
        emit2 = s11 * s22 - s12 ** 2
        emit = np.sqrt(emit2)
        beta = s11 / emit
        alpha = -s12 / emit
        
        self.sig_quad = [s11, s12, s22]

        # Matrix form for propagation
        sigma0 = np.array([[s11, s12], [s12, s22]])

        # Propagate forward to the screen
        s11_screen = []
        s12_screen = []
        s22_screen = []
        for kL, mat2 in zip(kLlist, mat2s):
            sigma1 = self.propagate_sigma(sigma0, mat2)
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