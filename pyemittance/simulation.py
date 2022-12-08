import numpy as np

from pyemittance.optics import quad_mat2, propagate_sigma, sigma_from_twiss, mec2
from pyemittance.optics import machine_value_from_kL, kL_from_machine_value

import matplotlib.pyplot as plt


def generate_particles(n_particles, x_mean=0, x_std=1, y_mean=0, y_std=1, x_y_cov=0):
    """Generate a set of particles with a given mean and standard deviation in x and y,
    and a given covariance between x and y.

    Parameters
    ----------
    n_particles : int
        Number of particles to generate
    x_mean : float
        Mean of x distribution
    x_std : float
        Standard deviation of x distribution
    y_mean : float
        Mean of y distribution
    y_std : float
        Standard deviation of y distribution
    x_y_cov : float
        Covariance between x and y

    Returns
    -------
    x : ndarray
        x positions of particles
    y : ndarray
        y positions of particles
    """
    x, y = np.random.multivariate_normal(
        mean=[x_mean, y_mean], cov=[[x_std ** 2, x_y_cov], [x_y_cov, y_std ** 2]], size=n_particles
    ).T
    return x, y

class Screen:
    def __init__(self, nrow=1040,
                 ncol=1392,
                 resolution = 4e-6,
                 dtype = np.int16,
                 noise = 10):
                
        self.nrow = nrow
        self.ncol = ncol
        self.resolution = resolution
        self.noise = noise
        self.dtype=dtype
    """
    Simulated screen
    """
        
    @property
    def width(self):
        return self.ncol * self.resolution
    @property
    def height(self):
        return self.nrow * self.resolution    
    @property
    def xmin(self):
        return -self.width/2
    @property
    def xmax(self):
        return  self.width/2 
    @property
    def ymin(self):
        return -self.height/2
    @property
    def ymax(self):
        return  self.height/2      
        
    def background(self): 
        im = np.random.rand(self.nrow, self.ncol) * self.noise
        return im.astype(self.dtype)
    
    def spot(self, sigma_x=100e-6, sigma_y=200e-6, 
             mean_x = 0,
             mean_y = 0,
             total_charge=100e-12, n_particle=10_000):
        
        """
        Generate a spot due to a beam of particles.
        """
        x, y = generate_particles(n_particle, x_std=sigma_x,
                                  y_std=sigma_y,
                                  x_mean = mean_x,
                                  y_mean = mean_y,
                                  x_y_cov=0)  
        w = np.full(n_particle, total_charge/n_particle)
        H, xedges, yedges = np.histogram2d(x, y, weights=w,
                                 range = ( (self.xmin, self.xmax), 
                                           (self.ymin, self.ymax)
                                         ),
                                 bins=(self.ncol, self.nrow),
                                        density=False,
                                         )        
        H *= 1e17 # TEMP
        
        return (np.flipud(H.T)).astype(self.dtype)
    
    
    
class BeamSim:
    def __init__(self,
                 bunch_params=None,
                 beamline_info=None,
                 screen=None) -> None:
        self.bunch_params = bunch_params.copy()
        self.beamline_info = beamline_info.copy()
        if screen is None:
            screen = Screen() # use defaults
        self.screen = screen
        
        self._quad_value = 0.0  # machine units
        self.configure()
        
    @property
    def quad_value(self):
        return self._quad_value
    
    @quad_value.setter
    def quad_value(self, value):
        self._quad_value = value
        
        
    def configure(self):
        self.energy = self.beamline_info['energy']
        self.rmatx = np.array(self.beamline_info['rMatx']).reshape(2,2)
        self.rmaty = np.array(self.beamline_info['rMaty']).reshape(2,2)
        self.Lquad = self.beamline_info['Lquad']

    def initial_sigma_matrix2(self, dim=''):
        norm_emit = self.bunch_params[f'norm_emit_{dim}']
        gammabeta = np.sqrt((self.bunch_params['energy']/mec2)**2 -1)
        emit = norm_emit / gammabeta
        beta = self.bunch_params[f'beta_{dim}']
        alpha = self.bunch_params[f'alpha_{dim}']
        return sigma_from_twiss(emit, beta, alpha)
    
    def screen_beam_sizes(self):
        kL = kL_from_machine_value(self.quad_value, self.energy)
        
        # X
        sigma0 = self.initial_sigma_matrix2('x')
        mat2 = self.rmatx @ quad_mat2(kL, L=self.Lquad)
        sigma1 = propagate_sigma(sigma0, mat2)
        meas_sigma_x =  np.sqrt(sigma1[0,0])
        # Y
        sigma0 =  self.initial_sigma_matrix2('y')
        mat2 = self.rmaty @ quad_mat2(-kL, L=self.Lquad)
        sigma1 = propagate_sigma(sigma0, mat2) 
        meas_sigma_y =  np.sqrt(sigma1[0,0])
        return meas_sigma_x, meas_sigma_y  
    
    def beam_size_meas(self, quad_value):
        self.quad_value = quad_value
        return self.screen_beam_sizes()
    
    
    def screen_image(self):
        S = self.screen
        bg = S.background()
        sigma_x, sigma_y = self.screen_beam_sizes()
        total_charge = self.bunch_params['total_charge']
        im = S.spot(sigma_x=sigma_x,
                    sigma_y=sigma_y,
                    n_particle=100_000,
                    total_charge=total_charge)
        return im + bg
    
    def plot_screen(self, return_figure=False, vmax=128):
        """
        
        """
        fig, ax = plt.subplots()
        im = self.screen_image()
        S = self.screen
        ax.imshow(im, extent=1e3*np.array([S.xmin, S.xmax, S.ymin, S.ymax]), vmax=vmax)
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_aspect('auto')
        if return_figure:
            return figure
        
        
        
        