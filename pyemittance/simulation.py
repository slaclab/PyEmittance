import numpy as np

from pyemittance.optics import quad_mat2, propagate_sigma, sigma_from_twiss, mec2
from pyemittance.optics import machine_value_from_kL, kL_from_machine_value

import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)

BUNCH_PARAMS = {
    'LCLS2_OTR0H04': {
        'total_charge': 50e-12,
        'norm_emit_x': 1e-6,
        'norm_emit_y': 2e-6,
        'beta_x': 10,
        'alpha_x': -2,
        'beta_y': 11,
        'alpha_y': 20,
        'energy': 80e6,
        'species':'electron'
    }
}

SCREEN_PARAMS = {
    'LCLS2_OTR0H04': {
     'nrow':1040,
     'ncol':1392,
     'resolution': 20.2e-6,
     'noise': 10,
    }    
}

DEFAULT_CONFIG_NAME = 'LCLS2_OTR0H04'
DEFAULT_BUNCH_PARAMS = BUNCH_PARAMS[DEFAULT_CONFIG_NAME]
DEFAULT_SCREEN_PARAMS = SCREEN_PARAMS[DEFAULT_CONFIG_NAME]



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
             total_charge=100e-12, n_particle=10_000,
            brightness = 1e17,
            ):
        
        """
        Generate a spot due to a beam of particles.
        
        brightness: enhancement factor
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
        H *= brightness 
        
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
        
        # Internal state
        self._quad_value = 0.0  # machine units
        self._beamon = True # Beam on
        
        self.output = {}
    
    @property
    def beamon(self):
        return self._beamon
    @beamon.setter
    def beamon(self, value):
        self._beamon = bool(value)
        
    @property
    def quad_value(self):
        return self._quad_value  
    @quad_value.setter
    def quad_value(self, value):
        self._quad_value = value
        
    @property
    def energy(self):
        return self.beamline_info['energy']
    @energy.setter
    def energy(self, value):
        self.beamline_info['energy'] = value
    @property
    def rmatx(self):
        return np.array(self.beamline_info['rMatx']).reshape(2,2)
    @rmatx.setter
    def rmatx(self, value):
        self.beamline_info['rMatx'] = np.array(value).reshape(2,2)    
    @property
    def rmaty(self):
        return np.array(self.beamline_info['rMaty']).reshape(2,2)
    @rmaty.setter
    def rmaty(self, value):
        self.beamline_info['rMaty'] = np.array(value).reshape(2,2)            
    
    @property 
    def Lquad(self):
        return self.beamline_info['Lquad']
    @Lquad.setter 
    def Lquad(self, value):
        self.beamline_info['Lquad'] = value

    def initial_sigma_matrix2(self, dim=''):
        norm_emit = self.bunch_params[f'norm_emit_{dim}']
        gammabeta = np.sqrt((self.bunch_params['energy']/mec2)**2 -1)
        emit = norm_emit / gammabeta
        beta = self.bunch_params[f'beta_{dim}']
        alpha = self.bunch_params[f'alpha_{dim}']
        return sigma_from_twiss(emit, beta, alpha)
    
    def screen_sigma(self, dim='x'):
        kL = kL_from_machine_value(self.quad_value, self.energy)
        sigma0 = self.initial_sigma_matrix2(dim)
        if dim == 'x':
            sign = 1
            rmat = self.rmatx
        elif dim == 'y':
            sign = -1
            rmat = self.rmaty
        else:
            raise ValueError(f"dim = {dim} not in ('x', 'y')")
        mat2 = rmat @ quad_mat2(sign*kL, L=self.Lquad)
        sigma1 = propagate_sigma(sigma0, mat2)
        meas_sigma =  np.sqrt(sigma1[0,0])
        
        return meas_sigma
      
    def screen_beam_sizes(self):
        return self.screen_sigma('x'), self.screen_sigma('y')
    
    def beam_size_meas(self, quad_value):
        self.quad_value = quad_value
        return self.screen_beam_sizes()
    
    
    def screen_image(self):
        S = self.screen
        bg = S.background()
        if not self.beamon:
            return bg
        
        
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
        
        
        
        
        