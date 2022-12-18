from pcaspy import Driver, SimpleServer

from pyemittance.load_json_configs import load_configs
from pyemittance.optics import drift_mat2
from pyemittance.simulation import BeamSim, Screen, BUNCH_PARAMS, SCREEN_PARAMS
import epics

import logging
logger = logging.getLogger(__name__)

def get_all_params(config_name='LCLS2_OTR0H04'):
    assert config_name == 'LCLS2_OTR0H04'
        
    config_dict = load_configs(config_name)
    # Conveniences
    beamline_info = config_dict['beamline_info']
    meas_pv_info = config_dict['meas_pv_info']
    

    #beamline_info['rMatx'] = drift_mat2(2.2)
    #beamline_info['rMaty'] = drift_mat2(2.2)
    #beamline_info['Lquad'] = 0.108
    
    bunch_params = BUNCH_PARAMS[config_name]
    screen_params = SCREEN_PARAMS[config_name]

    
    pvmap = {
     'quadval': meas_pv_info['meas_device']['pv']['cntrl'],
     'quadval_rbv': meas_pv_info['meas_device']['pv']['read'],
     'image_array': meas_pv_info['diagnostic']['pv']['image'],
     'nrow': meas_pv_info['diagnostic']['pv']['nrow'],
     'ncol': meas_pv_info['diagnostic']['pv']['ncol'],
     'resolution': meas_pv_info['diagnostic']['pv']['resolution'],    
    }
    
    return bunch_params, screen_params, pvmap, beamline_info
    
def make_pvdb(pvmap, screen_params):
    pvdb = {
        pvmap['quadval']: {
            'value': 0.0,
            'prec' : 5,
        },
        pvmap['quadval_rbv']: {
            'value': 0.0,
            'prec' : 5,
        },   
        pvmap['image_array']: {
            'type': 'int',
            'count': screen_params['nrow'] *  screen_params['ncol']
        },    
        pvmap['nrow']: {
            'type': 'int',
            'value': screen_params['nrow']
        },
        pvmap['ncol']: {
            'type': 'int',
            'value': screen_params['ncol'],
        },
        pvmap['resolution']: {
            'value': screen_params['resolution'],
            'unit': 'um',
        },
        # Simulation values
        'sim_screen_sigma_x': {
            'value': 0.0,
            'prec': 5,
        },
        'sim_screen_sigma_y': {
            'value': 0.0
        },        
    }    
    return pvdb
    
    
    
class BeamSimDriver(Driver):
    def __init__(self, beamsim, pvmap):
        super().__init__()
        self.sim = beamsim
        self.pvmap = pvmap

    def read(self, reason):  
        logger.debug(f'BeamSimDriver read reason: {reason}')
        if reason in (self.pvmap['quadval_rbv'], self.pvmap['quadval']):
             value = self.sim.quad_value          
        elif reason == self.pvmap['image_array']:
            value = list(self.sim.screen_image().flatten())
        elif reason == self.pvmap['nrow']:
            value = self.sim.screen.nrow
        elif reason == self.pvmap['ncol']:
            value = self.sim.screen.ncol          
        elif reason == self.pvmap['resolution']:
            value = self.sim.screen.resolution * 1e6 # um  
        else:
            value = self.getParam(reason)
        return value
    
    def write(self, reason, value):
        status = False
        if reason == self.pvmap['quadval']:
            try:
                self.sim.quad_value = value 
                self.setParam('sim_screen_sigma_x', self.sim.screen_sigma('x'))
                self.setParam('sim_screen_sigma_y', self.sim.screen_sigma('y'))
                self.updatePVs()
                status = True
            except Exception as ex:
                logger.info(ex)
                logger.debug(f'something wrong setting {reason} =  {value}')

        # store the values
        if status:
            logger.info(f'Setting {reason} =  {value}')
            self.setParam(reason, value)     
    
class BeamSimServer:
    def __init__(self,
                 *,
                 bunch_params = None,
                 screen_params,
                 beamline_info,
                 pvdb,
                ):
        pass
        
        
        
def start_server(config_name='LCLS2_OTR0H04', prefix=''):
    bunch_params, screen_params, pvmap, beamline_info = get_all_params(config_name=config_name)
    pvdb = make_pvdb(pvmap, screen_params)
    
    # Try a pv
    #pvtest = list(pvdb)[0]
    #logger.info(f'Checking that another server is not serving {pvtest}')    
    #val = epics.caget(pvtest)
    #if val is not None:
    #    raise ValueError(f'Another server is serving {pvtest}, aborting')
    
    beamsim = BeamSim(bunch_params=bunch_params,
              beamline_info=beamline_info,
             screen=Screen(**screen_params),
                 )
    logger.info(f'Initialized config: {config_name}')
    
    server = SimpleServer()
    server.createPV(prefix, pvdb)
    driver = BeamSimDriver(beamsim, pvmap)
    
    logger.info(f'Serving: {list(pvdb)}')

    # process CA transactions
    while True:
        server.process(0.1)      
        
        

if __name__ == '__main__':
    start_server()
  
    