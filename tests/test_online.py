from pyemittance import PyEmittance
import numpy as np


def test_emittance_calc_online():
    meas = PyEmittance(config_name='LCLS2_OTR0H04', online=True)
    meas.config_dict['img_proc']['subtract_bg'] = False
    meas.quad_init = list(np.linspace(0,5,7))
    meas.config_dict['img_proc']['n_to_acquire']=1
    result = meas.measure_emittance()
    assert result is not None