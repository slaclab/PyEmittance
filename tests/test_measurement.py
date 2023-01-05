from pyemittance import PyEmittance


def test_emittance_calc():
    meas = PyEmittance(config_name='LCLS2_OTR0H04')
    result = meas.measure_emittance()
    assert result is not None