from pyemittance import PyEmittance


def test_emittance_calc():
    meas = PyEmittance()
    result = meas.measure_emittance()
    assert result is not None