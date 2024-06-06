from numpy.testing import assert_array_equal
import pytest


def test_PhaseCurve_init():
    import numpy as np
    import astropy.units as u
    from adler.science.PhaseCurve import PhaseCurve

    # test creating a model
    H = 18.9
    G = 0.12
    pc = PhaseCurve(H=H * u.mag, phase_parameter_1=G, model_name="HG")

    assert pc.H.value == 18.9
    assert pc.H.unit == u.mag
    assert pc.phase_parameter_1 == 0.12
    assert pc.model_name == "HG"


def test_PhaseCurve_ReducedMag():
    import numpy as np
    import astropy.units as u
    from adler.science.PhaseCurve import PhaseCurve

    # define the phase angles
    alpha = np.array([0, 10]) * u.deg

    # linear phase curve model
    pc_lin = PhaseCurve(model_name="LinearPhaseFunc", H=18 * u.mag, phase_parameter_1=0.1 * (u.mag / u.deg))

    # find the reduced mag
    red_mag = pc_lin.ReducedMag(alpha)

    assert red_mag.unit == u.mag
    assert_array_equal(red_mag.value, np.array([18.0, 19.0]))


def test_PhaseCurve_FitModel():
    import numpy as np
    import astropy.units as u
    from adler.science.PhaseCurve import PhaseCurve

    # define the observations
    alpha = np.array([0, 10]) * u.deg
    red_mag = np.array([18.0, 19.0]) * u.mag

    # empty linear phase curve model
    pc_lin = PhaseCurve(model_name="LinearPhaseFunc")

    # fit the model to the data
    pc_fit = pc_lin.FitModel(alpha, red_mag)

    assert pc_fit.H.unit == u.mag
    assert pc_fit.H.value == 18.0
    assert pc_fit.S.unit == u.mag / u.deg
    assert pc_fit.S.value == 0.1


### TODO - add more phase curve tests, test bounds, fixed values and fit uncertainties
