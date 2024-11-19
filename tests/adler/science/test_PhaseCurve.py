from adler.science.PhaseCurve import PhaseCurve
from numpy.testing import assert_array_equal
import pytest
import numpy as np
import astropy.units as u


def test_PhaseCurve_init():
    """Test intialising the model."""

    # test creating a model
    H = 18.9
    G = 0.12
    pc = PhaseCurve(H=H * u.mag, phase_parameter_1=G, model_name="HG")

    assert pc.H.value == 18.9
    assert pc.H.unit == u.mag
    assert pc.phase_parameter_1 == 0.12
    assert pc.model_name == "HG"


def test_PhaseCurve_ReducedMag():
    """Test calculating the reduced magnitude from a PhaseCurve object."""

    # define the phase angles
    alpha = np.array([0, 10]) * u.deg

    # linear phase curve model
    pc_lin = PhaseCurve(model_name="LinearPhaseFunc", H=18 * u.mag, phase_parameter_1=0.1 * (u.mag / u.deg))

    # find the reduced mag
    red_mag = pc_lin.ReducedMag(alpha)

    assert red_mag.unit == u.mag
    assert_array_equal(red_mag.value, np.array([18.0, 19.0]))


def test_PhaseCurve_FitModel_linear():
    """Test fitting a linear phase function to two data points."""

    # define the observations
    alpha = np.array([0, 10]) * u.deg
    red_mag = np.array([18.0, 19.0]) * u.mag

    # empty linear phase curve model
    pc_lin = PhaseCurve(model_name="LinearPhaseFunc")

    # fit the model to the data - this returns an sbpy object
    pc_fit = pc_lin.FitModel(alpha, red_mag)

    assert pc_fit.H.unit == u.mag
    assert pc_fit.H.value == 18.0
    assert pc_fit.S.unit == u.mag / u.deg
    assert pc_fit.S.value == 0.1
    assert (
        hasattr(pc_fit, "H_err") == False
    )  # with only two data points the covariance matrix is None and no uncertainties are stored


def test_PhaseCurve_FitModel_HG():
    """Test fitting a HG model to generated data."""

    # generate some model data
    pc1 = PhaseCurve(H=18.0 * u.mag, phase_parameter_1=0.15, model_name="HG")
    alpha = np.linspace(0, 30) * u.deg
    red_mag = pc1.ReducedMag(alpha)

    # fit the same phase curve model to the data
    pc_fit = pc1.FitModel(alpha, red_mag)
    # convert from sbpy to adler PhaseCurve object
    pc2 = pc1.InitModelSbpy(pc_fit)

    # the new fitted model should have the same parameters as the input model
    assert pc2.H == pc1.H
    assert pc2.phase_parameter_1 == pc1.phase_parameter_1
    assert pc2.phase_parameter_2 is None
    assert pc1.phase_parameter_1_err is None  # the first model had no uncertainties
    assert pc2.phase_parameter_1_err is not None  # the fitted model has some uncertainties


def test_PhaseCurve_FitModel_HG_no_units():
    """Test fitting a HG model to generated data, but without units.
    If units are not provided, the phase angles must be in radians!"""

    # generate some model data
    pc1 = PhaseCurve(H=18.0, phase_parameter_1=0.15, model_name="HG")
    alpha = np.radians(np.linspace(0, 30))
    red_mag = pc1.ReducedMag(alpha)

    # fit the same phase curve model to the data
    pc_fit = pc1.FitModel(alpha, red_mag)
    # convert from sbpy to adler PhaseCurve object
    pc2 = pc1.InitModelSbpy(pc_fit)

    # the new fitted model should have the same parameters as the input model
    assert pc2.H == pc1.H
    assert pc2.phase_parameter_1 == pc1.phase_parameter_1
    assert pc2.phase_parameter_2 is None
    assert pc1.phase_parameter_1_err is None  # the first model had no uncertainties
    assert pc2.phase_parameter_1_err is not None  # the fitted model has some uncertainties


def test_PhaseCurve_FitModel_HG_fixed():
    """Test fitting a just H whilst keeping G fixed."""

    # generate some model data
    pc1 = PhaseCurve(H=18.0 * u.mag, phase_parameter_1=0.15, model_name="HG")
    alpha = np.linspace(0, 30) * u.deg
    red_mag = pc1.ReducedMag(alpha)

    # fix phase_parameter_1
    pc1.FixParam("phase_parameter_1")

    # fit the same phase curve model to the data
    pc_fit = pc1.FitModel(alpha, red_mag)
    # convert from sbpy to adler PhaseCurve object
    pc2 = pc1.InitModelSbpy(pc_fit)

    # the new fitted model should have the same parameters as the input model, but G is fixed
    assert pc_fit._fixed["G"] is True
    assert pc2.H == pc1.H
    assert pc2.phase_parameter_1 == pc1.phase_parameter_1
    assert pc2.phase_parameter_1_err is None  # the fitted model has no uncertainties when param is fixed


def test_PhaseCurve_FitModel_HG_bounds():
    """Test fitting a just H whilst keeping G fixed."""

    # generate some model data
    pc1 = PhaseCurve(H=18.0 * u.mag, phase_parameter_1=0.15, model_name="HG")
    alpha = np.linspace(0, 30) * u.deg
    red_mag = pc1.ReducedMag(alpha)

    # set bounds on phase parameter
    pc1.SetModelBounds("phase_parameter_1", (0.0, 0.1))

    # fit the same phase curve model to the data
    pc_fit = pc1.FitModel(alpha, red_mag)
    # convert from sbpy to adler PhaseCurve object
    pc2 = pc1.InitModelSbpy(pc_fit)

    # the new fitted model should have the same parameters as the input model, but G is fixed
    assert pc_fit.G.bounds == (0.0, 0.1)
    assert pc2.phase_parameter_1 == 0.1
    assert pc2.phase_parameter_1_err is not None


# TODO: test absmag
