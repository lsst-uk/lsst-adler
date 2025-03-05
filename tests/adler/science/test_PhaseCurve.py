import numpy as np
import astropy.units as u
from numpy.testing import assert_array_equal, assert_array_less, assert_almost_equal

from adler.science.PhaseCurve import PhaseCurve, ADLER_SBPY_DICT
from adler.utilities.tests_utilities import get_test_data_filepath
from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid


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
    assert pc_fit.fixed["G"] is True
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

def test_PhaseCurve_ReducedMagBounds():
    """Test calculating the reduced magnitude from a PhaseCurve object."""

    # define the phase angles
    alpha = np.radians(np.array([0, 10]))

    # Test the HG12 model
    pc = PhaseCurve(model_name="HG12", H_err=0.1, phase_parameter_1_err=0.1)

    # find the reduced mag and the bounds
    red_mag = pc.ReducedMag(alpha)
    pc_bounds = pc.ReducedMagBounds(alpha)

    assert_almost_equal(pc_bounds["mag_min"][0], pc.H - pc.H_err)
    assert_almost_equal(pc_bounds["mag_max"][0], pc.H + pc.H_err)
    assert_almost_equal(pc_bounds["mag_mean"][0], pc.H)
    assert len(pc_bounds["PhaseCurves"]) == 4
    assert_array_less(pc_bounds["mag_min"], red_mag)
    assert_array_less(red_mag, pc_bounds["mag_max"])

    # also test the HG1G2 model with all uncertainties
    pc = PhaseCurve(model_name="HG1G2", H_err=0.1, phase_parameter_1_err=0.1, phase_parameter_2_err=0.1)
    # find the reduced mag and the bounds
    red_mag = pc.ReducedMag(alpha)
    pc_bounds = pc.ReducedMagBounds(alpha)

    assert_almost_equal(pc_bounds["mag_min"][0], pc.H - pc.H_err)
    assert_almost_equal(pc_bounds["mag_max"][0], pc.H + pc.H_err)
    assert_almost_equal(pc_bounds["mag_mean"][0], pc.H)
    assert len(pc_bounds["PhaseCurves"]) == 8
    assert_array_less(pc_bounds["mag_min"], red_mag)
    assert_array_less(red_mag, pc_bounds["mag_max"])

    # also test the HG1G2 model with only 2/3 uncertainties
    pc = PhaseCurve(model_name="HG1G2", H_err=0.1, phase_parameter_1_err=0.1)
    # find the reduced mag and the bounds
    red_mag = pc.ReducedMag(alpha)
    pc_bounds = pc.ReducedMagBounds(alpha)

    assert_almost_equal(pc_bounds["mag_min"][0], pc.H - pc.H_err)
    assert_almost_equal(pc_bounds["mag_max"][0], pc.H + pc.H_err)
    assert_almost_equal(pc_bounds["mag_mean"][0], pc.H)
    assert len(pc_bounds["PhaseCurves"]) == 4
    assert_array_less(pc_bounds["mag_min"], red_mag)
    assert_array_less(red_mag, pc_bounds["mag_max"])

def test_PhaseCurve_FitModel_resample():

    np.random.seed(0)  # set the seed to ensure reproducibility
    resample = 100  # number of resamples

    # these are the exact values that should be recovered with 100 resampled fits, for the random seed of 0
    resample_compare_vals = {
        "H": 16.29648544,
        "phase_parameter_1": 0.6225681900053228,
        "H_err": 0.00774547,
        "phase_parameter_1_err": 0.051412070779357485,
    }

    # load a test object
    ssoid = "6098332225018"  # good MBA test object
    filt = "r"
    test_db_path = get_test_data_filepath("testing_database.db")
    planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path, filter_list=[filt])
    sso = planetoid.SSObject_in_filter(filt)
    obs = planetoid.observations_in_filter(filt)

    # set up the initial phasecurve model
    H = sso.H
    G12 = sso.G12
    pc = PhaseCurve(H=H * u.mag, phase_parameter_1=G12, model_name="HG12_Pen16")

    # get the observations
    alpha = np.array(getattr(obs, "phaseAngle")) * u.deg
    mag_err = np.array(getattr(obs, "magErr")) * u.mag
    red_mag = np.array(getattr(obs, "reduced_mag")) * u.mag

    # do a single fit
    pc_fit = pc.FitModel(alpha, red_mag, mag_err)
    pc_fit = PhaseCurve().InitModelSbpy(pc_fit)
    print(pc_fit.__dict__)

    # use the resample function within FitModel
    pc_fit_resamp = pc.FitModel(alpha, red_mag, mag_err, resample=resample)
    pc_fit_resamp = PhaseCurve().InitModelSbpy(pc_fit_resamp)
    print(pc_fit_resamp.__dict__)

    # check that the fits are different
    tolerance = 0.001
    for x in ["H", "H_err", "phase_parameter_1", "phase_parameter_1_err"]:
        x1 = getattr(pc_fit_resamp, x)
        x2 = getattr(pc_fit, x)

        if hasattr(x1, "unit") and (x1.unit is not None):
            x1 = x1.value
        if hasattr(x2, "unit") and (x2.unit is not None):
            x2 = x2.value

        print(x, np.abs(x1 - x2))
        assert np.abs(x1 - x2) > tolerance

        # the resampled fit should have larger uncertainties
        if "err" in x:
            assert x1 > x2

    # check the exact values of each fit
    for x in ["H", "H_err", "phase_parameter_1", "phase_parameter_1_err"]:
        x1 = getattr(pc_fit_resamp, x)
        x2 = resample_compare_vals[x]

        if hasattr(x1, "unit") and (x1.unit is not None):
            x1 = x1.value

        print(x, x1, x2)
        assert_almost_equal(x1, x2)


def test_set_models():

    for model in ["HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc"]:

        # create a model
        pc = PhaseCurve(model_name=model)

        # check which phase parameters the model should have
        phase_params = ADLER_SBPY_DICT[pc.model_name]
        for p in phase_params:
            assert getattr(pc.model_function, phase_params[p])


def test_ReturnInit():

    # define a PhaseCurve object
    pc = PhaseCurve(H=18.0, phase_parameter_1=0.15, model_name="HG")

    # convert the PhaseCurve object to dict
    pc_dict = pc.ReturnModelDict()
    print(pc_dict)

    assert pc_dict["H"] == 18.0
    assert pc_dict["phase_parameter_1"] == 0.15
    assert pc_dict["model_name"] == "HG"

    # create PhaseCurve from dict
    print(type(pc_dict))
    pc2 = PhaseCurve().InitModelDict(pc_dict)

    assert pc2.H == 18.0
    assert pc2.phase_parameter_1 == 0.15
    assert pc2.model_name == "HG"

    # create sbpy model from PhaseCurve
    pc_sbpy = pc.model_function
    print(pc_sbpy)

    assert pc_sbpy.H == 18.0
    assert pc_sbpy.G == 0.15

    # create PhaseCurve from sbpy model
    pc3 = PhaseCurve().InitModelSbpy(pc_sbpy)

    assert pc3.H == 18.0
    assert pc3.phase_parameter_1 == 0.15
    assert pc3.model_name == "HG"


# TODO: test absmag and units
# TODO: test FitModel with weighted data
