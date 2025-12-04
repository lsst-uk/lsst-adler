import pytest
import numpy as np
from numpy.testing import assert_almost_equal

from adler.utilities.tests_utilities import get_test_data_filepath
from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid


ssoid = "8268570668335894776"
test_db_path = get_test_data_filepath("testing_database.db")


def test_construct_from_SQL():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(
        ssoid, test_db_path, filter_list=["u", "g", "r", "i", "z", "y"]
    )

    # testing just a few values here to ensure correct setup: these objects have their own unit tests
    assert test_planetoid.MPCORB.mpcH == 19.8799991607666
    assert test_planetoid.SSObject.discoverySubmissionDate == 60218.0
    assert_almost_equal(
        test_planetoid.observations_by_filter[0].mag,
        [
            21.33099937,
            22.67099953,
            23.5359993,
            22.85000038,
            22.97599983,
            22.94499969,
            23.13599968,
            23.19400024,
            23.1609993,
        ],
    )

    # did we pick up all the filters? note we ask for ugrizy but u and y are unpopulated in DP0.3, so the code should eliminate them
    assert len(test_planetoid.observations_by_filter) == 4
    assert len(test_planetoid.SSObject.filter_dependent_values) == 4
    assert test_planetoid.filter_list == ["g", "r", "i", "z"]

    # checking the date range to ensure it's the default
    assert test_planetoid.date_range == [60000.0, 67300.0]


def test_construct_with_single_filter():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path, filter_list=["g"])

    # should only be one filter in here now
    assert len(test_planetoid.observations_by_filter) == 1
    assert len(test_planetoid.SSObject.filter_dependent_values) == 1
    assert test_planetoid.filter_list == ["g"]

    assert_almost_equal(
        test_planetoid.observations_by_filter[0].mag,
        [
            21.33099937,
            22.67099953,
            23.5359993,
            22.85000038,
            22.97599983,
            22.94499969,
            23.13599968,
            23.19400024,
            23.1609993,
        ],
    )


def test_construct_with_date_range():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(
        ssoid, test_db_path, filter_list=["r"], date_range=[61000.0, 62000.0]
    )

    expected_dates = np.array(
        [
            61294.15865,
            61322.07319,
            61323.00925,
            61326.03134,
            61329.00043,
            61330.01524,
            61355.02232,
            61355.02277,
            61253.97055,
            61253.96744,
            61253.96433,
            61052.13729,
        ]
    )

    assert_almost_equal(test_planetoid.observations_by_filter[0].midPointMjdTai, expected_dates)

    with pytest.raises(ValueError) as error_info_1:
        test_planetoid = AdlerPlanetoid.construct_from_SQL(
            ssoid, test_db_path, date_range=[61000.0, 62000.0, 63000.0]
        )

    assert error_info_1.value.args[0] == "date_range attribute must be of length 2."


def test_observations_in_filter():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path)

    # Python dataclasses create an __eq__ for you so object-to-object comparison just works, isn't that nice?
    assert test_planetoid.observations_in_filter("g") == test_planetoid.observations_by_filter[0]
    assert test_planetoid.observations_in_filter("r") == test_planetoid.observations_by_filter[1]
    assert test_planetoid.observations_in_filter("i") == test_planetoid.observations_by_filter[2]
    assert test_planetoid.observations_in_filter("z") == test_planetoid.observations_by_filter[3]

    with pytest.raises(ValueError) as error_info_1:
        test_planetoid.observations_in_filter("f")

    assert error_info_1.value.args[0] == "Filter f is not in AdlerPlanetoid.filter_list."


def test_SSObject_in_filter():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path)

    assert test_planetoid.SSObject_in_filter("g") == test_planetoid.SSObject.filter_dependent_values[0]
    assert test_planetoid.SSObject_in_filter("r") == test_planetoid.SSObject.filter_dependent_values[1]
    assert test_planetoid.SSObject_in_filter("i") == test_planetoid.SSObject.filter_dependent_values[2]
    assert test_planetoid.SSObject_in_filter("z") == test_planetoid.SSObject.filter_dependent_values[3]

    with pytest.raises(ValueError) as error_info_1:
        test_planetoid.SSObject_in_filter("f")

    assert error_info_1.value.args[0] == "Filter f is not in AdlerPlanetoid.filter_list."


def test_no_observations():
    with pytest.raises(Exception) as error_info:
        test_planetoid = AdlerPlanetoid.construct_from_SQL(826857066833589477, test_db_path)

    assert (
        error_info.value.args[0]
        == "No observations found for this object in the given filter(s). Check SSOID and try again."
    )


def test_failed_SQL_queries():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(
        ssoid, test_db_path, filter_list=["u", "g", "r", "i", "z", "y"]
    )

    with pytest.raises(Exception) as error_info_1:
        test_planetoid.populate_MPCORB("826857066833589477", sql_filename=test_db_path, schema="")

    assert error_info_1.value.args[0] == "No MPCORB data for this object could be found for this SSObjectId."

    with pytest.raises(Exception) as error_info_2:
        test_planetoid.populate_SSObject(
            "826857066833589477", filter_list=["u"], sql_filename=test_db_path, schema=""
        )

    assert (
        error_info_2.value.args[0] == "No SSObject data for this object could be found for this SSObjectId."
    )


def test_attach_previous_adlerdata():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path, filter_list=["g", "r"])

    # TODO: this setup is currently a bit dodgy. Because AdlerData write_row_to_database appends a new line to the db the most recent model for a given object may be nan in that row
    # as such this test depends on the order/number of times the adler commands have been run to make it

    # the test database can be recreated by running the Adler commands in the tests/data dir:
    # adler -s 8268570668335894776 -i testing_database.db -n test_AdlerData_database.db
    # adler -s 8268570668335894776 -i testing_database.db -n test_AdlerData_database.db -m HG
    db_location = get_test_data_filepath("test_AdlerData_database.db")
    print(db_location)

    test_planetoid.attach_previous_adler_data(db_location)

    test_output = test_planetoid.PreviousAdlerData.get_phase_parameters_in_filter("r", "HG12_Pen16")
    print(test_output.__dict__)

    expected_output = {
        "filter_name": "r",
        "phaseAngle_min": 2.553332567214966,
        "phaseAngle_range": 124.23803400993347,
        "nobs": 38,
        "arc": 3338.0655999999944,
        "model_name": "HG12_Pen16",
        "H": 19.92863542616601,
        "H_err": 0.018525355171274356,
        "phase_parameter_1": 1.0,
        "phase_parameter_1_err": 0.05300829494059732,
        "phase_parameter_2": np.nan,
        "phase_parameter_2_err": np.nan,
    }
    print(expected_output)

    # assert test_output.__dict__ == pytest.approx(expected_output) # TODO: why isn't this working now?

    for x in test_output.__dict__.keys():
        print(x)
        test_val = test_output.__dict__[x]
        expect_val = expected_output[x]
        if type(expect_val) == str:
            assert test_val == expect_val
        else:
            assert_almost_equal(test_val, expect_val)


mpc_ssoid = "2025 MS22"
mpc_test_db_path = get_test_data_filepath("mpc_obs_sbn_testing_database.sqlite")


def test_construct_from_mpc_obs_sbn():
    test_planetoid = AdlerPlanetoid.construct_from_mpc_obs_sbn(
        mpc_ssoid, mpc_test_db_path, filter_list=["u", "g", "r", "i", "z", "y"]
    )

    # testing just a few values here to ensure correct setup: these objects have their own unit tests
    assert test_planetoid.MPCORB.mpcH == 18.54
    assert test_planetoid.SSObject.numObs == 73
    assert_almost_equal(
        test_planetoid.observations_by_filter[0].mag,
        [24.366, 24.173, 24.529, 24.48, 24.205, 24.058],
    )

    # did we pick up all the filters? note we ask for ugrizy but this object was not observed in uzy, so the code should eliminate them
    assert len(test_planetoid.observations_by_filter) == 3
    assert len(test_planetoid.SSObject.filter_dependent_values) == 3
    assert test_planetoid.filter_list == ["g", "r", "i"]

    # checking the date range to ensure it's the default
    assert test_planetoid.date_range == [60000.0, 67300.0]


# TODO these tests
def test_construct_from_mpc_with_single_filter():
    test_planetoid = AdlerPlanetoid.construct_from_mpc_obs_sbn(mpc_ssoid, mpc_test_db_path, filter_list=["g"])

    # should only be one filter in here now
    assert len(test_planetoid.observations_by_filter) == 1
    assert len(test_planetoid.SSObject.filter_dependent_values) == 1
    assert test_planetoid.filter_list == ["g"]

    assert_almost_equal(
        test_planetoid.observations_by_filter[0].mag,
        [24.366, 24.173, 24.529, 24.48, 24.205, 24.058],
    )


def test_construct_from_mpc_with_date_range():
    test_planetoid = AdlerPlanetoid.construct_from_mpc_obs_sbn(
        mpc_ssoid, mpc_test_db_path, filter_list=["g"], date_range=[60795.0, 60798.0]
    )

    expected_dates = np.array(
        [60797.11380400463, 60797.121616296296, 60797.12747140046, 60797.12846829861, 60797.130895902774]
    )

    assert_almost_equal(test_planetoid.observations_by_filter[0].midPointMjdTai, expected_dates)

    with pytest.raises(ValueError) as error_info_1:
        test_planetoid = AdlerPlanetoid.construct_from_mpc_obs_sbn(
            mpc_ssoid, mpc_test_db_path, date_range=[61000.0, 62000.0, 63000.0]
        )

    assert error_info_1.value.args[0] == "date_range argument must be of length 2."


def test_mpc_no_observations():
    with pytest.raises(Exception) as error_info:
        test_planetoid = AdlerPlanetoid.construct_from_mpc_obs_sbn("2025 FakeId", mpc_test_db_path)

    assert (
        error_info.value.args[0]
        == "No observations found for this object in the given filter(s). Check SSOID and try again."
    )


def test_mpc_failed_SQL_queries():
    test_planetoid = AdlerPlanetoid.construct_from_mpc_obs_sbn(
        mpc_ssoid, mpc_test_db_path, filter_list=["u", "g", "r", "i", "z", "y"]
    )

    with pytest.raises(Exception) as error_info_1:
        test_planetoid.populate_MPCORB_from_mpc_obs_sbn("2025 FakeId", sql_filename=mpc_test_db_path)

    assert (
        error_info_1.value.args[0] == "No mpc_orbits data for this object could be found for this SSObjectId."
    )

    with pytest.raises(Exception) as error_info_2:
        test_planetoid.populate_SSObject_from_mpc_obs_sbn(
            "2025 FakeId", filter_list=["u"], sql_filename=mpc_test_db_path
        )

    assert (
        error_info_2.value.args[0] == "No SSObject data for this object could be found for this SSObjectId."
    )
