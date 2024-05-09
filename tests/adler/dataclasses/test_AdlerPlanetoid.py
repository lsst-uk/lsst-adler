import pytest
import numpy as np
from numpy.testing import assert_almost_equal

from adler.utilities.tests_utilities import get_test_data_filepath
from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid


ssoid = 8268570668335894776
test_db_path = get_test_data_filepath("testing_database.db")


def test_construct_from_SQL():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path, filter_list=["u", "g", "r", "i", "z", "y"])

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

    assert error_info.value.args[0] == "No observations found for this object in the given filter(s). Check SSOID and try again."


def test_for_warnings(capsys):

    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path, filter_list=["u", "g"])
    captured = capsys.readouterr()

    expected = ("WARNING: No observations found in u filter for this object. Skipping this filter.\n"
    + "WARNING: n unpopulated in MPCORB table for this object. Storing NaN instead.\n"
    + "WARNING: uncertaintyParameter unpopulated in MPCORB table for this object. Storing NaN instead.\n")

    assert captured.out == expected
    