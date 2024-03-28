import numpy as np
from numpy.testing import assert_almost_equal

from adler.utilities.tests_utilities import get_test_data_filepath
from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid


ssoid = 8268570668335894776
test_db_path = get_test_data_filepath("testing_database.db")


def test_construct_from_SQL():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path)

    # testing just a few values here to ensure correct setup: these objects have their own unit tests
    assert test_planetoid.MPCORB.mpcH == 19.8799991607666
    assert test_planetoid.SSObject.discoverySubmissionDate == 60218.0
    assert_almost_equal(
        test_planetoid.observations_by_filter[1].mag,
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

    # did we pick up all the filters?
    assert len(test_planetoid.observations_by_filter) == 6
    assert len(test_planetoid.SSObject.H) == 6
    assert test_planetoid.filter_list == ["u", "g", "r", "i", "z", "y"]

    # checking the date range to ensure it's the default
    assert test_planetoid.date_range == [60000.0, 67300.0]


def test_construct_with_single_filter():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path, filter_list=["g"])

    # should only be one filter in here now
    assert len(test_planetoid.observations_by_filter) == 1
    assert len(test_planetoid.SSObject.H) == 1
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

    assert_almost_equal(test_planetoid.observations_by_filter[0].midpointMjdTai, expected_dates)


def test_observations_in_filter():
    test_planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, test_db_path)

    assert_almost_equal(
        test_planetoid.observations_in_filter("g").mag,
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
