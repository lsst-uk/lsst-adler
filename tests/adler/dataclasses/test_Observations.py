from numpy.testing import assert_almost_equal
import numpy as np

from adler.dataclasses.Observations import Observations
from adler.dataclasses.dataclass_utilities import get_data_table
from adler.utilities.tests_utilities import get_test_data_filepath


def test_populate_from_data_table():
    ssoid = 8268570668335894776
    test_db_path = get_test_data_filepath("testing_database.db")
    schema = ""
    filter_name = "r"
    date_range = [61000.0, 62000.0]

    test_query = f"""
                    SELECT
                        ssObject.ssObjectId, mag, magErr, band, midPointMjdTai, ra, dec, phaseAngle,
                        topocentricDist, heliocentricDist
                    FROM
                        ssObject
                        JOIN {schema}diaSource ON {schema}ssObject.ssObjectId   = {schema}diaSource.ssObjectId
                        JOIN {schema}ssSource  ON {schema}diaSource.diaSourceId = {schema}ssSource.diaSourceId
                    WHERE
                        ssObject.ssObjectId = {ssoid} AND band = '{filter_name}' AND midPointMjdTai BETWEEN {date_range[0]} AND {date_range[1]}
                    """

    data_table = get_data_table(test_query, sql_filename=test_db_path)
    test_observations = Observations.construct_from_data_table(ssoid, filter_name, data_table)

    expected_mag = np.array(
        [
            20.4470005,
            22.24799919,
            22.27599907,
            22.29100037,
            22.66200066,
            22.69799995,
            23.56800079,
            23.70199966,
            22.22900009,
            22.04899979,
            22.17700005,
            23.52199936,
        ]
    )
    expected_magerr = np.array(
        [0.011, 0.046, 0.06, 0.082, 0.083, 0.063, 0.17200001, 0.199, 0.264, 0.28299999, 0.31600001, 0.163]
    )
    expected_mjd = np.array(
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
    expected_ra = np.array(
        [
            301.5682143,
            315.2799036,
            315.6390597,
            316.7793075,
            317.8804823,
            318.2527231,
            327.0967385,
            327.0968451,
            170.4323598,
            170.4247983,
            170.4172297,
            62.6767601,
        ]
    )
    expected_dec = np.array(
        [
            -17.709034,
            -13.7156598,
            -13.6098197,
            -13.2712987,
            -12.9420707,
            -12.8300092,
            -10.035292,
            -10.0351984,
            -3.6550378,
            -3.6513877,
            -3.6478445,
            27.0535373,
        ]
    )
    expected_phaseangle = np.array(
        [
            33.56268311,
            32.00559616,
            31.97466469,
            31.86707687,
            31.74268913,
            31.69605255,
            29.76596832,
            29.76592445,
            126.7827301,
            126.78705597,
            126.79136658,
            18.63666534,
        ]
    )
    expected_heliodist = np.array(
        [
            1.35723567,
            1.66025329,
            1.66984999,
            1.70058703,
            1.73042095,
            1.74053574,
            1.97701716,
            1.97702122,
            0.87671232,
            0.87667543,
            0.87663853,
            2.21820784,
        ]
    )
    expected_topodist = np.array(
        [
            0.45957258,
            0.93472457,
            0.95206416,
            1.00859618,
            1.06491303,
            1.08432984,
            1.58487654,
            1.58488584,
            0.20776999,
            0.20780498,
            0.20783997,
            1.4200809,
        ]
    )
    expected_reduced_mag = np.array(
        [
            21.47195362,
            21.29370916,
            21.2692807,
            21.11941945,
            21.33467112,
            21.31877818,
            21.08797141,
            21.22195309,
            25.92680044,
            25.74652589,
            25.87425196,
            21.03042275,
        ]
    )

    assert_almost_equal(test_observations.mag, expected_mag)
    assert_almost_equal(test_observations.magErr, expected_magerr)
    assert_almost_equal(test_observations.midpointMjdTai, expected_mjd)
    assert_almost_equal(test_observations.ra, expected_ra)
    assert_almost_equal(test_observations.dec, expected_dec)
    assert_almost_equal(test_observations.phaseAngle, expected_phaseangle)
    assert_almost_equal(test_observations.heliocentricDist, expected_heliodist)
    assert_almost_equal(test_observations.topocentricDist, expected_topodist)
    assert_almost_equal(test_observations.reduced_mag, expected_reduced_mag)


def test_calculate_reduced_mag():
    test_object = Observations(
        mag=np.array([20]), topocentricDist=np.array([0.5]), heliocentricDist=np.array([1.5])
    )
    reduced_mag = test_object.calculate_reduced_mag(
        test_object.mag, test_object.topocentricDist, test_object.heliocentricDist
    )[0]

    expected_reduced_mag = 20.0 - 5 * np.log10(0.5 * 1.5)

    assert_almost_equal(reduced_mag, expected_reduced_mag)
