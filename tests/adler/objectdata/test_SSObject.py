from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
import numpy as np

from adler.utilities.tests_utilities import get_test_data_filepath
from adler.objectdata.objectdata_utilities import get_data_table
from adler.objectdata.SSObject import SSObject


def test_construct_SSObject_from_data_table():
    ssoid = 8268570668335894776
    test_db_path = get_test_data_filepath("testing_database.db")

    filter_list = ["r", "g"]

    filter_dependent_columns = ""

    for filter_name in filter_list:
        filter_string = "{}_H, {}_G12, {}_HErr, {}_G12Err, {}_Ndata, ".format(
            filter_name, filter_name, filter_name, filter_name, filter_name
        )

        filter_dependent_columns += filter_string

    test_query = f"""
        SELECT
            discoverySubmissionDate, firstObservationDate, arc, numObs, 
            {filter_dependent_columns}
            maxExtendedness, minExtendedness, medianExtendedness
        FROM
            SSObject
        WHERE
            ssObjectId = {ssoid}
    """

    data_table = get_data_table(test_query, sql_filename=test_db_path)
    test_SSObject = SSObject.construct_from_data_table(ssoid, filter_list, data_table)

    test_values_r = [
        test_SSObject.filter_dependent_values[0].H,
        test_SSObject.filter_dependent_values[0].G12,
        test_SSObject.filter_dependent_values[0].Herr,
        test_SSObject.filter_dependent_values[0].G12err,
        test_SSObject.filter_dependent_values[0].nData,
    ]

    test_values_g = [
        test_SSObject.filter_dependent_values[1].H,
        test_SSObject.filter_dependent_values[1].G12,
        test_SSObject.filter_dependent_values[1].Herr,
        test_SSObject.filter_dependent_values[1].G12err,
        test_SSObject.filter_dependent_values[1].nData,
    ]

    expected_values_g = [20.292325973510742, 1.7233933210372925, 0.030210301280021667, 0.0404973067343235, 9]

    expected_values_r = [19.805892944335938, 1.52932608127594, 0.01974303089082241, 0.05071713775396347, 38]

    assert test_SSObject.ssObjectId == 8268570668335894776
    assert test_SSObject.filter_list == filter_list
    assert_almost_equal(test_SSObject.discoverySubmissionDate, 60218.0, decimal=6)
    assert_almost_equal(test_SSObject.firstObservationDate, 60220.01958, decimal=6)
    assert_almost_equal(test_SSObject.arc, 3342.05859375, decimal=6)
    assert test_SSObject.numObs == 94
    assert_equal(test_SSObject.maxExtendedness, 0.0)
    assert_equal(test_SSObject.minExtendedness, 0.0)
    assert_equal(test_SSObject.medianExtendedness, 0.0)

    assert_almost_equal(test_values_r, expected_values_r)
    assert_almost_equal(test_values_g, expected_values_g)
