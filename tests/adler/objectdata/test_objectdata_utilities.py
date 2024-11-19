import pytest
import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal
from numpy.testing import assert_equal

from adler.objectdata.objectdata_utilities import get_data_table
from adler.objectdata.objectdata_utilities import get_from_table
from adler.objectdata.objectdata_utilities import check_value_populated
from adler.utilities.tests_utilities import get_test_data_filepath


def test_get_data_table():
    ssoid = 8268570668335894776
    test_db_path = get_test_data_filepath("testing_database.db")
    filter_name = "r"
    date_range = [61000.0, 62000.0]

    test_query = f"""
                    SELECT
                        ssObject.ssObjectId, mag, magErr, band, midPointMjdTai, ra, dec, phaseAngle,
                        topocentricDist, heliocentricDist
                    FROM
                        ssObject
                        JOIN diaSource ON ssObject.ssObjectId   = diaSource.ssObjectId
                        JOIN ssSource  ON diaSource.diaSourceId = ssSource.diaSourceId
                    WHERE
                        ssObject.ssObjectId = {ssoid} AND band = '{filter_name}' AND midPointMjdTai BETWEEN {date_range[0]} AND {date_range[1]}
                    """

    data_table = get_data_table(test_query, sql_filename=test_db_path)

    expected_table = pd.read_csv(get_test_data_filepath("test_dataclass_utilities_table.csv"))
    assert_frame_equal(data_table, expected_table)


def test_get_from_table():
    test_table = pd.DataFrame(
        {"string_col": "a test string", "int_col": 4, "float_col": 4.5, "array_col": [5, 6]}
    )

    assert get_from_table(test_table, "string_col", str) == "a test string"
    assert get_from_table(test_table, "int_col", int) == 4
    assert get_from_table(test_table, "float_col", float) == 4.5
    assert_equal(get_from_table(test_table, "array_col", np.ndarray), [5, 6])

    with pytest.raises(ValueError) as error_info_1:
        get_from_table(test_table, "string_col", int)

    assert error_info_1.value.args[0] == "Could not cast column name to type."

    with pytest.raises(TypeError) as error_info_2:
        get_from_table(test_table, "string_col", "fake")

    assert (
        error_info_2.value.args[0]
        == "Type for argument data_type not recognised for column string_col in table default: must be str, float, int or np.ndarray."
    )


def test_check_value_populated():
    populated_value = check_value_populated(3, int, "column", "table")
    assert populated_value == 3

    array_length_zero = check_value_populated(np.array([]), np.ndarray, "column", "table")
    number_is_nan = check_value_populated(np.nan, float, "column", "table")
    str_is_empty = check_value_populated("", str, "column", "table")

    assert np.isnan(array_length_zero)
    assert np.isnan(number_is_nan)
    assert np.isnan(str_is_empty)
