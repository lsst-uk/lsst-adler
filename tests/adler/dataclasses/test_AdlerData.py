import os
import pytest
import numpy as np
import pandas as pd
import sqlite3

from numpy.testing import assert_array_equal

from adler.dataclasses.AdlerData import AdlerData
from adler.utilities.tests_utilities import get_test_data_filepath

# setting up the AdlerData object to be used for testing

test_object = AdlerData("8268570668335894776", ["u", "g", "r"])

u_model_1 = {
    "model_name": "model_1",
    "phaseAngle_min": 11.0,
    "phaseAngle_range": 12.0,
    "nobs": 13,
    "arc": 14.0,
    "H": 15.0,
    "H_err": 16.0,
    "phase_parameter_1": 17.0,
    "phase_parameter_1_err": 18.0,
}

u_model_2 = {
    "model_name": "model_2",
    "H": 25.0,
    "H_err": 26.0,
    "phase_parameter_1": 27.0,
    "phase_parameter_1_err": 28.0,
    "phase_parameter_2": 29.0,
    "phase_parameter_2_err": 30.0,
}

g_model_1 = {
    "model_name": "model_1",
    "phaseAngle_min": 31.0,
    "phaseAngle_range": 32.0,
    "nobs": 33,
    "arc": 34.0,
    "H": 35.0,
    "H_err": 36.0,
    "phase_parameter_1": 37.0,
    "phase_parameter_1_err": 38.0,
}

r_model_2 = {
    "model_name": "model_2",
    "phaseAngle_min": 41.0,
    "phaseAngle_range": 42.0,
    "nobs": 43,
    "arc": 44.0,
    "H": 45.0,
    "H_err": 46.0,
    "phase_parameter_1": 47.0,
    "phase_parameter_1_err": 48.0,
    "phase_parameter_2": 49.0,
    "phase_parameter_2_err": 50.0,
}

test_object.populate_phase_parameters("u", **u_model_1)
test_object.populate_phase_parameters("u", **u_model_2)
test_object.populate_phase_parameters("g", **g_model_1)
test_object.populate_phase_parameters("r", **r_model_2)


def test_populate_phase_parameters():
    # test to make sure the object is correctly populated
    assert test_object.filter_list == ["u", "g", "r"]

    assert test_object.filter_dependent_values[0].model_list == ["model_1", "model_2"]
    assert test_object.filter_dependent_values[1].model_list == ["model_1"]
    assert test_object.filter_dependent_values[2].model_list == ["model_2"]

    assert_array_equal([a.phaseAngle_min for a in test_object.filter_dependent_values], [11.0, 31.0, 41.0])
    assert_array_equal([a.phaseAngle_range for a in test_object.filter_dependent_values], [12.0, 32.0, 42.0])
    assert_array_equal([a.nobs for a in test_object.filter_dependent_values], [13, 33, 43])
    assert_array_equal([a.arc for a in test_object.filter_dependent_values], [14.0, 34.0, 44.0])

    assert test_object.filter_dependent_values[0].model_dependent_values[0].__dict__ == {
        "filter_name": "u",
        "model_name": "model_1",
        "H": 15.0,
        "H_err": 16.0,
        "phase_parameter_1": 17.0,
        "phase_parameter_1_err": 18.0,
        "phase_parameter_2": np.nan,
        "phase_parameter_2_err": np.nan,
    }

    assert test_object.filter_dependent_values[0].model_dependent_values[1].__dict__ == {
        "filter_name": "u",
        "model_name": "model_2",
        "H": 25.0,
        "H_err": 26.0,
        "phase_parameter_1": 27.0,
        "phase_parameter_1_err": 28.0,
        "phase_parameter_2": 29.0,
        "phase_parameter_2_err": 30.0,
    }

    assert test_object.filter_dependent_values[1].model_dependent_values[0].__dict__ == {
        "filter_name": "g",
        "model_name": "model_1",
        "H": 35.0,
        "H_err": 36.0,
        "phase_parameter_1": 37.0,
        "phase_parameter_1_err": 38.0,
        "phase_parameter_2": np.nan,
        "phase_parameter_2_err": np.nan,
    }

    assert test_object.filter_dependent_values[2].model_dependent_values[0].__dict__ == {
        "filter_name": "r",
        "model_name": "model_2",
        "H": 45.0,
        "H_err": 46.0,
        "phase_parameter_1": 47.0,
        "phase_parameter_1_err": 48.0,
        "phase_parameter_2": 49.0,
        "phase_parameter_2_err": 50.0,
    }

    # check to make sure model-dependent parameter is correctly updated (then return it to previous)
    test_object.populate_phase_parameters("u", model_name="model_1", H=99.0)
    assert test_object.filter_dependent_values[0].model_dependent_values[0].H == 99.0
    test_object.populate_phase_parameters("u", model_name="model_1", H=15.0)

    # check to make sure filter-dependent parameter is correctly updated (then return it to previous)
    test_object.populate_phase_parameters("u", nobs=99)
    assert test_object.filter_dependent_values[0].nobs == 99
    test_object.populate_phase_parameters("u", nobs=13)

    # testing to make sure the correct error messages trigger
    with pytest.raises(ValueError) as error_info_1:
        test_object.populate_phase_parameters("y")

    assert error_info_1.value.args[0] == "Filter y does not exist in AdlerData.filter_list."

    with pytest.raises(NameError) as error_info_2:
        test_object.populate_phase_parameters("u", H=4.0)

    assert error_info_2.value.args[0] == "No model name given. Cannot update model-specific phase parameters."


def test_get_phase_parameters_in_filter():
    assert test_object.get_phase_parameters_in_filter("u", "model_1").__dict__ == {
        "filter_name": "u",
        "phaseAngle_min": 11.0,
        "phaseAngle_range": 12.0,
        "nobs": 13,
        "arc": 14.0,
        "model_name": "model_1",
        "H": 15.0,
        "H_err": 16.0,
        "phase_parameter_1": 17.0,
        "phase_parameter_1_err": 18.0,
        "phase_parameter_2": np.nan,
        "phase_parameter_2_err": np.nan,
    }

    assert test_object.get_phase_parameters_in_filter("u", "model_2").__dict__ == {
        "filter_name": "u",
        "phaseAngle_min": 11.0,
        "phaseAngle_range": 12.0,
        "nobs": 13,
        "arc": 14.0,
        "model_name": "model_2",
        "H": 25.0,
        "H_err": 26.0,
        "phase_parameter_1": 27.0,
        "phase_parameter_1_err": 28.0,
        "phase_parameter_2": 29.0,
        "phase_parameter_2_err": 30.0,
    }

    assert test_object.get_phase_parameters_in_filter("g", "model_1").__dict__ == {
        "filter_name": "g",
        "phaseAngle_min": 31.0,
        "phaseAngle_range": 32.0,
        "nobs": 33,
        "arc": 34.0,
        "model_name": "model_1",
        "H": 35.0,
        "H_err": 36.0,
        "phase_parameter_1": 37.0,
        "phase_parameter_1_err": 38.0,
        "phase_parameter_2": np.nan,
        "phase_parameter_2_err": np.nan,
    }

    assert test_object.get_phase_parameters_in_filter("r", "model_2").__dict__ == {
        "filter_name": "r",
        "phaseAngle_min": 41.0,
        "phaseAngle_range": 42.0,
        "nobs": 43,
        "arc": 44.0,
        "model_name": "model_2",
        "H": 45.0,
        "H_err": 46.0,
        "phase_parameter_1": 47.0,
        "phase_parameter_1_err": 48.0,
        "phase_parameter_2": 49.0,
        "phase_parameter_2_err": 50.0,
    }

    # checking the error messages
    with pytest.raises(ValueError) as error_info_1:
        error_dict = test_object.get_phase_parameters_in_filter("f", model_name="model_2")

    assert error_info_1.value.args[0] == "Filter f does not exist in AdlerData.filter_list."

    with pytest.raises(ValueError) as error_info_2:
        error_dict_2 = test_object.get_phase_parameters_in_filter("r", model_name="model_1")

    assert error_info_2.value.args[0] == "Model model_1 does not exist for filter r in AdlerData.model_lists."


# here the capsys fixture captures any output to the terminal
def test_print_data(capsys):
    test_object.print_data()

    # get what was printed to the terminal
    captured = capsys.readouterr()

    expected = "Filter: u\nPhase angle minimum: 11.0\nPhase angle range: 12.0\nNumber of observations: 13\nArc: 14.0\nModel: model_1.\n\tH: 15.0\n\tH error: 16.0\n\tPhase parameter 1: 17.0\n\tPhase parameter 1 error: 18.0\n\tPhase parameter 2: nan\n\tPhase parameter 2 error: nan\nModel: model_2.\n\tH: 25.0\n\tH error: 26.0\n\tPhase parameter 1: 27.0\n\tPhase parameter 1 error: 28.0\n\tPhase parameter 2: 29.0\n\tPhase parameter 2 error: 30.0\n\n\nFilter: g\nPhase angle minimum: 31.0\nPhase angle range: 32.0\nNumber of observations: 33\nArc: 34.0\nModel: model_1.\n\tH: 35.0\n\tH error: 36.0\n\tPhase parameter 1: 37.0\n\tPhase parameter 1 error: 38.0\n\tPhase parameter 2: nan\n\tPhase parameter 2 error: nan\n\n\nFilter: r\nPhase angle minimum: 41.0\nPhase angle range: 42.0\nNumber of observations: 43\nArc: 44.0\nModel: model_2.\n\tH: 45.0\n\tH error: 46.0\n\tPhase parameter 1: 47.0\n\tPhase parameter 1 error: 48.0\n\tPhase parameter 2: 49.0\n\tPhase parameter 2 error: 50.0\n\n\n"

    assert captured.out == expected


def test_write_row_to_database(tmp_path):
    db_location = os.path.join(tmp_path, "test_AdlerData_database.db")
    test_object.write_row_to_database(db_location)

    con = sqlite3.connect(db_location)
    written_data = pd.read_sql_query("SELECT * from AdlerData", con)
    con.close()

    expected_data_filepath = get_test_data_filepath("test_SQL_database_table.csv")
    expected_data = pd.read_csv(expected_data_filepath)

    # we don't expect the timestamp column to be the same, obviously
    expected_data = expected_data.drop(columns="timestamp")
    written_data = written_data.drop(columns="timestamp")

    # note that because I'm using Pandas there's some small dtype and np.nan/None stuff to clear up
    # but this makes for a quick streamlined test anyway
    expected_data = expected_data.replace({np.nan: None})
    expected_data = expected_data.astype({"ssObjectId": str})
    pd.testing.assert_frame_equal(expected_data, written_data, check_dtype=False)


def test_read_row_from_database():
    # NOTE: the test database here has two rows, one with an earlier timestamp and different data
    # So this test also ensures that only the most recent data for the object is pulled.

    db_location = get_test_data_filepath("test_AdlerData_database.db")

    new_object = AdlerData("8268570668335894776", ["u", "g", "r"])
    new_object.populate_from_database(db_location)

    assert new_object.__dict__ == test_object.__dict__

    with pytest.raises(ValueError) as error_info_1:
        empty_data = AdlerData("pretend_object", ["u", "g", "r"])
        empty_data.populate_from_database(db_location)

    assert error_info_1.value.args[0] == "No data found in this database for the supplied ssObjectId."

    with pytest.raises(ValueError) as error_info_2:
        bad_filter = AdlerData("8268570668335894776", ["u", "g", "h"])
        bad_filter.populate_from_database(db_location)

    assert (
        error_info_2.value.args[0]
        == "Data does not exist for some of the requested filters in this database. Filters in database for this object: ['u', 'g', 'r']"
    )

    with pytest.raises(ValueError) as error_info_3:
        bad_filter = AdlerData("8268570668335894776", ["u", "g", "h"])
        bad_filter.populate_from_database("./dummy_location.db")

    assert error_info_3.value.args[0] == "Database cannot be found at given filepath."
