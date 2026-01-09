import os
import pytest
import numpy as np
import pandas as pd
import sqlite3

from numpy.testing import assert_array_equal, assert_almost_equal

from adler.objectdata.AdlerData import AdlerData
from adler.utilities.tests_utilities import get_test_data_filepath

# setting up the AdlerData object to be used for testing

test_object = AdlerData("8268570668335894776", ["g", "r", "i"])
model_1 = "HG12_Pen16"
model_2 = "HG"

# test data was generated in the tests/data dir using:
# python adler_run_test_data.py
# the output is copied here
g_model_1 = {
    "filter_name": "g",
    "phaseAngle_min": 27.426668167114258,
    "phaseAngle_range": 36.17388343811035,
    "nobs": 9,
    "arc": 3306.0079999999944,
    "model_name": "HG12_Pen16",
    "H": 20.719465178426468,
    "H_err": 0.018444134023977106,
    "phase_parameter_1": 0.62,
    "phase_parameter_1_err": np.nan,
    "phase_parameter_2": np.nan,
    "phase_parameter_2_err": np.nan,
}
r_model_1 = {
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
i_model_1 = {
    "filter_name": "i",
    "phaseAngle_min": 10.0595064163208,
    "phaseAngle_range": 101.74150371551514,
    "nobs": 32,
    "arc": 3342.0585599999977,
    "model_name": "HG12_Pen16",
    "H": 19.653155995892117,
    "H_err": 0.01415178691419005,
    "phase_parameter_1": 1.0,
    "phase_parameter_1_err": 0.01516569963597136,
    "phase_parameter_2": np.nan,
    "phase_parameter_2_err": np.nan,
}
g_model_2 = {
    "filter_name": "g",
    "phaseAngle_min": 27.426668167114258,
    "phaseAngle_range": 36.17388343811035,
    "nobs": 9,
    "arc": 3306.0079999999944,
    "model_name": "HG",
    "H": 20.72954332871528,
    "H_err": 0.018444134023977106,
    "phase_parameter_1": 0.15,
    "phase_parameter_1_err": np.nan,
    "phase_parameter_2": np.nan,
    "phase_parameter_2_err": np.nan,
}
r_model_2 = {
    "filter_name": "r",
    "phaseAngle_min": 2.553332567214966,
    "phaseAngle_range": 124.23803400993347,
    "nobs": 38,
    "arc": 3338.0655999999944,
    "model_name": "HG",
    "H": 18.879873513575124,
    "H_err": 0.007456882346021157,
    "phase_parameter_1": -0.253,
    "phase_parameter_1_err": 1.6701987135262256e-05,
    "phase_parameter_2": np.nan,
    "phase_parameter_2_err": np.nan,
}
i_model_2 = {
    "filter_name": "i",
    "phaseAngle_min": 10.0595064163208,
    "phaseAngle_range": 101.74150371551514,
    "nobs": 32,
    "arc": 3342.0585599999977,
    "model_name": "HG",
    "H": 18.628894992991583,
    "H_err": 0.01716803531553185,
    "phase_parameter_1": -0.253,
    "phase_parameter_1_err": 0.00014324314956736168,
    "phase_parameter_2": np.nan,
    "phase_parameter_2_err": np.nan,
}


_g_model_1 = g_model_1.copy()
del _g_model_1["filter_name"]
_g_model_2 = g_model_2.copy()
del _g_model_2["filter_name"]
_r_model_1 = r_model_1.copy()
del _r_model_1["filter_name"]
_i_model_2 = i_model_2.copy()
del _i_model_2["filter_name"]

test_object.populate_phase_parameters("g", **_g_model_1)
test_object.populate_phase_parameters("g", **_g_model_2)
test_object.populate_phase_parameters("r", **_r_model_1)
test_object.populate_phase_parameters("i", **_i_model_2)


def test_populate_phase_parameters():
    # test to make sure the object is correctly populated
    assert test_object.filter_list == ["g", "r", "i"]

    assert test_object.filter_dependent_values[0].model_list == [model_1, model_2]
    assert test_object.filter_dependent_values[1].model_list == [model_1]
    assert test_object.filter_dependent_values[2].model_list == [model_2]

    assert_almost_equal(
        [a.phaseAngle_min for a in test_object.filter_dependent_values],
        [g_model_1["phaseAngle_min"], r_model_1["phaseAngle_min"], i_model_1["phaseAngle_min"]],
    )
    assert_array_equal(
        [a.phaseAngle_range for a in test_object.filter_dependent_values],
        [g_model_1["phaseAngle_range"], r_model_1["phaseAngle_range"], i_model_1["phaseAngle_range"]],
    )
    assert_array_equal(
        [a.nobs for a in test_object.filter_dependent_values],
        [g_model_1["nobs"], r_model_1["nobs"], i_model_1["nobs"]],
    )
    assert_array_equal(
        [a.arc for a in test_object.filter_dependent_values],
        [g_model_1["arc"], r_model_1["arc"], i_model_1["arc"]],
    )

    test_dict_list = [
        test_object.get_phase_parameters_in_filter("g", model_1).__dict__,
        test_object.get_phase_parameters_in_filter("g", model_2).__dict__,
        test_object.get_phase_parameters_in_filter("r", model_1).__dict__,
        test_object.get_phase_parameters_in_filter("i", model_2).__dict__,
    ]
    expect_dict_list = [g_model_1, g_model_2, r_model_1, i_model_2]

    for test_dict, expect_dict in zip(test_dict_list, expect_dict_list):
        for x in test_dict.keys():
            # print(x)
            test_val = test_dict[x]
            expect_val = expect_dict[x]
            if type(expect_val) == str:
                assert test_val == expect_val
            else:
                assert_almost_equal(test_val, expect_val)

    # check to make sure model-dependent parameter is correctly updated (then return it to previous)
    test_object.populate_phase_parameters("g", model_name=model_1, H=99.0)
    assert test_object.get_phase_parameters_in_filter("g", model_1).H == 99.0
    test_object.populate_phase_parameters("g", model_name=model_1, H=g_model_1["H"])

    # check to make sure filter-dependent parameter is correctly updated (then return it to previous)

    test_object.populate_phase_parameters("r", model_name=model_1, nobs=99)
    assert test_object.get_phase_parameters_in_filter("r", model_1).nobs == 99
    test_object.populate_phase_parameters("r", model_name=model_1, nobs=r_model_1["nobs"])

    # testing to make sure the correct error messages trigger
    test_filt = "y"
    with pytest.raises(ValueError) as error_info_1:
        test_object.populate_phase_parameters(test_filt)

    assert error_info_1.value.args[0] == "Filter {} does not exist in AdlerData.filter_list.".format(
        test_filt
    )

    with pytest.raises(NameError) as error_info_2:
        test_object.populate_phase_parameters("r", H=4.0)

    assert error_info_2.value.args[0] == "No model name given. Cannot update model-specific phase parameters."


def test_get_phase_parameters_in_filter():

    test_dict_list = [
        test_object.get_phase_parameters_in_filter("g", model_1).__dict__,
        test_object.get_phase_parameters_in_filter("g", model_2).__dict__,
        test_object.get_phase_parameters_in_filter("r", model_1).__dict__,
        test_object.get_phase_parameters_in_filter("i", model_2).__dict__,
    ]
    expect_dict_list = [g_model_1, g_model_2, r_model_1, i_model_2]

    for test_dict, expect_dict in zip(test_dict_list, expect_dict_list):
        for x in test_dict.keys():
            # print(x)
            test_val = test_dict[x]
            expect_val = expect_dict[x]
            if type(expect_val) == str:
                assert test_val == expect_val
            else:
                assert_almost_equal(test_val, expect_val)

    # checking the error messages
    test_filt = "f"
    with pytest.raises(ValueError) as error_info_1:
        error_dict = test_object.get_phase_parameters_in_filter(test_filt, model_name=model_2)

    assert error_info_1.value.args[0] == "Filter {} does not exist in AdlerData.filter_list.".format(
        test_filt
    )

    test_filt = "r"
    with pytest.raises(ValueError) as error_info_2:
        error_dict_2 = test_object.get_phase_parameters_in_filter(test_filt, model_name=model_2)

    print(error_info_2.value.args[0])
    print("Model {} does not exist for filter {} in AdlerData.model_lists.".format(model_2, test_filt))
    assert error_info_2.value.args[
        0
    ] == "Model {} does not exist for filter {} in AdlerData.model_lists.".format(model_2, test_filt)


# here the capsys fixture captures any output to the terminal
def test_print_data(capsys):
    test_object.print_data()

    # get what was printed to the terminal
    captured = capsys.readouterr()

    expected = "Filter: g\nPhase angle minimum: 27.426668167114258\nPhase angle range: 36.17388343811035\nNumber of observations: 9\nArc: 3306.0079999999944\nModel: HG12_Pen16.\n\tH: 20.719465178426468\n\tH error: 0.018444134023977106\n\tPhase parameter 1: 0.62\n\tPhase parameter 1 error: nan\n\tPhase parameter 2: nan\n\tPhase parameter 2 error: nan\nModel: HG.\n\tH: 20.72954332871528\n\tH error: 0.018444134023977106\n\tPhase parameter 1: 0.15\n\tPhase parameter 1 error: nan\n\tPhase parameter 2: nan\n\tPhase parameter 2 error: nan\n\n\nFilter: r\nPhase angle minimum: 2.553332567214966\nPhase angle range: 124.23803400993347\nNumber of observations: 38\nArc: 3338.0655999999944\nModel: HG12_Pen16.\n\tH: 19.92863542616601\n\tH error: 0.018525355171274356\n\tPhase parameter 1: 1.0\n\tPhase parameter 1 error: 0.05300829494059732\n\tPhase parameter 2: nan\n\tPhase parameter 2 error: nan\n\n\nFilter: i\nPhase angle minimum: 10.0595064163208\nPhase angle range: 101.74150371551514\nNumber of observations: 32\nArc: 3342.0585599999977\nModel: HG.\n\tH: 18.628894992991583\n\tH error: 0.01716803531553185\n\tPhase parameter 1: -0.253\n\tPhase parameter 1 error: 0.00014324314956736168\n\tPhase parameter 2: nan\n\tPhase parameter 2 error: nan\n\n\n"

    # print(captured)
    # print(expected)
    assert captured.out == expected


def test_write_row_to_database(tmp_path):
    db_location = os.path.join(tmp_path, "test_AdlerData_database.db")
    test_object.write_row_to_database(db_location)

    con = sqlite3.connect(db_location)
    written_data = pd.read_sql_query("SELECT * from AdlerData", con)
    con.close()

    expected_data_filepath = get_test_data_filepath("test_SQL_database_table.csv")
    expected_data = pd.read_csv(expected_data_filepath, index_col=0)

    # we don't expect the timestamp column to be the same, obviously
    drop_cols = [x for x in written_data if (x == "timestamp") or ("modelFitMjd" in x)]
    expected_data = expected_data.drop(columns=drop_cols)
    written_data = written_data.drop(columns=drop_cols)

    # note that because I'm using Pandas there's some small dtype and np.nan/None stuff to clear up
    # but this makes for a quick streamlined test anyway
    # expected_data = expected_data.replace({np.nan: None}) # TODO: fix weirdness between nan and None?
    written_data = written_data.replace({None: np.nan})  # TODO: fix weirdness between nan and None?
    expected_data = expected_data.astype({"ssObjectId": str})
    written_data = written_data.astype({"ssObjectId": str})

    for x in written_data.columns:
        if type(written_data.iloc[0][x]) == str:
            assert expected_data.iloc[0][x] == written_data.iloc[0][x]
        else:
            assert_almost_equal(expected_data.iloc[0][x], written_data.iloc[0][x])


def test_overwriting_rows_in_database(tmp_path):
    # a test to ensure we're overwriting correctly

    # write the initial object
    db_location = os.path.join(tmp_path, "test_AdlerData_database.db")
    test_object.write_row_to_database(db_location)

    # make a change to a value, then create a new AdlerData object
    # and write that to the database
    test_object_2 = AdlerData("8268570668335894776", ["g"])
    _g_model_1["nobs"] = 666
    test_object_2.populate_phase_parameters("g", **_g_model_1)
    test_object_2.write_row_to_database(db_location)

    con = sqlite3.connect(db_location)
    written_data = pd.read_sql_query("SELECT * from AdlerData", con)
    con.close()

    # should only have one row in the database
    assert len(written_data) == 1

    # we should have a new value for g_nobs but r_nobs should be unchanged
    assert written_data["g_nobs"][0] == 666
    assert written_data["r_nobs"][0] == 38


def test_read_row_from_database():
    # NOTE: the test database here has two rows, one with an earlier timestamp and different data
    # So this test also ensures that only the most recent data for the object is pulled.

    db_location = get_test_data_filepath("test_AdlerData_database.db")

    new_object = AdlerData("8268570668335894776", ["g", "r", "i"])
    new_object.populate_from_database(db_location)

    # print(new_object.__dict__)
    # print(test_object.__dict__)
    # assert new_object.__dict__ == test_object.__dict__

    for filt, model in zip(["g", "g", "r", "i"], [model_1, model_2, model_1, model_2]):
        # print(filt,model)
        new = new_object.get_phase_parameters_in_filter(filt, model)
        # print(new.__dict__)
        test = test_object.get_phase_parameters_in_filter(filt, model)
        # print(test.__dict__)
        # assert new.__dict__ == test.__dict__
        # assert new.__dict__ == pytest.approx(test.__dict__)

        for x in new.__dict__.keys():
            print(x)
            new_val = new.__dict__[x]
            test_val = test.__dict__[x]
            print(new_val, test_val)
            if type(new_val) == str:
                assert new_val == test_val
            else:
                assert_almost_equal(new_val, test_val)

    with pytest.raises(ValueError) as error_info_1:
        empty_data = AdlerData("pretend_object", ["g", "r", "i"])
        empty_data.populate_from_database(db_location)

    assert error_info_1.value.args[0] == "No data found in this database for the supplied ssObjectId."

    with pytest.raises(ValueError) as error_info_2:
        bad_filter = AdlerData("8268570668335894776", ["g", "r", "h"])
        bad_filter.populate_from_database(db_location)

    assert error_info_2.value.args[
        0
    ] == "Data does not exist for some of the requested filters in this database. Filters in database for this object: {}".format(
        test_object.filter_list
    )

    with pytest.raises(ValueError) as error_info_3:
        bad_filter = AdlerData("8268570668335894776", ["g", "r", "h"])
        bad_filter.populate_from_database("./dummy_location.db")

    assert error_info_3.value.args[0] == "Database cannot be found at given filepath."


def test_write_db_dtypes():
    # Tests creating and reading an AdlerData SQL database
    # Importantly, this checks that different python dtypes (np.float32, np.float64) do not cause problems with SQL types (See: https://github.com/lsst-uk/lsst-adler/issues/188)
    db_location = get_test_data_filepath("test_write_AdlerData_database.db")

    # make a test AdlerData object
    filt = "r"
    model_1 = "HG12_Pen16"
    objid = "6098332225018"
    test_object = AdlerData(objid, [filt])
    # results for "adler -s 6098332225018 -n adler_data_61562.db -d 60000.0 61561.5 -np" (see notebooks/adler_demo_rsp)
    filt_model_1 = {
        # 'filter_name': str(filt),
        "phaseAngle_min": np.float32(2.0016675),
        "phaseAngle_range": np.float32(14.333714),
        "nobs": int(29),
        "arc": np.float64(1328.0008499999967),
        "model_name": str(model_1),
        "H": np.float64(16.298209121561705),
        "H_err": np.float64(0.01159924566937738),
        "phase_parameter_1": np.float64(0.634134437903631),
        "phase_parameter_1_err": np.float64(0.06199827430819154),
        "phase_parameter_2": float(0.2),
        "phase_parameter_2_err": None,
    }
    test_object.populate_phase_parameters(filt, **filt_model_1)
    ad = test_object.get_phase_parameters_in_filter(filt, model_1).__dict__

    # delete database file so that it is made from scratch
    if os.path.isfile(db_location):
        os.remove(db_location)

    # create an AdlerData database for this object
    test_object.write_row_to_database(db_location)

    # Check that the database file exists
    assert os.path.isfile(db_location)

    # read the database
    con = sqlite3.connect(db_location)
    written_data = pd.read_sql_query("SELECT * from AdlerData", con)
    con.close()

    # we expect the following data types for each field
    expected_dtypes = {
        "ssObjectId": np.int64,
        "timestamp": np.float64,
        "r_phaseAngle_min": np.float64,
        "r_phaseAngle_range": np.float64,
        "r_nobs": np.int64,
        "r_arc": np.float64,
        "r_HG12_Pen16_H": np.float64,
        "r_HG12_Pen16_H_err": np.float64,
        "r_HG12_Pen16_phase_parameter_1": np.float64,
        "r_HG12_Pen16_phase_parameter_1_err": np.float64,
        "r_HG12_Pen16_phase_parameter_2": np.float64,
        "r_HG12_Pen16_phase_parameter_2_err": object,
        "r_HG12_Pen16_modelFitMjd": object,
    }

    # check the type of each field
    for x in expected_dtypes:
        assert isinstance(written_data.iloc[0][x], expected_dtypes[x])
