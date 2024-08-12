import os
import pytest
from adler.utilities.AdlerCLIArguments import AdlerCLIArguments
from adler.utilities.tests_utilities import get_test_data_filepath


# AdlerCLIArguments object takes an object as input, so we define a quick one here
class args:
    def __init__(
        self,
        ssObjectId,
        ssObjectId_list,
        filter_list,
        colour_list,
        date_range,
        outpath,
        db_name,
        sql_filename,
    ):
        self.ssObjectId = ssObjectId
        self.ssObjectId_list = ssObjectId_list
        self.filter_list = filter_list
        self.colour_list = colour_list
        self.date_range = date_range
        self.outpath = outpath
        self.db_name = db_name
        self.sql_filename = sql_filename


def test_AdlerCLIArguments_population():
    # test correct population
    good_input_dict = {
        "ssObjectId": "666",
        "ssObjectId_list": None,
        "filter_list": ["g", "r", "i"],
        "colour_list": ["g-r", "r-i"],
        "date_range": [60000.0, 67300.0],
        "outpath": "./",
        "db_name": "output",
        "sql_filename": None,
    }
    good_arguments = args(**good_input_dict)
    good_arguments_object = AdlerCLIArguments(good_arguments)

    good_input_dict["outpath"] = os.path.abspath("./")

    assert good_arguments_object.__dict__ == good_input_dict

    # and double-check that ssObjectID_list works too

    good_input_dict["ssObjectId_list"] = get_test_data_filepath("test_SSOIDs.txt")
    good_arguments = args(**good_input_dict)
    good_arguments_object = AdlerCLIArguments(good_arguments)
    assert good_arguments_object.ssObjectId_list == get_test_data_filepath("test_SSOIDs.txt")

    # also test setting sql_filename

    good_input_dict["sql_filename"] = get_test_data_filepath("testing_database.db")
    good_arguments = args(**good_input_dict)
    good_arguments_object = AdlerCLIArguments(good_arguments)
    assert good_arguments_object.sql_filename == get_test_data_filepath("testing_database.db")


def test_AdlerCLIArguments_badSSOID():
    # test that a bad ssObjectId triggers the right error
    bad_ssoid_arguments = args(
        "hello!", None, ["g", "r", "i"], ["g-r", "r-i"], [60000.0, 67300.0], "./", "output", "None"
    )

    with pytest.raises(ValueError) as bad_ssoid_error:
        bad_ssoid_object = AdlerCLIArguments(bad_ssoid_arguments)

    assert (
        bad_ssoid_error.value.args[0]
        == "--ssObjectId command-line argument does not appear to be a valid ssObjectId."
    )


def test_AdlerCLIArguments_badfilters():
    # test that non-LSST or unexpected filters trigger the right error
    bad_filter_arguments = args(
        "666", None, ["g", "r", "i", "m"], ["g-r", "r-i"], [60000.0, 67300.0], "./", "output", None
    )

    with pytest.raises(ValueError) as bad_filter_error:
        bad_filter_object = AdlerCLIArguments(bad_filter_arguments)

    assert (
        bad_filter_error.value.args[0]
        == "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
    )

    bad_filter_arguments_2 = args(
        "666", None, ["pony"], ["g-r", "r-i"], [60000.0, 67300.0], "./", "output", None
    )

    with pytest.raises(ValueError) as bad_filter_error_2:
        bad_filter_object = AdlerCLIArguments(bad_filter_arguments_2)

    assert (
        bad_filter_error_2.value.args[0]
        == "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
    )


# TODO: bad colour filters, e.g. ["g-r", "r-j"] or ["g-r", "r - y"] or ["g-r", "rxi"]
# also test that correct filter_list has been requested


def test_AdlerCLIArguments_baddates():
    # test that overly-large dates trigger the right error
    big_date_arguments = args(
        "666", None, ["g", "r", "i"], ["g-r", "r-i"], [260000.0, 267300.0], "./", "output", None
    )

    with pytest.raises(ValueError) as big_date_error:
        big_date_object = AdlerCLIArguments(big_date_arguments)

    assert (
        big_date_error.value.args[0]
        == "Dates for --date_range command-line argument seem rather large. Did you input JD instead of MJD?"
    )

    # test that unexpected date values trigger the right error
    bad_date_arguments = args(
        "666", None, ["g", "r", "i"], ["g-r", "r-i"], [60000.0, "cheese"], "./", "output", None
    )

    with pytest.raises(ValueError) as bad_date_error:
        bad_date_object = AdlerCLIArguments(bad_date_arguments)

    assert (
        bad_date_error.value.args[0]
        == "One or both of the values for the --date_range command-line argument do not seem to be valid numbers."
    )


def test_AdlerCLIArguments_badoutput():
    bad_output_arguments = args(
        "666",
        None,
        ["g", "r", "i"],
        ["g-r", "r-i"],
        [60000.0, 67300.0],
        "./definitely_fake_folder/",
        "output",
        None,
    )

    with pytest.raises(ValueError) as bad_output_error:
        bad_output_object = AdlerCLIArguments(bad_output_arguments)

    assert (
        bad_output_error.value.args[0]
        == "The output path for the command-line argument --outpath cannot be found."
    )


def test_AdlerCLIArguments_badlist():
    bad_list_arguments = args(
        None,
        "./fake_input/here.txt",
        ["g", "r", "i"],
        ["g-r", "r-i"],
        [60000.0, 67300.0],
        "./",
        "output",
        "./",
    )

    with pytest.raises(ValueError) as bad_list_error:
        bad_list_object = AdlerCLIArguments(bad_list_arguments)

    assert (
        bad_list_error.value.args[0]
        == "The file supplied for the command-line argument --ssObjectId_list cannot be found."
    )


def test_AdlerCLIArguments_badsql():
    bad_sql_arguments = args(
        "666",
        None,
        ["g", "r", "i"],
        ["g-r", "r-i"],
        [60000.0, 67300.0],
        "./",
        "output",
        "./dummy_database.db",
    )

    with pytest.raises(ValueError) as bad_sql_error:
        bad_sql_object = AdlerCLIArguments(bad_sql_arguments)

    assert (
        bad_sql_error.value.args[0]
        == "The file supplied for the command-line argument --sql_filename cannot be found."
    )
