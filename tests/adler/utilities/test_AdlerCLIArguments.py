import os
import pytest
from adler.utilities.AdlerCLIArguments import AdlerCLIArguments


# AdlerCLIArguments object takes an object as input, so we define a quick one here
class args:
    def __init__(self, ssObjectId, filter_list, date_range, outpath, db_name):
        self.ssObjectId = ssObjectId
        self.filter_list = filter_list
        self.date_range = date_range
        self.outpath = outpath
        self.db_name = db_name


def test_AdlerCLIArguments_population():
    # test correct population
    good_input_dict = {
        "ssObjectId": "666",
        "filter_list": ["g", "r", "i"],
        "date_range": [60000.0, 67300.0],
        "outpath": "./",
        "db_name": "output",
    }
    good_arguments = args(**good_input_dict)
    good_arguments_object = AdlerCLIArguments(good_arguments)

    good_input_dict["outpath"] = os.path.abspath("./")

    assert good_arguments_object.__dict__ == good_input_dict


def test_AdlerCLIArguments_badSSOID():
    # test that a bad ssObjectId triggers the right error
    bad_ssoid_arguments = args("hello!", ["g", "r", "i"], [60000.0, 67300.0], "./", "output")

    with pytest.raises(ValueError) as bad_ssoid_error:
        bad_ssoid_object = AdlerCLIArguments(bad_ssoid_arguments)

    assert (
        bad_ssoid_error.value.args[0]
        == "--ssObjectId command-line argument does not appear to be a valid ssObjectId."
    )


def test_AdlerCLIArguments_badfilters():
    # test that non-LSST or unexpected filters trigger the right error
    bad_filter_arguments = args("666", ["g", "r", "i", "m"], [60000.0, 67300.0], "./", "output")

    with pytest.raises(ValueError) as bad_filter_error:
        bad_filter_object = AdlerCLIArguments(bad_filter_arguments)

    assert (
        bad_filter_error.value.args[0]
        == "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
    )

    bad_filter_arguments_2 = args("666", ["pony"], [60000.0, 67300.0], "./", "output")

    with pytest.raises(ValueError) as bad_filter_error_2:
        bad_filter_object = AdlerCLIArguments(bad_filter_arguments_2)

    assert (
        bad_filter_error_2.value.args[0]
        == "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
    )


def test_AdlerCLIArguments_baddates():
    # test that overly-large dates trigger the right error
    big_date_arguments = args("666", ["g", "r", "i"], [260000.0, 267300.0], "./", "output")

    with pytest.raises(ValueError) as big_date_error:
        big_date_object = AdlerCLIArguments(big_date_arguments)

    assert (
        big_date_error.value.args[0]
        == "Dates for --date_range command-line argument seem rather large. Did you input JD instead of MJD?"
    )

    # test that unexpected date values trigger the right error
    bad_date_arguments = args("666", ["g", "r", "i"], [60000.0, "cheese"], "./", "output")

    with pytest.raises(ValueError) as bad_date_error:
        bad_date_object = AdlerCLIArguments(bad_date_arguments)

    assert (
        bad_date_error.value.args[0]
        == "One or both of the values for the --date_range command-line argument do not seem to be valid numbers."
    )


def test_AdlerCLIArguments_badoutput():
    bad_output_arguments = args(
        "666", ["g", "r", "i"], [60000.0, 67300.0], "./definitely_fake_folder/", "output"
    )

    with pytest.raises(ValueError) as bad_output_error:
        bad_output_object = AdlerCLIArguments(bad_output_arguments)

    assert (
        bad_output_error.value.args[0]
        == "The output path for the command-line argument --outpath cannot be found."
    )
