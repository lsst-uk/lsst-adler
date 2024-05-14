import pytest

from adler.utilities.tests_utilities import get_test_data_filepath


def test_read_in_SSObjectID_file():

    from adler.utilities.readin_utilities import read_in_SSObjectID_file

    good_ssoids = read_in_SSObjectID_file(get_test_data_filepath("test_SSOIDs.txt"))

    expected = ["2150553186630", "3369984299447", "5992863104062", "6098332225018", "6102997768245"]

    assert good_ssoids == expected

    with pytest.raises(ValueError) as error_info:
        bad_ssoids = read_in_SSObjectID_file(get_test_data_filepath("test_SSOIDs_bad.txt"))

    assert (
        error_info.value.args[0]
        == "One or more of the SSObjectIDs in the supplied list does not seem to be valid."
    )
