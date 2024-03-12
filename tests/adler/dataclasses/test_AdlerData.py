from numpy.testing import assert_array_equal
import pytest

from adler.dataclasses.AdlerData import AdlerData


# setting up the AdlerData object to be used for testing

test_object = AdlerData(["u", "g", "r"])

test_object.populate_phase_parameters("u", "model_1", 11.0, 12.0, 13, 14.0, 15.0, 16.0, [17.0], [18.0])
test_object.populate_phase_parameters(
    "u", "model_2", 21.0, 22.0, 23, 24.0, 25.0, 26.0, [27.0, 27.5], [28.0, 28.5]
)
test_object.populate_phase_parameters("g", "model_1", 31.0, 32.0, 33, 34.0, 35.0, 36.0, [37.0], [38.0])
test_object.populate_phase_parameters(
    "r", "model_2", 41.0, 42.0, 43, 44.0, 45.0, 46.0, [47.0, 47.5], [48.0, 48.5]
)


def test_populate_phase_parameters():
    # test to make sure the object is correctly populated
    assert test_object.filter_list == ["u", "g", "r"]
    assert test_object.model_lists == [["model_1", "model_2"], ["model_1"], ["model_2"]]

    assert_array_equal(test_object.phaseAngle_min_adler, [21.0, 31.0, 41.0])
    assert_array_equal(test_object.phaseAngle_range_adler, [22.0, 32.0, 42.0])
    assert_array_equal(test_object.nobs_adler, [23, 33, 43])
    assert_array_equal(test_object.arc_adler, [24.0, 34.0, 44.0])

    assert test_object.H_adler == [[15.0, 25.0], [35.0], [45.0]]
    assert test_object.H_err_adler == [[16.0, 26.0], [36.0], [46.0]]
    assert test_object.phase_parameters == [[[17.0], [27.0, 27.5]], [[37.0]], [[47.0, 47.5]]]
    assert test_object.phase_parameters_err == [[[18.0], [28.0, 28.5]], [[38.0]], [[48.0, 48.5]]]

    # testing to make sure the correct error messages trigger
    with pytest.raises(TypeError) as error_info_1:
        test_object.populate_phase_parameters("u", "model_1", 11.0, 12.0, 13, 14.0, 15.0, 16.0, 17.0, 18.0)

    assert error_info_1.value.args[0] == "Both parameters and parameters_err arguments must be lists."

    with pytest.raises(Exception) as error_info_2:
        test_object.populate_phase_parameters("y", "model_1", 11.0, 12.0, 13, 14.0, 15.0, 16.0, 17.0, 18.0)

    assert error_info_2.value.args[0] == "Filter y does not exist in AdlerData.filter_list."


def test_get_phase_parameters_in_filter():
    u_model2 = {
        "phaseAngle_min": 21.0,
        "phaseAngle_range": 22.0,
        "nobs": 23,
        "arc": 24.0,
        "H": 25.0,
        "H_err": 26.0,
        "phase_parameters": [27.0, 27.5],
        "phase_parameters_err": [28.0, 28.5],
    }

    u_model1 = {
        "phaseAngle_min": 21.0,
        "phaseAngle_range": 22.0,
        "nobs": 23,
        "arc": 24.0,
        "H": 15.0,
        "H_err": 16.0,
        "phase_parameters": [17.0],
        "phase_parameters_err": [18.0],
    }

    # making sure the correct parameters are retreived
    assert test_object.get_phase_parameters_in_filter("u", model_name="model_2") == u_model2
    assert test_object.get_phase_parameters_in_filter("u") == u_model1
    assert test_object.get_phase_parameters_in_filter(
        "r", "model_2"
    ) == test_object.get_phase_parameters_in_filter("r")

    # checking the error messages
    with pytest.raises(ValueError) as error_info_1:
        error_dict = test_object.get_phase_parameters_in_filter("f", model_name="model_2")

    assert error_info_1.value.args[0] == "Filter f does not exist in AdlerData.filter_list."

    with pytest.raises(ValueError) as error_info_2:
        error_dict_2 = test_object.get_phase_parameters_in_filter("r", model_name="model_1")

    assert error_info_2.value.args[0] == "Model model_1 does not exist for filter r in AdlerData.model_lists."
