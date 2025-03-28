import numpy as np
import pytest

from numpy.testing import assert_array_equal

import adler.utilities.science_utilities as sci_utils

# set up test data residuals
y_res = np.array(
    [
        0.32766394,
        0.15295882,
        0.25896465,
        0.18974566,
        -0.09082071,
        -0.03174331,
        -0.03881347,
        0.04204286,
        -0.24224314,
        -0.17745899,
        -0.00492505,
        -0.0862188,
        -0.14134951,
        -0.2535722,
        -0.07715194,
        -0.10041373,
        -0.24621894,
        -0.02628162,
        -0.04041881,
        0.09704904,
        -0.31718808,
        -0.06039929,
        -0.04072827,
        -0.07548697,
        -0.15773742,
        -0.18076436,
        -0.11133244,
        -0.19906023,
        -0.06507692,
        0.19938055,
        1.97134538,
        1.79108103,
        1.91881725,
        0.20952169,
        0.57060468,
        -0.12014891,
        0.15038993,
        0.63377464,
    ]
)
# uncertainties on y axis parameters
y_err = np.array(
    [
        0.095,
        0.049,
        0.07,
        0.047,
        0.036,
        0.038,
        0.034,
        0.011,
        0.052,
        0.06,
        0.042,
        0.098,
        0.1,
        0.116,
        0.046,
        0.06,
        0.082,
        0.083,
        0.063,
        0.35100001,
        0.2,
        0.057,
        0.034,
        0.036,
        0.054,
        0.051,
        0.046,
        0.17200001,
        0.199,
        0.17,
        0.264,
        0.28299999,
        0.31600001,
        0.163,
        0.248,
        0.182,
        0.211,
        0.38299999,
    ]
)
# provide x data in case residual as a function of x parameter needs tested
x = np.array(
    [
        56.70839691,
        50.07988358,
        72.76689148,
        48.80239868,
        40.30448914,
        39.98022079,
        39.97434235,
        33.56268311,
        38.53084564,
        38.03194427,
        37.81433105,
        36.34851074,
        36.34802628,
        36.50101089,
        32.00559616,
        31.97466469,
        31.86707687,
        31.74268913,
        31.69605255,
        31.88496971,
        30.16068649,
        26.89824295,
        26.91216087,
        27.17892265,
        27.42977715,
        27.42983437,
        27.63732719,
        29.76596832,
        29.76592445,
        27.73990822,
        126.7827301,
        126.78705597,
        126.79136658,
        18.63666534,
        5.73195124,
        7.0533433,
        2.55333257,
        9.17568874,
    ]
)
# flag True for outlying data points
outliers_y = np.array(
    [
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
        True,
        True,
        False,
        False,
        False,
        False,
        False,
    ]
)


def test_outlier_diff():
    outliers = sci_utils.outlier_diff(y_res, diff_cut=1.0)
    assert_array_equal(outliers, outliers_y)


# check that a single value gets processed as a np.array
def test_outlier_diff_single_val():
    outliers = sci_utils.outlier_diff(y_res[0], diff_cut=1.0)
    assert_array_equal(outliers, outliers_y[:1])


# check that a list gets processed as a np.array
def test_outlier_diff_single_val():
    outliers = sci_utils.outlier_diff(list(y_res[:3]), diff_cut=1.0)
    assert_array_equal(outliers, outliers_y[:3])


# test if a single data point is outlying the std residuals of the previous data points
def test_outlier_std():
    outliers = sci_utils.outlier_std(y_res[-8], y_res[:-8])
    assert_array_equal(outliers, outliers_y[-8])


# test the std residual of a list of data points
def test_outlier_std():
    outliers = sci_utils.outlier_std(y_res[-8:-5], y_res[:-8])
    assert_array_equal(outliers, outliers_y[-8:-5])


def test_zero_func():
    assert sci_utils.zero_func(y_res) == 0


def test_sigma_clip():
    outliers = sci_utils.sigma_clip(y_res, sigma=3, maxiters=1, cenfunc=sci_utils.zero_func)
    assert_array_equal(outliers, outliers_y)


def test_outlier_sigma_diff():
    outliers = sci_utils.outlier_sigma_diff(y_res, y_err, std_sigma=5.0)
    assert_array_equal(outliers, outliers_y)


# TODO: test apparition_gap_finder, correct number of groups identified? correct boundaries?
# check boundaires are in ascending order

# TODO: test get_df_obs_filt, check filter and sorting/time range constraints,
# check the column list behaviour
# check AbsMag calculation

# TODO: test running_stats
