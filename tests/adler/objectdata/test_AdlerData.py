import os
import sqlite3
import pandas as pd
import numpy as np
import pytest

from numpy.testing import assert_equal, assert_array_equal, assert_almost_equal, assert_array_almost_equal

from adler.objectdata.AdlerData import AdlerData, AdlerSourceFlags, VALID_PHASE_MODELS, VALID_AVG_MAG_MODELS
from adler.utilities.tests_utilities import get_test_data_filepath

mpc_model_name_in = "median"
mpc_process_mjd_in = 60799.5
mpc_data_timespan_in = 7
mpc_n_new_nights_in = 3
db_path_mpc = get_test_data_filepath(
    f"adler_output_{mpc_model_name_in}_{mpc_process_mjd_in:.1f}_{mpc_data_timespan_in}n_{mpc_n_new_nights_in}n.sqlite"
)
ssoid_mpc = "2025 MX40Fake"
filter_list_mpc = ["g", "r"]
mpc_expected_modelId = "median_60799.5_7n_3n"

g_model_mpc = {
    "filter_name": "g",
    "phaseAngle_min": 10.4569,
    "phaseAngle_range": 0.5461,
    "observationTime_max": 60795.154395300924,
    "nobs": 20,
    "arc": 2.018855601854739,
    "n_outliers": 13,
    "n_std_outliers": 16,
    "sustained_outliers": 1.6170569944178368,
    "model_name": "median",
    "avg_mag": 19.905139516804923,
    "std_mag": 0.1293180601923424,
}
r_model_mpc = {
    "filter_name": "r",
    "phaseAngle_min": 10.4703,
    "phaseAngle_range": 0.5175,
    "observationTime_max": 60795.09814260416,
    "nobs": 20,
    "arc": 1.9146494097221876,
    "n_outliers": 14,
    "n_std_outliers": 42,
    "sustained_outliers": np.nan,
    "model_name": "median",
    "avg_mag": 19.371997673855667,
    "std_mag": 0.14717178763567726,
}

g_df_mpc = pd.DataFrame(
    {
        "diaSourceId": [
            "LorDdACC0000Gugo010000qGt",
            "LorDdACC0000Gugo010000qGu",
            "LorDdACC0000Gugo010000qGv",
            "LorDdACC0000Gugo010000qGx",
            "LorDdACC0000Gugo010000qGy",
            "LorDdACC0000Gugo010000qGz",
            "LorDdACC0000Gugo010000qH0",
            "LorDdACC0000Gugo010000qH1",
            "LorDdACC0000Gugo010000qH2",
            "LorDdACC0000Gugo010000qH3",
            "LorDdACC0000Gugo010000qH7",
            "LorDdACC0000Gugo010000qH8",
            "LorDdACC0000Gugo010000qH9",
            "LorDdACC0000Gugo010000qHT",
            "LorDdACC0000Gugo010000qHU",
            "LorDdACC0000Gugo010000qHV",
        ],
        "midpointMjdTai": [
            60797.090522,
            60797.095330,
            60797.097267,
            60797.100167,
            60797.101139,
            60797.102597,
            60797.110885,
            60797.111853,
            60797.122101,
            60797.123073,
            60798.105471,
            60798.108366,
            60798.112740,
            60799.150561,
            60799.152043,
            60799.153052,
        ],
        "mag_diff": [
            -1.668184,
            -1.500223,
            -1.537238,
            0.0,
            -1.568270,
            0.0,
            0.0,
            -1.758356,
            -3.543438,
            -3.634446,
            -1.920408,
            -1.745432,
            -2.092468,
            -1.576040,
            -1.587053,
            -1.647061,
        ],
        "std_diff": [
            -10.974898,
            -7.463796,
            -8.540214,
            -6.429907,
            -10.053010,
            -7.144345,
            -7.656432,
            -12.126591,
            -18.171478,
            -19.860360,
            -9.650293,
            -8.432038,
            -11.560596,
            -7.540863,
            -7.557394,
            -8.034445,
        ],
    }
)

r_df_mpc = pd.DataFrame(
    {
        "diaSourceId": [
            "LorDdACC0000Gugo010000qGX",
            "LorDdACC0000Gugo010000qGY",
            "LorDdACC0000Gugo010000qGZ",
            "LorDdACC0000Gugo010000qGa",
            "LorDdACC0000Gugo010000qGb",
            "LorDdACC0000Gugo010000qGc",
            "LorDdACC0000Gugo010000qGd",
            "LorDdACC0000Gugo010000qGe",
            "LorDdACC0000Gugo010000qGf",
            "LorDdACC0000Gugo010000qGg",
            "LorDdACC0000Gugo010000qGh",
            "LorDdACC0000Gugo010000qGi",
            "LorDdACC0000Gugo010000qGj",
            "LorDdACC0000Gugo010000qGk",
            "LorDdACC0000Gugo010000qGl",
            "LorDdACC0000Gugo010000qGm",
            "LorDdACC0000Gugo010000qGn",
            "LorDdACC0000Gugo010000qGo",
            "LorDdACC0000Gugo010000qGp",
            "LorDdACC0000Gugo010000qGq",
            "LorDdACC0000Gugo010000qGr",
            "LorDdACC0000Gugo010000qGs",
            "LorDdACC0000Gugo010000qH5",
            "LorDdACC0000Gugo010000qH6",
            "LorDdACC0000Gugo010000qHA",
            "LorDdACC0000Gugo010000qHB",
            "LorDdACC0000Gugo010000qHC",
            "LorDdACC0000Gugo010000qHD",
            "LorDdACC0000Gugo010000qHE",
            "LorDdACC0000Gugo010000qHF",
            "LorDdACC0000Gugo010000qHG",
            "LorDdACC0000Gugo010000qHH",
            "LorDdACC0000Gugo010000qHJ",
            "LorDdACC0000Gugo010000qHK",
            "LorDdACC0000Gugo010000qHL",
            "LorDdACC0000Gugo010000qHM",
            "LorDdACC0000Gugo010000qHN",
            "LorDdACC0000Gugo010000qHO",
            "LorDdACC0000Gugo010000qHP",
            "LorDdACC0000Gugo010000qHQ",
            "LorDdACC0000Gugo010000qHR",
            "LorDdACC0000Gugo010000qHS",
        ],
        "midpointMjdTai": [
            60797.055184,
            60797.060023,
            60797.061476,
            60797.062450,
            60797.063412,
            60797.064858,
            60797.065812,
            60797.066294,
            60797.067287,
            60797.068249,
            60797.069794,
            60797.071229,
            60797.072190,
            60797.073152,
            60797.073638,
            60797.074592,
            60797.075564,
            60797.076049,
            60797.077037,
            60797.078029,
            60797.080465,
            60797.080953,
            60798.081340,
            60798.089184,
            60799.088375,
            60799.093265,
            60799.095768,
            60799.096772,
            60799.098220,
            60799.100829,
            60799.101795,
            60799.110593,
            60799.115549,
            60799.120395,
            60799.121356,
            60799.123822,
            60799.125279,
            60799.131254,
            60799.133702,
            60799.135169,
            60799.137626,
            60799.138589,
        ],
        "mag_diff": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -1.596888,
            0.0,
            -1.660904,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -1.794069,
            -1.633133,
            -1.674378,
            -1.572419,
            -1.678440,
            0.0,
            -1.514460,
            -1.702482,
            -1.537490,
            -1.652564,
            0.0,
            -1.526646,
            -1.519654,
            0.0,
            0.0,
            0.0,
            -1.517757,
            0.0,
            0.0,
            0.0,
        ],
        "std_diff": [
            -8.021991,
            -7.880467,
            -9.132070,
            -6.933433,
            -8.403962,
            -8.752083,
            -8.922084,
            -9.349993,
            -5.755856,
            -9.319376,
            -7.538622,
            -12.097639,
            -9.056304,
            -11.072691,
            -7.809539,
            -6.918899,
            -7.073666,
            -9.564460,
            -9.519307,
            -7.488570,
            -6.590224,
            -11.014782,
            -10.678983,
            -7.814034,
            -9.250710,
            -8.408659,
            -9.273149,
            -7.036157,
            -8.460673,
            -9.728470,
            -8.401586,
            -10.393482,
            -7.439000,
            -9.196660,
            -9.266180,
            -8.718299,
            -6.616110,
            -6.706031,
            -8.824168,
            -6.857420,
            -7.009274,
            -7.490596,
        ],
    }
)


dp03_model_name_in = "HG12_Pen16"
dp03_process_mjd_in = 63335.5
dp03_data_timespan_in = 400
dp03_n_new_nights_in = 31
db_path_dp03 = get_test_data_filepath(
    f"adler_output_{dp03_model_name_in}_{dp03_process_mjd_in:.1f}_{dp03_data_timespan_in}n_{dp03_n_new_nights_in}n.sqlite"
)
ssoid_dp03 = "6098332225018000"
filter_list_dp03 = ["r", "i"]
dp03_expected_modelId = "HG12_Pen16_63335.5_400n_31n"

r_model_dp03 = {
    "filter_name": "r",
    "phaseAngle_min": 9.478214263916016,
    "phaseAngle_range": 10.772890090942383,
    "observationTime_max": 63001.97873,
    "nobs": 9,
    "arc": 62.78039,
    "n_outliers": 1,
    "n_std_outliers": 6,
    "sustained_outliers": np.nan,
    "model_name": "HG12_Pen16",
    "H": 16.299738958850803,
    "H_err": 0.008126418528902073,
    "phase_parameter_1": 0.7074497288054502,
    "phase_parameter_1_err": 0.09699260052621479,
    "phase_parameter_2": np.nan,
    "phase_parameter_2_err": np.nan,
}
i_model_dp03 = {
    "filter_name": "i",
    "phaseAngle_min": 8.388969421386719,
    "phaseAngle_range": 11.862516403198242,
    "observationTime_max": 63019.97406,
    "nobs": 16,
    "arc": 83.8726,
    "n_outliers": 2,
    "n_std_outliers": 2,
    "sustained_outliers": np.nan,
    "model_name": "HG12_Pen16",
    "H": 16.180738503756075,
    "H_err": 0.007805159809346133,
    "phase_parameter_1": 0.524322463314371,
    "phase_parameter_1_err": 0.08594588017371281,
    "phase_parameter_2": np.nan,
    "phase_parameter_2_err": np.nan,
}

r_df_dp03 = pd.DataFrame(
    {
        "diaSourceId": [
            -7402261812818583560,
            -6064878124717393944,
            -7827113376855357760,
            -42385174921712352,
            6104500740852729128,
            5185229085827158712,
        ],
        "midpointMjdTai": [
            63305.34722,
            63305.34774,
            63308.28513,
            63308.28558,
            63316.24186,
            63335.34367,
        ],
        "mag_diff": [
            -1.506460,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        "std_diff": [
            -48.595488,
            -44.559145,
            -39.589294,
            -41.937080,
            -29.603367,
            -70.110469,
        ],
    }
)

i_df_dp03 = pd.DataFrame(
    {
        "diaSourceId": [
            -7564747308461840760,
            3495766367455214856,
        ],
        "midpointMjdTai": [
            63316.26713,
            63335.36783,
        ],
        "mag_diff": [
            -1.570527,
            -1.503535,
        ],
        "std_diff": [
            -39.263165,
            -60.141398,
        ],
    }
)

# Create versions with no model name
_g_model_mpc = g_model_mpc.copy()
del _g_model_mpc["filter_name"]
_r_model_mpc = r_model_mpc.copy()
del _r_model_mpc["filter_name"]
_r_model_dp03 = r_model_dp03.copy()
del _r_model_dp03["filter_name"]
_i_model_dp03 = i_model_dp03.copy()
del _i_model_dp03["filter_name"]


def test_set_modelId():
    ad_mpc = AdlerData(ssoid_mpc, filter_list_mpc)
    ad_mpc.set_modelId(
        model_name=mpc_model_name_in,
        end_mjd=mpc_process_mjd_in,
        data_timespan=mpc_data_timespan_in,
        n_new_nights=mpc_n_new_nights_in,
    )

    assert ad_mpc.modelId == mpc_expected_modelId


mpc_test_obj = AdlerData(ssoid_mpc, filter_list_mpc)
mpc_test_obj.populate_filter_dependent_parameters("g", **_g_model_mpc)
mpc_test_obj.populate_filter_dependent_parameters("r", **_r_model_mpc)


def test_populate_filter_dependent_parameters():
    # test to make sure the object is correctly populated
    assert mpc_test_obj.filter_list == ["g", "r"]

    assert_almost_equal(
        [a.phaseAngle_min for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["phaseAngle_min"], r_model_mpc["phaseAngle_min"]],
    )

    assert_almost_equal(
        [a.phaseAngle_range for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["phaseAngle_range"], r_model_mpc["phaseAngle_range"]],
    )

    assert_almost_equal(
        [a.observationTime_max for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["observationTime_max"], r_model_mpc["observationTime_max"]],
    )

    assert_equal(
        [a.nobs for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["nobs"], r_model_mpc["nobs"]],
    )

    assert_almost_equal(
        [a.arc for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["arc"], r_model_mpc["arc"]],
    )

    assert_equal(
        [a.n_outliers for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["n_outliers"], r_model_mpc["n_outliers"]],
    )

    assert_equal(
        [a.n_std_outliers for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["n_std_outliers"], r_model_mpc["n_std_outliers"]],
    )

    assert_almost_equal(
        [a.sustained_outliers for a in mpc_test_obj.filter_dependent_values],
        [g_model_mpc["sustained_outliers"], r_model_mpc["sustained_outliers"]],
    )

    # check to make sure filter-dependent parameter is correctly updated (then return it to previous)
    mpc_test_obj.populate_avg_mag_parameters("g", model_name=mpc_model_name_in, nobs=999)
    assert mpc_test_obj.get_avg_mag_parameters_in_filter("g", mpc_model_name_in).nobs == 999
    mpc_test_obj.populate_avg_mag_parameters("g", model_name=mpc_model_name_in, nobs=g_model_mpc["nobs"])

    # testing to make sure the correct error messages trigger
    test_filt = "y"
    with pytest.raises(ValueError) as error_info_1:
        mpc_test_obj.populate_phase_parameters(test_filt)

    assert error_info_1.value.args[0] == "Filter {} does not exist in AdlerData.filter_list.".format(
        test_filt
    )


dp03_test_obj = AdlerData(ssoid_dp03, ["r", "i"])
dp03_test_obj.set_modelId(
    model_name=dp03_model_name_in,
    end_mjd=dp03_process_mjd_in,
    data_timespan=dp03_data_timespan_in,
    n_new_nights=dp03_n_new_nights_in,
)
dp03_test_obj.populate_phase_parameters("r", **_r_model_dp03)
dp03_test_obj.populate_phase_parameters("i", **_i_model_dp03)


def test_populate_phase_parameters():
    test_dict_list = [
        dp03_test_obj.get_phase_parameters_in_filter("r", dp03_model_name_in).__dict__,
        dp03_test_obj.get_phase_parameters_in_filter("i", dp03_model_name_in).__dict__,
    ]
    expect_dict_list = [r_model_dp03, i_model_dp03]

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
    dp03_test_obj.populate_phase_parameters("r", model_name=dp03_model_name_in, H=99.0)
    assert dp03_test_obj.get_phase_parameters_in_filter("r", dp03_model_name_in).H == 99.0
    dp03_test_obj.populate_phase_parameters("r", model_name=dp03_model_name_in, H=r_model_dp03["H"])

    # testing to make sure the correct error messages trigger
    test_filt = "y"
    with pytest.raises(ValueError) as error_info_1:
        dp03_test_obj.populate_phase_parameters(test_filt)

    assert error_info_1.value.args[0] == "Filter {} does not exist in AdlerData.filter_list.".format(
        test_filt
    )

    with pytest.raises(NameError) as error_info_2:
        dp03_test_obj.populate_phase_parameters("r", H=4.0)

    assert error_info_2.value.args[0] == "No model name given. Cannot update model-specific phase parameters."


mpc_test_obj = AdlerData(ssoid_mpc, filter_list_mpc)
mpc_test_obj.populate_avg_mag_parameters("g", **_g_model_mpc)
mpc_test_obj.populate_avg_mag_parameters("r", **_r_model_mpc)


def test_populate_avg_mag_parameters():
    test_dict_list = [
        mpc_test_obj.get_avg_mag_parameters_in_filter("g", mpc_model_name_in).__dict__,
        mpc_test_obj.get_avg_mag_parameters_in_filter("r", mpc_model_name_in).__dict__,
    ]
    expect_dict_list = [g_model_mpc, r_model_mpc]

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
    mpc_test_obj.populate_avg_mag_parameters("r", model_name=mpc_model_name_in, avg_mag=99.0)
    assert mpc_test_obj.get_avg_mag_parameters_in_filter("r", mpc_model_name_in).avg_mag == 99.0
    mpc_test_obj.populate_avg_mag_parameters(
        "r", model_name=mpc_model_name_in, avg_mag=r_model_mpc["avg_mag"]
    )

    # testing to make sure the correct error messages trigger
    test_filt = "y"
    with pytest.raises(ValueError) as error_info_1:
        mpc_test_obj.populate_avg_mag_parameters(test_filt)

    assert error_info_1.value.args[0] == "Filter {} does not exist in AdlerData.filter_list.".format(
        test_filt
    )

    with pytest.raises(NameError) as error_info_2:
        mpc_test_obj.populate_avg_mag_parameters("r", avg_mag=4.0)

    assert (
        error_info_2.value.args[0]
        == "No model name given. Cannot update model-specific average magnitude parameters."
    )


mpc_source_flags_obj = AdlerSourceFlags.construct_source_flags_from_data_table(
    ssoid_mpc, "g", mpc_expected_modelId, g_df_mpc
)


def test_construct_source_flags_from_data_table():
    # Test that the constructed AdlerSourceFlags object has the correct attributes
    assert mpc_source_flags_obj.ssObjectId == ssoid_mpc
    assert mpc_source_flags_obj.filter_name == "g"
    assert mpc_source_flags_obj.modelId == mpc_expected_modelId
    assert mpc_source_flags_obj.n_outliers == g_model_mpc["n_outliers"]
    assert mpc_source_flags_obj.n_std_outliers == g_model_mpc["n_std_outliers"]

    # Test that the arrays are correctly populated from the DataFrame
    assert_array_equal(mpc_source_flags_obj.diaSourceId, g_df_mpc["diaSourceId"].values)
    assert_array_almost_equal(mpc_source_flags_obj.midpointMjdTai, g_df_mpc["midpointMjdTai"].values)
    assert_array_almost_equal(mpc_source_flags_obj.mag_diff, g_df_mpc["mag_diff"].values)
    assert_array_almost_equal(mpc_source_flags_obj.std_diff, g_df_mpc["std_diff"].values)


mpc_test_obj.set_modelId(
    model_name=mpc_model_name_in,
    end_mjd=mpc_process_mjd_in,
    data_timespan=mpc_data_timespan_in,
    n_new_nights=mpc_n_new_nights_in,
)
mpc_test_obj.populate_source_flags("r", mpc_expected_modelId, r_df_mpc)


def test_populate_source_flags():
    # Test that source flags are correctly populated in the AdlerData object

    # Verify that the source_flags object was created and stored in the correct filter
    assert mpc_test_obj.filter_dependent_values[1].source_flags is not None
    assert isinstance(mpc_test_obj.filter_dependent_values[1].source_flags, AdlerSourceFlags)

    # Retrieve the source flags from the r filter
    r_source_flags = mpc_test_obj.filter_dependent_values[1].source_flags

    # Test the basic attributes match expected values
    assert r_source_flags.ssObjectId == ssoid_mpc
    assert r_source_flags.filter_name == "r"
    assert r_source_flags.modelId == mpc_expected_modelId

    # Calculate expected counts from DataFrame
    expected_n_outliers = len(r_df_mpc.loc[r_df_mpc.mag_diff != 0])
    expected_n_std_outliers = len(r_df_mpc.loc[r_df_mpc.std_diff != 0])

    assert r_source_flags.n_outliers == expected_n_outliers
    assert r_source_flags.n_std_outliers == expected_n_std_outliers

    # Test that the arrays are correctly populated from the DataFrame
    assert_array_equal(r_source_flags.diaSourceId, r_df_mpc["diaSourceId"].values)
    assert_array_almost_equal(r_source_flags.midpointMjdTai, r_df_mpc["midpointMjdTai"].values)
    assert_array_almost_equal(r_source_flags.mag_diff, r_df_mpc["mag_diff"].values)
    assert_array_almost_equal(r_source_flags.std_diff, r_df_mpc["std_diff"].values)

    # Test that filter-dependent parameters were also updated with correct outlier counts
    assert mpc_test_obj.filter_dependent_values[1].n_outliers == expected_n_outliers
    assert mpc_test_obj.filter_dependent_values[1].n_std_outliers == expected_n_std_outliers

    # testing to make sure the correct error messages trigger
    test_modelId = "wrong_modelId"
    with pytest.raises(ValueError) as error_info_1:
        mpc_test_obj.populate_source_flags("r", test_modelId, r_df_mpc)

    assert (
        error_info_1.value.args[0]
        == f"modelId {test_modelId} does not match the modelId in AdlerData.modelId: {mpc_test_obj.modelId}"
    )

    # testing to make sure the correct error messages trigger
    test_filt = "y"
    with pytest.raises(ValueError) as error_info_2:
        mpc_test_obj.populate_source_flags(test_filt, mpc_expected_modelId, r_df_mpc)

    assert error_info_2.value.args[0] == "Filter {} does not exist in AdlerData.filter_list.".format(
        test_filt
    )


# here the capsys fixture captures any output to the terminal
def test_print_data(capsys):
    dp03_test_obj.print_data()

    # get what was printed to the terminal
    captured = capsys.readouterr()

    expected = "Filter: r\nPhase angle minimum: 9.478214263916016\nPhase angle range: 10.772890090942383\nMaximum observation time: 63001.97873\nNumber of observations: 9\nArc: 62.78039\nNumber of outliers detected: 1\nNumber of outliers in sigma-space detected: 6\nMagnitude change of sustained outliers: nan\nModel: HG12_Pen16.\nH: 16.299738958850803\nH error: 0.008126418528902073\nPhase parameter 1: 0.7074497288054502\nPhase parameter 1 error: 0.09699260052621479\nPhase parameter 2: nan\nPhase parameter 2 error: nan\n\n\nFilter: i\nPhase angle minimum: 8.388969421386719\nPhase angle range: 11.862516403198242\nMaximum observation time: 63019.97406\nNumber of observations: 16\nArc: 83.8726\nNumber of outliers detected: 2\nNumber of outliers in sigma-space detected: 2\nMagnitude change of sustained outliers: nan\nModel: HG12_Pen16.\nH: 16.180738503756075\nH error: 0.007805159809346133\nPhase parameter 1: 0.524322463314371\nPhase parameter 1 error: 0.08594588017371281\nPhase parameter 2: nan\nPhase parameter 2 error: nan\n\n\n"

    # print(captured)
    # print(expected)
    assert captured.out == expected


def test_get_model_name():
    ad_mpc = AdlerData(ssoid_mpc, filter_list_mpc)
    ad_mpc.set_modelId(
        model_name=mpc_model_name_in,
        end_mjd=mpc_process_mjd_in,
        data_timespan=mpc_data_timespan_in,
        n_new_nights=mpc_n_new_nights_in,
    )

    model_name_mpc = ad_mpc._get_model_name()

    assert model_name_mpc == mpc_model_name_in

    ad_dp03 = AdlerData(ssoid_dp03, filter_list_dp03)
    ad_dp03.set_modelId(
        model_name=dp03_model_name_in,
        end_mjd=dp03_process_mjd_in,
        data_timespan=dp03_data_timespan_in,
        n_new_nights=dp03_n_new_nights_in,
    )

    model_name_dp03 = ad_dp03._get_model_name()

    assert model_name_dp03 == dp03_model_name_in

    # test correct error is raised
    ad_dp03_wrongmodel = AdlerData(ssoid_dp03, filter_list_dp03)
    dp03_wrongmodel_name = "foo"
    ad_dp03_wrongmodel.set_modelId(
        model_name=dp03_wrongmodel_name,
        end_mjd=dp03_process_mjd_in,
        data_timespan=dp03_data_timespan_in,
        n_new_nights=dp03_n_new_nights_in,
    )

    with pytest.raises(ValueError) as error_info_1:
        ad_dp03_wrongmodel._get_model_name()

    assert (
        error_info_1.value.args[0]
        == f"Unknown model in string: {dp03_wrongmodel_name}_{dp03_process_mjd_in:.1f}_{dp03_data_timespan_in}n_{dp03_n_new_nights_in}n"
    )


def generate_expected_adlerdata_csv(adler_obj):
    """
    Generate expected CSV data from an AdlerData object for comparison during database tests.

    Parameters
    -----------
    adler_obj : AdlerData
        The AdlerData object to generate expected data from.

    Returns
    -----------
    pd.DataFrame
        DataFrame with expected columns and data from the AdlerData object.
    """
    row_dict = {
        "ssObjectId": adler_obj.ssObjectId,
        "modelId": adler_obj.modelId,
        "updatedMJD": adler_obj.updatedMJD,
    }

    # Add filter-dependent columns
    for f, filter_name in enumerate(adler_obj.filter_list):
        fda = adler_obj.filter_dependent_values[f]
        row_dict[f"{filter_name}_phaseAngle_min"] = fda.phaseAngle_min
        row_dict[f"{filter_name}_phaseAngle_range"] = fda.phaseAngle_range
        row_dict[f"{filter_name}_observationTime_max"] = fda.observationTime_max
        row_dict[f"{filter_name}_arc"] = fda.arc
        row_dict[f"{filter_name}_nobs"] = fda.nobs
        row_dict[f"{filter_name}_n_outliers"] = fda.n_outliers
        row_dict[f"{filter_name}_n_std_outliers"] = fda.n_std_outliers
        row_dict[f"{filter_name}_sustained_outliers"] = fda.sustained_outliers

        # Add model-dependent columns if model exists
        if fda.model_dependent_values is not None:
            model_name = fda.model_name

            # Check if it's a phase model or avg mag model
            if model_name in VALID_PHASE_MODELS:
                # Phase model dependent
                row_dict[f"{filter_name}_{model_name}_H"] = fda.model_dependent_values.H
                row_dict[f"{filter_name}_{model_name}_H_err"] = fda.model_dependent_values.H_err
                row_dict[f"{filter_name}_{model_name}_phase_parameter_1"] = (
                    fda.model_dependent_values.phase_parameter_1
                )
                row_dict[f"{filter_name}_{model_name}_phase_parameter_1_err"] = (
                    fda.model_dependent_values.phase_parameter_1_err
                )
                row_dict[f"{filter_name}_{model_name}_phase_parameter_2"] = (
                    fda.model_dependent_values.phase_parameter_2
                )
                row_dict[f"{filter_name}_{model_name}_phase_parameter_2_err"] = (
                    fda.model_dependent_values.phase_parameter_2_err
                )
            elif model_name in VALID_AVG_MAG_MODELS:
                # Avg mag model dependent
                row_dict[f"{filter_name}_{model_name}_avg_mag"] = fda.model_dependent_values.avg_mag
                row_dict[f"{filter_name}_{model_name}_std_mag"] = fda.model_dependent_values.std_mag

    return pd.DataFrame([row_dict])


# Re-initialize the AdlerData objects to ensure they're correctly populated from the values in this script
mpc_test_obj = AdlerData(ssoid_mpc, filter_list_mpc)
mpc_test_obj.set_modelId(
    model_name=mpc_model_name_in,
    end_mjd=mpc_process_mjd_in,
    data_timespan=mpc_data_timespan_in,
    n_new_nights=mpc_n_new_nights_in,
)
mpc_test_obj.populate_avg_mag_parameters("g", **_g_model_mpc)
mpc_test_obj.populate_avg_mag_parameters("r", **_r_model_mpc)
mpc_test_obj.populate_source_flags("g", mpc_expected_modelId, g_df_mpc)
mpc_test_obj.populate_source_flags("r", mpc_expected_modelId, r_df_mpc)

dp03_test_obj = AdlerData(ssoid_dp03, filter_list_dp03)
dp03_test_obj.set_modelId(
    model_name=dp03_model_name_in,
    end_mjd=dp03_process_mjd_in,
    data_timespan=dp03_data_timespan_in,
    n_new_nights=dp03_n_new_nights_in,
)
dp03_test_obj.populate_phase_parameters("i", **_i_model_dp03)
dp03_test_obj.populate_phase_parameters("r", **_r_model_dp03)
dp03_test_obj.populate_source_flags("i", dp03_expected_modelId, i_df_dp03)
dp03_test_obj.populate_source_flags("r", dp03_expected_modelId, r_df_dp03)


def test_write_to_database_avg_mag(tmp_path):
    db_location = os.path.join(tmp_path, "test_AdlerData_database.db")

    # Write the test object to database
    mpc_test_obj.write_to_database(db_location, write_model_data=True)

    # Read back the written data
    con = sqlite3.connect(db_location)
    written_AdlerData = pd.read_sql_query("SELECT * from AdlerData", con)
    written_FilterDependentAdler = pd.read_sql_query("SELECT * from FilterDependentAdler", con)
    written_AvgMagModelDependentAdler = pd.read_sql_query("SELECT * from AvgMagModelDependentAdler", con)
    con.close()

    # Generate expected data
    expected_data = generate_expected_adlerdata_csv(mpc_test_obj)

    # We don't expect the timestamp column to be the same
    drop_cols = "updatedMJD"
    written_AdlerData = written_AdlerData.drop(columns=drop_cols)
    written_FilterDependentAdler = written_FilterDependentAdler.drop(columns=drop_cols)
    written_AvgMagModelDependentAdler = written_AvgMagModelDependentAdler.drop(columns=drop_cols)

    # Handle dtype and NaN/None conversions
    written_AdlerData = written_AdlerData.where(pd.notna(written_AdlerData), np.nan)
    written_FilterDependentAdler = written_FilterDependentAdler.where(
        pd.notna(written_FilterDependentAdler), np.nan
    )
    written_AvgMagModelDependentAdler = written_AvgMagModelDependentAdler.where(
        pd.notna(written_AvgMagModelDependentAdler), np.nan
    )

    # Compare all columns
    for col in written_AdlerData.columns:
        if col in ["ssObjectId", "modelId"]:
            # String comparison for modelId
            assert expected_data.iloc[0][col] == written_AdlerData.iloc[0][col]
        else:
            # Numeric columns - use assert_almost_equal to handle floating point precision
            assert_almost_equal(expected_data.iloc[0][col], written_AdlerData.iloc[0][col])

    # Compare all columns
    for col in written_FilterDependentAdler.columns:
        if col in ["ssObjectId", "modelId"]:
            # String comparison for modelId
            assert expected_data.iloc[0][col] == written_FilterDependentAdler.iloc[0][col]
        else:
            # Numeric columns - use assert_almost_equal to handle floating point precision
            assert_almost_equal(expected_data.iloc[0][col], written_FilterDependentAdler.iloc[0][col])

    # Compare all columns
    for col in written_AvgMagModelDependentAdler.columns:
        if col in ["ssObjectId", "modelId"]:
            # String comparison for modelId
            assert expected_data.iloc[0][col] == written_AvgMagModelDependentAdler.iloc[0][col]
        else:
            # Numeric columns - use assert_almost_equal to handle floating point precision
            assert_almost_equal(expected_data.iloc[0][col], written_AvgMagModelDependentAdler.iloc[0][col])


def test_write_to_database_phase(tmp_path):
    db_location = os.path.join(tmp_path, "test_AdlerData_database.db")

    # Write the test object to database
    dp03_test_obj.write_to_database(db_location, write_model_data=True)

    # Read back the written data
    con = sqlite3.connect(db_location)
    written_AdlerData = pd.read_sql_query("SELECT * from AdlerData", con)
    written_FilterDependentAdler = pd.read_sql_query("SELECT * from FilterDependentAdler", con)
    written_PhaseModelDependentAdler = pd.read_sql_query("SELECT * from PhaseModelDependentAdler", con)
    con.close()

    # Generate expected data
    expected_data = generate_expected_adlerdata_csv(dp03_test_obj)

    # We don't expect the timestamp column to be the same
    drop_cols = "updatedMJD"
    written_AdlerData = written_AdlerData.drop(columns=drop_cols)
    written_FilterDependentAdler = written_FilterDependentAdler.drop(columns=drop_cols)
    written_PhaseModelDependentAdler = written_PhaseModelDependentAdler.drop(columns=drop_cols)

    # Handle dtype and NaN/None conversions
    written_AdlerData = written_AdlerData.where(pd.notna(written_AdlerData), np.nan)
    written_FilterDependentAdler = written_FilterDependentAdler.where(
        pd.notna(written_FilterDependentAdler), np.nan
    )
    written_PhaseModelDependentAdler = written_PhaseModelDependentAdler.where(
        pd.notna(written_PhaseModelDependentAdler), np.nan
    )

    # Compare all columns
    for col in written_AdlerData.columns:
        if col in ["ssObjectId", "modelId"]:
            # String comparison for modelId
            assert expected_data.iloc[0][col] == written_AdlerData.iloc[0][col]
        else:
            # Numeric columns - use assert_almost_equal to handle floating point precision
            assert_almost_equal(expected_data.iloc[0][col], written_AdlerData.iloc[0][col])

    # Compare all columns
    for col in written_FilterDependentAdler.columns:
        if col in ["ssObjectId", "modelId"]:
            # String comparison for modelId
            assert expected_data.iloc[0][col] == written_FilterDependentAdler.iloc[0][col]
        else:
            # Numeric columns - use assert_almost_equal to handle floating point precision
            assert_almost_equal(expected_data.iloc[0][col], written_FilterDependentAdler.iloc[0][col])

    # Compare all columns
    for col in written_PhaseModelDependentAdler.columns:
        if col in ["ssObjectId", "modelId"]:
            # String comparison for modelId
            assert expected_data.iloc[0][col] == written_PhaseModelDependentAdler.iloc[0][col]
        else:
            # Numeric columns - use assert_almost_equal to handle floating point precision
            assert_almost_equal(expected_data.iloc[0][col], written_PhaseModelDependentAdler.iloc[0][col])


def test_write_flags_to_database(tmp_path):
    db_location = os.path.join(tmp_path, "test_AdlerSourceFlags_database.db")

    # Write the test object to database
    mpc_test_obj.filter_dependent_values[0].source_flags.write_flags_to_database(db_location)
    mpc_test_obj.filter_dependent_values[1].source_flags.write_flags_to_database(db_location)

    con = sqlite3.connect(db_location)
    for filter_name, _df in zip(filter_list_mpc, [g_df_mpc, r_df_mpc]):
        written_flags_data = pd.read_sql_query(
            f"SELECT * FROM AdlerSourceFlags WHERE filter_name='{filter_name}'", con
        )
        for key in _df:
            expected = _df[key].to_numpy()
            written = written_flags_data[key].to_numpy()

            if expected.dtype == object:
                assert_array_equal(expected, written)

            # Numeric column
            else:
                assert_array_almost_equal(expected, written)


def test_populate_from_database():
    mpc_pop_obj = AdlerData(ssoid_mpc, filter_list=filter_list_mpc)
    mpc_pop_obj.populate_from_database(db_path_mpc)

    # Compare basic attributes
    assert mpc_pop_obj.ssObjectId == mpc_test_obj.ssObjectId
    assert mpc_pop_obj.filter_list == mpc_test_obj.filter_list
    assert mpc_pop_obj.modelId == mpc_test_obj.modelId
    # updatedMJD will be different, so we skip it

    # Compare filter-dependent values for each filter
    for f, filter_name in enumerate(mpc_pop_obj.filter_list):
        pop_fda = mpc_pop_obj.filter_dependent_values[f]
        test_fda = mpc_test_obj.filter_dependent_values[f]

        # Filter-dependent attributes
        assert pop_fda.filter_name == test_fda.filter_name
        assert_almost_equal(pop_fda.phaseAngle_min, test_fda.phaseAngle_min)
        assert_almost_equal(pop_fda.phaseAngle_range, test_fda.phaseAngle_range)
        assert_almost_equal(pop_fda.observationTime_max, test_fda.observationTime_max)
        assert_almost_equal(pop_fda.arc, test_fda.arc)
        assert_equal(pop_fda.nobs, test_fda.nobs)
        assert_equal(pop_fda.n_outliers, test_fda.n_outliers)
        assert_equal(pop_fda.n_std_outliers, test_fda.n_std_outliers)
        assert_almost_equal(pop_fda.sustained_outliers, test_fda.sustained_outliers, decimal=5)

        # Model-dependent attributes
        assert pop_fda.model_name == test_fda.model_name

        if pop_fda.model_dependent_values is not None and test_fda.model_dependent_values is not None:
            pop_mdv = pop_fda.model_dependent_values
            test_mdv = test_fda.model_dependent_values

            # Check if it's an avg mag model or phase model
            if pop_fda.model_name in VALID_AVG_MAG_MODELS:
                assert_almost_equal(pop_mdv.avg_mag, test_mdv.avg_mag)
                assert_almost_equal(pop_mdv.std_mag, test_mdv.std_mag)

        # Compare source flags if they exist
        if pop_fda.source_flags is not None and test_fda.source_flags is not None:
            pop_sf = pop_fda.source_flags
            test_sf = test_fda.source_flags

            assert pop_sf.ssObjectId == test_sf.ssObjectId
            assert pop_sf.filter_name == test_sf.filter_name
            assert pop_sf.modelId == test_sf.modelId
            assert pop_sf.n_outliers == test_sf.n_outliers
            assert pop_sf.n_std_outliers == test_sf.n_std_outliers
            assert_array_equal(pop_sf.diaSourceId, test_sf.diaSourceId)
            assert_array_almost_equal(pop_sf.midpointMjdTai, test_sf.midpointMjdTai)
            assert_array_almost_equal(pop_sf.mag_diff, test_sf.mag_diff)
            assert_array_almost_equal(pop_sf.std_diff, test_sf.std_diff)

    # Test with dp03 data (phase model)
    dp03_pop_obj = AdlerData(ssoid_dp03, filter_list=filter_list_dp03)
    dp03_pop_obj.populate_from_database(db_path_dp03)

    # Compare basic attributes
    assert dp03_pop_obj.ssObjectId == dp03_test_obj.ssObjectId
    assert dp03_pop_obj.filter_list == dp03_test_obj.filter_list
    assert dp03_pop_obj.modelId == dp03_test_obj.modelId

    # Compare filter-dependent values for each filter
    for f, filter_name in enumerate(dp03_pop_obj.filter_list):
        pop_fda = dp03_pop_obj.filter_dependent_values[f]
        test_fda = dp03_test_obj.filter_dependent_values[f]

        # Filter-dependent attributes
        assert pop_fda.filter_name == test_fda.filter_name
        assert_almost_equal(pop_fda.phaseAngle_min, test_fda.phaseAngle_min)
        assert_almost_equal(pop_fda.phaseAngle_range, test_fda.phaseAngle_range)
        assert_almost_equal(pop_fda.observationTime_max, test_fda.observationTime_max)
        assert_almost_equal(pop_fda.arc, test_fda.arc)
        assert_equal(pop_fda.nobs, test_fda.nobs)
        assert_equal(pop_fda.n_outliers, test_fda.n_outliers)
        assert_equal(pop_fda.n_std_outliers, test_fda.n_std_outliers)
        assert_almost_equal(pop_fda.sustained_outliers, test_fda.sustained_outliers, decimal=5)

        # Model-dependent attributes
        assert pop_fda.model_name == test_fda.model_name

        if pop_fda.model_dependent_values is not None and test_fda.model_dependent_values is not None:
            pop_mdv = pop_fda.model_dependent_values
            test_mdv = test_fda.model_dependent_values

            # Phase model comparison
            assert_almost_equal(pop_mdv.H, test_mdv.H)
            assert_almost_equal(pop_mdv.H_err, test_mdv.H_err)
            assert_almost_equal(pop_mdv.phase_parameter_1, test_mdv.phase_parameter_1)
            assert_almost_equal(pop_mdv.phase_parameter_1_err, test_mdv.phase_parameter_1_err)
            # phase_parameter_2 may be NaN
            if np.isnan(test_mdv.phase_parameter_2):
                assert np.isnan(pop_mdv.phase_parameter_2)
            else:
                assert_almost_equal(pop_mdv.phase_parameter_2, test_mdv.phase_parameter_2)
            if np.isnan(test_mdv.phase_parameter_2_err):
                assert np.isnan(pop_mdv.phase_parameter_2_err)
            else:
                assert_almost_equal(pop_mdv.phase_parameter_2_err, test_mdv.phase_parameter_2_err)

        # Compare source flags if they exist
        if pop_fda.source_flags is not None and test_fda.source_flags is not None:
            pop_sf = pop_fda.source_flags
            test_sf = test_fda.source_flags

            assert pop_sf.ssObjectId == test_sf.ssObjectId
            assert pop_sf.filter_name == test_sf.filter_name
            assert pop_sf.modelId == test_sf.modelId
            assert pop_sf.n_outliers == test_sf.n_outliers
            assert pop_sf.n_std_outliers == test_sf.n_std_outliers
            assert_array_equal(pop_sf.diaSourceId, test_sf.diaSourceId)
            assert_array_almost_equal(pop_sf.midpointMjdTai, test_sf.midpointMjdTai)
            assert_array_almost_equal(pop_sf.mag_diff, test_sf.mag_diff)
            assert_array_almost_equal(pop_sf.std_diff, test_sf.std_diff)
