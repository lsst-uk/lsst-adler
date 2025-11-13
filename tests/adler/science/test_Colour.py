import numpy as np
import pandas as pd
from numpy.testing import assert_almost_equal, assert_array_equal, assert_array_almost_equal
import pytest
import os
import astropy.units as u

from adler.science.Colour import col_obs_ref
from adler.utilities.tests_utilities import get_test_data_filepath
from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.objectdata.AdlerData import AdlerData
from adler.science.PhaseCurve import PhaseCurve
from adler.utilities.science_utilities import get_df_obs_filt

# set up a good test object (main belt asteroid)
ssoid = 6098332225018
test_db_path = get_test_data_filepath("testing_database.db")
planetoid = AdlerPlanetoid.construct_from_SQL(ssoid, sql_filename=test_db_path)
adler_data = AdlerData(ssoid, planetoid.filter_list)

# Set time boundaries for a given apparition (apparition index 3 for this object)
test_app_i = 3
t_app1 = 61469
t_app2 = 61678

# define columns
column_dict = {
    "x_col": "midPointMjdTai",
    "y_col": "AbsMag",
    "y_col": "reduced_mag",
    "yerr_col": "magErr",
    "obsId_col": "diaSourceId",
    "filt_obs": "g",  # observation filter
    "filt_ref": "r",  # reference filter (we are calculating a filt_obs - filt_ref colour)
}

# Load the stored colour results for this apparition from a df csv in tests/data, see also lsst-adler/notebooks/colour_functions_testing.ipynb
# test_data_path = "/Users/jrobinson/lsst-adler/tests/data" # pytest does not seem to like relative file paths, use absolute
test_data_path = "/".join(os.path.abspath(__file__).split("/")[:-3]) + "/data"
df_col_fname = "{}/df_{}_{}_{}_app_{}".format(
    test_data_path, ssoid, column_dict["filt_obs"], column_dict["filt_ref"], test_app_i
)
df_N_ref_1 = pd.read_csv("{}_N_ref_{}.csv".format(df_col_fname, 1), index_col=0)
df_N_ref_3 = pd.read_csv("{}_N_ref_{}.csv".format(df_col_fname, 3), index_col=0)
df_N_ref_5 = pd.read_csv("{}_N_ref_{}.csv".format(df_col_fname, 5), index_col=0)

# fit adler phase curve models
# TODO: replace this with the stored alder_cols when they are implemented in the database
adler_data = AdlerData(ssoid, planetoid.filter_list)
for filt in [column_dict["filt_obs"], column_dict["filt_ref"]]:
    sso = planetoid.SSObject_in_filter(filt)
    obs = planetoid.observations_in_filter(filt)

    H = sso.H
    G12 = sso.G12
    pc = PhaseCurve(
        H=H,
        # H=H * u.mag,
        phase_parameter_1=G12,
        model_name="HG12_Pen16",
    )

    pc_fit = pc.FitModel(
        np.radians(np.array(getattr(obs, "phaseAngle"))),
        np.array(getattr(obs, "reduced_mag")),
        # np.array(getattr(obs, "phaseAngle")) * u.deg,
        # np.array(getattr(obs, "reduced_mag")) * u.mag,
    )
    pc = pc.InitModelSbpy(pc_fit)

    adler_data.populate_phase_parameters(filt, **pc.ReturnModelDict())


def test_col_obs_ref(
    planetoid=planetoid,
    adler_data=adler_data,
    column_dict=column_dict,
    N_ref_list=[1, 3, 5],
    df_ref_list=[df_N_ref_1, df_N_ref_3, df_N_ref_5],
):
    # gather the observations for the apparition
    ad_g = adler_data.get_phase_parameters_in_filter(column_dict["filt_obs"], "HG12_Pen16")
    pc_g = PhaseCurve().InitModelDict(ad_g.__dict__)
    df_obs = get_df_obs_filt(
        planetoid,
        column_dict["filt_obs"],
        x1=t_app1,
        x2=t_app2,
        col_list=[column_dict["obsId_col"], column_dict["y_col"], column_dict["yerr_col"]],
        pc_model=pc_g,
    )

    # test multiple N_ref
    for N_ref, df_ref in zip(N_ref_list, df_ref_list):
        col_dict_list = []

        # simulate stepping through each filt_obs observation
        x1 = t_app1
        for xi in range(len(df_obs)):
            x2 = df_obs.iloc[xi][column_dict["x_col"]]

            # do the colour finding function here
            col_dict = col_obs_ref(
                planetoid,
                adler_data,
                N_ref=N_ref,
                x1=x1,
                x2=x2,
                **column_dict,
            )
            col_dict_list.append(col_dict)

        # store results as a dataframe
        df_col = pd.DataFrame(col_dict_list)
        df_col = df_col.merge(
            df_obs[[x for x in list(df_obs) if x != column_dict["x_col"]]], on=column_dict["obsId_col"]
        )  # merge with observation data (avoiding duplicating x_col)

        # compare results df to stored df
        # NB: see also pytest.approx for comparing dictionaries
        for x in list(df_col):
            # print(x)
            assert_array_almost_equal(np.array(df_col[x]), np.array(df_ref[x]))
            # TODO: diasourceId failing on ubuntu tests, due to float? need to save as str/int?


# test_col_obs_ref(
#     planetoid=planetoid,
#     adler_data=adler_data,
#     column_dict=column_dict,
#     N_ref_list=[1, 3, 5],
#     df_ref_list=[df_N_ref_1, df_N_ref_3, df_N_ref_5],
# )
