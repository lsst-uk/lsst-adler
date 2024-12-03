import logging
import argparse
import astropy.units as u
from astropy.time import Time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sqlite3

from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.objectdata.AdlerData import AdlerData
from adler.science.PhaseCurve import PhaseCurve, ADLER_SBPY_DICT
from adler.science.Colour import col_obs_ref
from adler.utilities.AdlerCLIArguments import AdlerCLIArguments
from adler.utilities.adler_logging import setup_adler_logging
from adler.utilities.readin_utilities import read_in_SSObjectID_file
from adler.utilities.plotting_utilities import plot_errorbar
import adler.utilities.science_utilities as sci_utils

"""
Copy of adler_run.oy which can be used to make the AdlerData testing database
"""

ssObjectId = "8268570668335894776"
filter_list = ["g", "r", "i"]
colour_list = None
date_range = [60000.0, 67300.0]
outpath = "."
db_name = "test_AdlerData_database.db"
sql_filename = "testing_database.db"
plot_show = False
# phase_model="HG12_Pen16"
phase_model_list = ["HG12_Pen16", "HG"]

# adler parameters
N_pc_fit = 10  # minimum number of data points to fit phase curve
diff_cut = 1.0  # magnitude difference used to identify outliers
obs_cols = ["diaSourceId", "midPointMjdTai", "outlier"]  # observation columns to use


# Define colour parameters
# set number of reference observations to use for colour estimate
N_ref = 5

ssObjectId_list = [ssObjectId]

# consider each ssObjectId in the list separately
for i, ssObjectId in enumerate(ssObjectId_list):
    print(
        "Query data in the range: {} <= date <= {}".format(date_range[0], date_range[1])
    )  # the adler planetoid date_range is used in an SQL BETWEEN statement which is inclusive

    # load ssObjectId data
    if sql_filename:  # load from a local database, if provided
        msg = "query sql database {}".format(sql_filename)
        print(msg)
        planetoid = AdlerPlanetoid.construct_from_SQL(
            ssObjectId,
            filter_list=filter_list,
            date_range=date_range,
            sql_filename=sql_filename,
        )
    else:  # otherwise load from the Rubin Science Platform
        msg = "query RSP"
        print(msg)
        planetoid = AdlerPlanetoid.construct_from_RSP(ssObjectId, filter_list, date_range)

    # TODO: Here we would load the AdlerData object from our data tables
    adler_data = AdlerData(ssObjectId, planetoid.filter_list)
    print(adler_data.__dict__)

    # now let's do some phase curves!

    # make an empty figure
    fig = plot_errorbar(planetoid, filt_list=[])

    # do each model
    for phase_model in phase_model_list:

        print(phase_model)
        # get the name of the phase parameter
        # TODO: make an option for two parameter HG1G2
        phase_parameter_1 = ADLER_SBPY_DICT[phase_model]["phase_parameter_1"]
        print(phase_model, phase_parameter_1)

        # set default phase parameter values
        # TODO: make an option for two parameter HG1G2
        if phase_model == "HG":
            phase_param_1_default = 0.15
        else:
            phase_param_1_default = 0.62

        # operate on each filter in turn
        for filt in filter_list:

            # get the filter SSObject metadata
            try:
                sso = planetoid.SSObject_in_filter(filt)
            except:
                continue

            # get the observations
            obs = planetoid.observations_in_filter(filt)
            df_obs = pd.DataFrame(obs.__dict__)
            df_obs["outlier"] = [False] * len(df_obs)

            # Determine the reference phase curve model
            # TODO: We would load the best phase curve model available in AdlerData here

            # initial simple phase curve filter model with fixed phase_parameter
            # use the ssObject value for H as initial guess, this is in HG12_Pen16
            # TODO: use the ssObject value for phase parameter as initial guess?
            pc = PhaseCurve(
                # H=sso.H * u.mag,
                H=sso.H,
                phase_parameter_1=phase_param_1_default,
                model_name=phase_model,
            )

            # only fit phase_parameter when sufficient data is available
            if len(df_obs) < N_pc_fit:
                msg = "Do not fit {}, use {}={:.2f}".format(
                    phase_parameter_1, phase_parameter_1, pc.phase_parameter_1
                )
                # pc.model_function.G12.fixed = True
                getattr(pc.model_function, phase_parameter_1).fixed = True
            else:
                # pc.model_function.G12.fixed = False
                getattr(pc.model_function, phase_parameter_1).fixed = False

            # do a phase model fit to the past data
            pc_fit = pc.FitModel(
                # np.array(df_obs["phaseAngle"]) * u.deg,
                # np.array(df_obs["reduced_mag"]) * u.mag,
                # np.array(df_obs["magErr"]) * u.mag,
                np.radians(np.array(df_obs["phaseAngle"])),
                np.array(df_obs["reduced_mag"]),
                np.array(df_obs["magErr"]),
            )
            pc_fit = pc.InitModelSbpy(pc_fit)

            # Store the fitted values, and metadata, in an AdlerData object
            ad_params = pc_fit.__dict__
            ad_params["phaseAngle_min"] = np.amin(df_obs["phaseAngle"])  # * u.deg
            ad_params["phaseAngle_range"] = np.ptp(df_obs["phaseAngle"])  # * u.deg
            ad_params["arc"] = np.ptp(df_obs["midPointMjdTai"])  # * u.d
            ad_params["nobs"] = len(df_obs)
            ad_params["modelFitMjd"] = Time.now().mjd
            # adler_data.populate_phase_parameters(filt, **pc_fit.__dict__)
            # TODO: replace any None with np.nan? e.g. phase_parameter_2?
            adler_data.populate_phase_parameters(filt, **ad_params)

            # print the test data set
            print(adler_data.get_phase_parameters_in_filter(filt, phase_model).__dict__)

            # add to plot
            ax1 = fig.axes[0]
            # TODO: set colours based on filter
            ax1.scatter(df_obs["phaseAngle"], df_obs["reduced_mag"])
            alpha = np.linspace(0, np.amax(obs.phaseAngle)) * u.deg
            ax1.plot(
                alpha.value,
                # pc_fit.ReducedMag(alpha).value,
                # label="{}, H={:.2f}, G12={:.2f}".format(filt, pc_fit.H.value, pc_fit.phase_parameter_1),
                pc_fit.ReducedMag(alpha),
                label="{}, H={:.2f}, {}={:.2f}".format(
                    filt, pc_fit.H, phase_parameter_1, pc_fit.phase_parameter_1
                ),
            )
        ax1.legend()

        # TODO: Use a CLI arg flag to open figure interactively instead of saving?
        if plot_show:
            plt.show()
        # Save figures at the outpath location
        else:
            fig_file = "{}/phase_curve_{}_{}_{}.png".format(
                outpath, ssObjectId, phase_model, int(np.amax(df_obs["midPointMjdTai"]))
            )
            msg = "Save figure: {}".format(fig_file)
            print(msg)
            fig = plot_errorbar(planetoid, fig=fig, filename=fig_file)  # TODO: add titles with filter name?
            plt.close()

        # analyse colours for the filters provided

        # if requested cycle through the filters, calculating a colour relative to the next filter
        if not colour_list:
            colour_list = []
        else:
            colour_list = colour_list

        # note that the order in which cli_args.filter_list is passed will determine which colours are calculated
        for colour in colour_list:
            col_filts = colour.split("-")
            filt_obs = col_filts[0]
            filt_ref = col_filts[1]

            if ~np.isin([filt_obs, filt_ref], planetoid.filter_list).all():
                missing_filts = np.array([filt_obs, filt_ref])[
                    ~np.isin([filt_obs, filt_ref], planetoid.filter_list)
                ]
                continue

            # determine the filt_obs - filt_ref colour
            # generate a plot
            if plot_show:
                plot_dir = None
            else:
                plot_dir = outpath

            col_dict = col_obs_ref(
                planetoid,
                adler_data,
                phase_model=phase_model,
                filt_obs=filt_obs,
                filt_ref=filt_ref,
                N_ref=N_ref,
                # x1 = x1,
                # x2 = x2,
                plot_dir=plot_dir,
                plot_show=plot_show,
            )

            # print(col_dict)

    # Output adler values to a database if a db_name is provided
    print(adler_data.__dict__)
    if db_name:
        adler_db = "{}/{}".format(outpath, db_name)

        if os.path.isfile(adler_db):
            print("clear old db file")
            os.remove(adler_db)

        msg = "write to {}".format(adler_db)
        print(msg)
        adler_data.write_row_to_database(adler_db)

# create the test csv file
fname_csv = "{}/test_SQL_database_table.csv".format(outpath)
conn = sqlite3.Connection(adler_db)
df = pd.read_sql("SELECT * from AdlerData", conn)
print("write to {}".format(fname_csv))
df.to_csv(fname_csv)
