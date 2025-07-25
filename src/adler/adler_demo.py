import logging
import argparse
import astropy.units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.objectdata.AdlerData import AdlerData
from adler.science.PhaseCurve import PhaseCurve
from adler.science.Colour import col_obs_ref
from adler.utilities.AdlerCLIArguments import AdlerCLIArguments
from adler.utilities.adler_logging import setup_adler_logging
from adler.utilities.readin_utilities import read_in_SSObjectID_file
from adler.utilities.plotting_utilities import plot_errorbar
import adler.utilities.science_utilities as sci_utils

logger = logging.getLogger(__name__)


def runAdlerDemo(cli_args):
    """
    This function is an experimental version of adlerRun (adler_run.py).
    Here the observations are split into "past" and "present" and the phase curve model is fit to past observations.
    The "present" observations are tested to see if they are outlying the phase curve model.
    Colours and plots are also made.
    """
    logger.info("Beginning Adler.")

    # adler parameters
    N_pc_fit = 10  # minimum number of data points to fit phase curve
    diff_cut = 1.0  # magnitude difference used to identify outliers
    obs_cols = ["diaSourceId", "midPointMjdTai", "outlier"]  # observation columns to use
    phase_model = "HG12_Pen16"  # which phase curve model to fit

    # Define colour parameters
    # set number of reference observations to use for colour estimate
    N_ref = 5

    if cli_args.ssObjectId_list:
        ssObjectId_list = read_in_SSObjectID_file(cli_args.ssObjectId_list)
    else:
        ssObjectId_list = [cli_args.ssObjectId]

    # consider each ssObjectId in the list separately
    for i, ssObjectId in enumerate(ssObjectId_list):
        logger.info("Processing object {}/{}.".format(i + 1, len(ssObjectId_list)))
        logger.info("Ingesting all data for object {} from RSP...".format(cli_args.ssObjectId))
        logger.info(
            "Query data in the range: {} <= date <= {}".format(cli_args.date_range[0], cli_args.date_range[1])
        )  # the adler planetoid date_range is used in an SQL BETWEEN statement which is inclusive
        print(
            "Query data in the range: {} <= date <= {}".format(cli_args.date_range[0], cli_args.date_range[1])
        )  # the adler planetoid date_range is used in an SQL BETWEEN statement which is inclusive
        logger.info("Consider the filters: {}".format(cli_args.filter_list))

        # load ssObjectId data
        if cli_args.sql_filename:  # load from a local database, if provided
            msg = "query sql database {}".format(cli_args.sql_filename)
            logger.info(msg)
            print(msg)
            planetoid = AdlerPlanetoid.construct_from_SQL(
                ssObjectId,
                filter_list=cli_args.filter_list,
                date_range=cli_args.date_range,
                sql_filename=cli_args.sql_filename,
            )
        else:  # otherwise load from the Rubin Science Platform
            msg = "query RSP"
            logger.info(msg)
            print(msg)
            planetoid = AdlerPlanetoid.construct_from_RSP(
                ssObjectId, cli_args.filter_list, cli_args.date_range
            )

        # TODO: Here we would load the AdlerData object from our data tables
        adler_data = AdlerData(ssObjectId, planetoid.filter_list)
        print(adler_data.__dict__)

        logger.info("Data successfully ingested.")

        # now let's do some phase curves!
        logger.info("Calculating phase curves...")

        # operate on each filter in turn
        for filt in cli_args.filter_list:
            logger.info("fit {} filter data".format(filt))

            # get the filter SSObject metadata
            try:
                sso = planetoid.SSObject_in_filter(filt)
            except:
                logger.info("error loading SSObject in filter {}".format(filt))
                continue

            # get the observations
            obs = planetoid.observations_in_filter(filt)
            df_obs = pd.DataFrame(obs.__dict__)
            df_obs["outlier"] = [False] * len(df_obs)
            logger.info("{} observations retrieved".format(len(df_obs)))

            # load and merge the previous obs
            # TODO: replace this part with classifications loaded from adlerData
            save_file = "{}/df_outlier_{}_{}.csv".format(cli_args.outpath, cli_args.ssObjectId, filt)
            if os.path.isfile(save_file):
                logger.info("load previously classified observations: {}".format(save_file))
                _df_obs = pd.read_csv(save_file, index_col=0)
                df_obs = df_obs.merge(_df_obs, on=["diaSourceId", "midPointMjdTai"], how="left")
                df_obs.loc[pd.isnull(df_obs["outlier_y"]), "outlier_y"] = (
                    False  # ensure that classifications exist (nan entries can only be false?). Weird behaviour here for g filter, is it to do with when new g obs appear relative to r/i etc?
                )
                df_obs = df_obs.rename({"outlier_y": "outlier"}, axis=1)
                df_obs = df_obs.drop("outlier_x", axis=1)
            else:
                logger.info("no previously classified observations to load")

            # define the date range to for new observations taken in the night to be analysed
            logger.info(
                "Most recent {} filter observation in query: date = {}".format(
                    filt, np.amax(df_obs["midPointMjdTai"])
                )
            )
            t1 = int(np.amax(df_obs["midPointMjdTai"])) + 1
            t0 = t1 - 1

            # get all past observations
            # mask = df_obs["midPointMjdTai"] < t0
            mask = (df_obs["midPointMjdTai"] < t0) & (df_obs["outlier"] == False)  # reject any past outliers

            # split observations into "old" and "new"
            df_obs_old = df_obs[(mask)]
            df_obs_new = df_obs[~mask]
            logger.info("Previous observations (date < {}): {}".format(t0, len(df_obs_old)))
            logger.info("New observations ({} <= date < {}): {}".format(t0, t1, len(df_obs_new)))

            # Determine the reference phase curve model
            # TODO: We would load the best phase curve model available in AdlerData here

            # we need sufficient past observations to fit the phase curve model
            if len(df_obs_old) < 2:
                print("save {}".format(save_file))
                df_save = df_obs[obs_cols]
                df_save.to_csv(save_file)
                print("insufficient data, use default SSObject phase model and continue")
                logger.info("insufficient data, use default SSObject phase model and continue")

                # use the default SSObject phase parameter if there is no better information
                pc_dict = {
                    "H": sso.H * u.mag,
                    "H_err": sso.Herr * u.mag,
                    "phase_parameter_1": sso.G12,
                    "phase_parameter_1_err": sso.G12err,
                    "model_name": "HG12_Pen16",
                }
                adler_data.populate_phase_parameters(filt, **pc_dict)
                continue

            # initial simple phase curve filter model with fixed G12
            pc = PhaseCurve(
                H=sso.H * u.mag,
                phase_parameter_1=0.62,
                model_name="HG12_Pen16",
            )

            # only fit G12 when sufficient data is available
            if len(df_obs_old) < N_pc_fit:
                pc.model_function.G12.fixed = True
            else:
                pc.model_function.G12.fixed = False

            # do a HG12_Pen16 fit to the past data
            pc_fit = pc.FitModel(
                np.array(df_obs_old["phaseAngle"]) * u.deg,
                np.array(df_obs_old["reduced_mag"]) * u.mag,
                np.array(df_obs_old["magErr"]) * u.mag,
            )
            pc_fit = pc.InitModelSbpy(pc_fit)

            # TODO: Here the best fit should be pushed back to our AdlerData tables
            adler_data.populate_phase_parameters(filt, **pc_fit.__dict__)

            # find outliers in new data
            # calculate data - model residuals
            res = (np.array(df_obs_new["reduced_mag"]) * u.mag) - pc_fit.ReducedMag(
                np.array(df_obs_new["phaseAngle"]) * u.deg
            )
            outlier_flag = sci_utils.outlier_diff(res.value, diff_cut=diff_cut)
            df_obs.loc[~mask, "outlier"] = outlier_flag

            # save the df_obs subset with outlier classification
            df_save = df_obs[obs_cols]
            print("save classifications: {}".format(save_file))
            logger.info("save classifications: {}".format(save_file))
            df_save.to_csv(save_file)

            # make a plot
            fig = plot_errorbar(planetoid, filt_list=[])
            ax1 = fig.axes[0]
            ax1.scatter(df_obs_old["phaseAngle"], df_obs_old["reduced_mag"], c="C0")
            alpha = np.linspace(0, np.amax(obs.phaseAngle)) * u.deg
            ax1.plot(alpha.value, pc_fit.ReducedMag(alpha).value, label="t={}".format(int(t0)))
            ax1.scatter(
                df_obs_new["phaseAngle"], df_obs_new["reduced_mag"], edgecolor="r", facecolor="none", zorder=3
            )
            out_mask = df_obs["outlier"] == True
            ax1.scatter(
                df_obs.loc[out_mask]["phaseAngle"],
                df_obs.loc[out_mask]["reduced_mag"],
                c="r",
                marker="x",
                s=75,
                zorder=3,
            )
            fig_file = "{}/plots/phase_curve_{}_{}_{}.png".format(
                cli_args.outpath, cli_args.ssObjectId, filt, int(t0)
            )
            # TODO: make the plots folder if it does not already exist?
            print("Save figure: {}".format(fig_file))
            logger.info("Save figure: {}".format(fig_file))
            fig = plot_errorbar(planetoid, fig=fig, filename=fig_file)  # TODO: add titles with filter name?
            plt.close()

        # analyse colours for the filters provided
        logger.info("Calculate colours: {}".format(cli_args.colour_list))

        # if requested cycle through the filters, calculating a colour relative to the next filter
        if not cli_args.colour_list:
            colour_list = []
        else:
            colour_list = cli_args.colour_list

        # note that the order in which cli_args.filter_list is passed will determine which colours are calculated
        for colour in colour_list:
            col_filts = colour.split("-")
            filt_obs = col_filts[0]
            filt_ref = col_filts[1]

            if ~np.isin([filt_obs, filt_ref], planetoid.filter_list).all():
                missing_filts = np.array([filt_obs, filt_ref])[
                    ~np.isin([filt_obs, filt_ref], planetoid.filter_list)
                ]
                logger.info(
                    "Filter(s) {} are missing for determining {} colour".format(missing_filts, colour)
                )
                continue

            logger.info("Determine {} colour".format(colour))

            # TODO: replace this with a colour loaded from adlerData
            save_file_colour = "{}/df_colour_{}_{}.csv".format(cli_args.outpath, cli_args.ssObjectId, colour)
            if os.path.isfile(save_file_colour):
                print("load previous colours from file: {}".format(save_file_colour))
                df_col = pd.read_csv(save_file_colour, index_col=0)
                # Check the last colour calculation date (x_obs) to avoid recalculation
                obs = planetoid.observations_in_filter(filt_obs)
                df_obs = pd.DataFrame(obs.__dict__)
                if np.amax(df_col["midPointMjdTai"]) >= np.amax(df_obs["midPointMjdTai"]):
                    print("colour already calculated, skip")
                    continue

            else:
                df_col = pd.DataFrame()

            # determine the filt_obs - filt_ref colour
            # generate a plot
            col_dict = col_obs_ref(
                planetoid,
                adler_data,
                phase_model=phase_model,
                filt_obs=filt_obs,
                filt_ref=filt_ref,
                N_ref=N_ref,
                # x1 = x1,
                # x2 = x2,
                plot_dir="{}/plots".format(cli_args.outpath),
            )

            print(col_dict)

            # save the colour data
            print("Append new colour and save to file: {}".format(save_file_colour))
            df_col = pd.concat([df_col, pd.DataFrame([col_dict])])
            df_col = df_col.reset_index(drop=True)
            df_col.to_csv(save_file_colour)

            # TODO: reject unreliable colours, e.g. high colErr or delta_t_col
            # TODO: determine if colour is outlying
            # compare this new colour to previous colour(s)


def main():
    parser = argparse.ArgumentParser(description="Runs Adler for select planetoid(s) and given user input.")

    # the below group ensures that AT LEAST one of the below arguments is included, but NOT both
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-s", "--ssObjectId", help="SSObject ID of planetoid.", type=str, default=None)
    input_group.add_argument(
        "-sl",
        "--ssObjectId_list",
        help="Filepath to text file listing multiple SSObjectIds.",
        type=str,
        default=None,
    )

    optional_group = parser.add_argument_group("Optional arguments")
    optional_group.add_argument(
        "-f",
        "--filter_list",
        help="Filters to be analysed.",
        nargs="*",
        type=str,
        default=["u", "g", "r", "i", "z", "y"],
    )
    optional_group.add_argument(
        "-c",
        "--colour_list",
        help="Colours to be analysed.",
        nargs="*",
        type=str,
        # default=["g-r", "r-i"],
        default=None,
    )
    optional_group.add_argument(
        "-d",
        "--date_range",
        help="Minimum and maximum MJD(TAI) of required observations. Default is to pull all observations.",
        nargs=2,
        type=float,
        default=[60000.0, 67300.0],
    )
    optional_group.add_argument(
        "-o",
        "--outpath",
        help="Output path location. Default is current working directory.",  # TODO: make adler create the outpath directory on start up if it does not exist? Also the "plots" dir within?
        type=str,
        default="./",
    )
    optional_group.add_argument(
        "-n",
        "--db_name",
        help="Stem filename of output database. If this doesn't exist, it will be created. Default: adler_out.",
        type=str,
        default="adler_out",
    )
    optional_group.add_argument(
        "-i",
        "--sql_filename",
        help="Optional input path location of a sql database file containing observations",
        type=str,
        default=None,
    )

    args = parser.parse_args()

    cli_args = AdlerCLIArguments(args)

    adler_logger = setup_adler_logging(cli_args.outpath)

    cli_args.logger = adler_logger

    runAdlerDemo(cli_args)


if __name__ == "__main__":
    main()
