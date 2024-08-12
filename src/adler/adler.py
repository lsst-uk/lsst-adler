import logging
import argparse
import astropy.units as u
import numpy as np
import pandas as pd
import os

from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid
from adler.science.PhaseCurve import PhaseCurve
from adler.utilities.AdlerCLIArguments import AdlerCLIArguments
from adler.utilities.adler_logging import setup_adler_logging
from adler.utilities.readin_utilities import read_in_SSObjectID_file
from adler.utilities.plotting_utilities import plot_errorbar
import adler.utilities.science_utilities as sci_utils

logger = logging.getLogger(__name__)


def runAdler(cli_args):
    logger.info("Beginning Adler.")

    N_pc_fit = 10  # minimum number of data points to fit phase curve
    diff_cut = 1.0  # magnitude difference used to identify outliers
    obs_cols = ["diaSourceId", "midPointMjdTai", "outlier"]  # observation columns to use

    # Define colour parameters
    # set number of reference observations to use for colour estimate
    N_ref = 5

    # observation and filter field names
    x_col = "midPointMjdTai"
    y_col = "AbsMag"
    yerr_col = "magErr"

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

        logger.info("Data successfully ingested.")

        # now let's do some phase curves!
        logger.info("Calculating phase curves...")

        # operate on each filter in turn
        for filt in cli_args.filter_list:
            logger.info("fit {} filter data".format(filt))

            # get the filter SSObject metadata
            sso = planetoid.SSObject_in_filter(filt)

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
            mask = df_obs["midPointMjdTai"] < t0

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
                print("insufficient data, continue")
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
                cli_args.outpath, cli_args.ssObjectId, int(t0), filt
            )
            # TODO: make the plots folder if it does not already exist?
            print("Save figure: {}".format(fig_file))
            logger.info("Save figure: {}".format(fig_file))
            fig = plot_errorbar(planetoid, fig=fig, filename=fig_file)

        # analyse colours for the filters provided
        logger.info("Calculate colours: {}".format(cli_args.colour_list))

        # cycle through the filters, calculating a colour relative to the next filter
        # note that the order in which cli_args.filter_list is passed will determine which colours are calculated
        for colour in cli_args.colour_list:
            col_filts = colour.split("-")
            filt_obs = col_filts[0]
            filt_ref = col_filts[1]
            logger.info("Determine {} - {} colour".format(filt_obs, filt_ref))

            # define colour field names
            colour = "{}-{}".format(filt_obs, filt_ref)
            colErr = "{}-{}Err".format(filt_obs, filt_ref)
            delta_t_col = "delta_t_{}".format(colour)
            y_ref_col = "{}_{}".format(y_col, filt_ref)
            x1_ref_col = "{}1_{}".format(x_col, filt_ref)
            x2_ref_col = "{}2_{}".format(x_col, filt_ref)


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
        default=["g-r", "r-i"],
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
        help="Output path location. Default is current working directory.",
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

    runAdler(cli_args)


if __name__ == "__main__":
    main()
