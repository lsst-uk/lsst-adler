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
    diff_cut = 1.0

    if cli_args.ssObjectId_list:
        ssObjectId_list = read_in_SSObjectID_file(cli_args.ssObjectId_list)
    else:
        ssObjectId_list = [cli_args.ssObjectId]

    for i, ssObjectId in enumerate(ssObjectId_list):
        logger.info("Processing object {}/{}.".format(i + 1, len(ssObjectId_list)))
        logger.info("Ingesting all data for object {} from RSP...".format(cli_args.ssObjectId))

        # if cli_args.sql_filename != "None":
        if cli_args.sql_filename:
            msg = "query sql database {}".format(cli_args.sql_filename)
            logger.info(msg)
            print(msg)
            planetoid = AdlerPlanetoid.construct_from_SQL(
                ssObjectId,
                filter_list=cli_args.filter_list,
                date_range=cli_args.date_range,
                sql_filename=cli_args.sql_filename,
            )
        else:
            msg = "query RSP"
            logger.info(msg)
            print(msg)
            planetoid = AdlerPlanetoid.construct_from_RSP(
                ssObjectId, cli_args.filter_list, cli_args.date_range
            )

        logger.info("Data successfully ingested.")
        logger.info("Calculating phase curves...")

        # now let's do some phase curves!

        # # operate on each filter in turn
        for filt in planetoid.filter_list:
            print("fit {} filter data".format(filt))

            # get the filter SSObject metadata
            sso = planetoid.SSObject_in_filter(filt)

            # get the observations
            obs = planetoid.observations_in_filter(filt)
            df_obs = pd.DataFrame(obs.__dict__)
            df_obs["outlier"] = [False] * len(df_obs)
            print(len(df_obs))

            # load and merge the previous obs
            save_file = "{}/df_outlier_{}.csv".format(cli_args.outpath, cli_args.ssObjectId)
            if os.path.isfile(save_file):
                print("load {}".format(save_file))
                _df_obs = pd.read_csv(save_file, index_col=0)
                df_obs = df_obs.merge(_df_obs, on="midPointMjdTai", how="left")
                df_obs = df_obs.rename({"outlier_y": "outlier"}, axis=1)
                df_obs = df_obs.drop("outlier_x", axis=1)
            else:
                print("no previous obs to load")

            print(cli_args.date_range)
            print(np.amax(df_obs["midPointMjdTai"]))
            t1 = int(np.amax(df_obs["midPointMjdTai"])) + 1
            t0 = t1 - 1

            # t_mask = (df_obs["midPointMjdTai"]<t1)
            # _df_obs = df_obs[t_mask]
            mask = df_obs["midPointMjdTai"] < t0
            df_obs_old = df_obs[(mask)]
            df_obs_new = df_obs[~mask]
            print(t0, t1, len(df_obs_old), len(df_obs_new))

            if len(df_obs_old) < 2:
                print("save {}".format(save_file))
                df_save = df_obs[["midPointMjdTai", "outlier"]]
                df_save.to_csv(save_file)
                print("insufficient data, continue")
                continue

            # initial simple phase curve filter model with fixed G12
            pc = PhaseCurve(
                abs_mag=sso.H * u.mag,
                phase_param=0.62,
                model_name="HG12_Pen16",
            )

            if len(df_obs_old) < N_pc_fit:
                # use an assumed value of G12 until more data is available
                pc.model_function.G12.fixed = True
            else:
                pc.model_function.G12.fixed = False

            # do a simple HG12_Pen16 fit to the past data
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
            # print(res)
            outlier_flag = sci_utils.outlier_diff(res.value, diff_cut=diff_cut)
            print(outlier_flag)
            df_obs.loc[~mask, "outlier"] = outlier_flag

            # save the df_obs subset
            df_save = df_obs[["midPointMjdTai", "outlier"]]
            print("save {}".format(save_file))
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
            # ax1.scatter(df_obs_new.loc[outlier_flag]["phaseAngle"], df_obs_new.loc[outlier_flag]["reduced_mag"], c = "r", marker = "x", s= 75, zorder = 3)
            out_mask = df_obs["outlier"] == True
            ax1.scatter(
                df_obs.loc[out_mask]["phaseAngle"],
                df_obs.loc[out_mask]["reduced_mag"],
                c="r",
                marker="x",
                s=75,
                zorder=3,
            )
            fig_file = "{}/plots/phase_curve_{}_{}.png".format(cli_args.outpath, cli_args.ssObjectId, int(t0))
            print(fig_file)
            fig = plot_errorbar(planetoid, fig=fig, filename=fig_file)

            # # when fitting, consider only the previous observations (not tonight's)
            # mjd = obs.midPointMjdTai
            # obs_mask = mjd < int(np.amax(mjd))
            # # also drop any past outlying observations?

            # print("total N obs = {}".format(len(mjd)))
            # print("number of past obs = {}".format(sum(obs_mask)))
            # print("number of tonight's obs = {}".format(sum(~obs_mask)))

            # if sum(obs_mask) < N_pc_fit:
            #     # use an assumed value of G12 until more data is available
            #     pc_fit = PhaseCurve(abs_mag=sso.H * u.mag, phase_param=0.62, model_name="HG12_Pen16")
            # else:
            #     # do a simple HG12_Pen16 fit to the past data
            #     pc_fit = pc.FitModel(alpha[obs_mask], red_mag[obs_mask], mag_err[obs_mask])
            #     pc_fit = pc.InitModelSbpy(pc_fit)

            # print(pc_fit)
            # print(pc_fit.abs_mag)
            # print(pc_fit.phase_param)

            # # now check if the new observations are outlying
            # alpha_new = alpha[~obs_mask]
            # red_mag_new = red_mag[~obs_mask]
            # mag_err_new = mag_err[~obs_mask]

            # # calculate data - model residuals
            # res = red_mag_new - pc_fit.ReducedMag(alpha_new)
            # print(res)

            # # also check for past observations that are outlying?

            # # output results:
            # # flag outlying observations
            # # output the phase curve model parameters
            # # make a plot?


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
        help="Filters required.",
        nargs="*",
        type=str,
        default=["u", "g", "r", "i", "z", "y"],
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
