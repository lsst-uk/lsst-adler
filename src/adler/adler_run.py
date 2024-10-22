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


def runAdler(cli_args):
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

        # make an empty figure
        fig = plot_errorbar(planetoid, filt_list=[])

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

            # Determine the reference phase curve model
            # TODO: We would load the best phase curve model available in AdlerData here

            # initial simple phase curve filter model with fixed G12
            pc = PhaseCurve(
                H=sso.H * u.mag,
                phase_parameter_1=0.62,
                model_name="HG12_Pen16",
            )

            # only fit G12 when sufficient data is available
            if len(df_obs) < N_pc_fit:
                msg = "Do not fit G12, use G12={:.2f}".format(pc.phase_parameter_1)
                logger.info(msg)
                pc.model_function.G12.fixed = True
            else:
                pc.model_function.G12.fixed = False

            # do a HG12_Pen16 fit to the past data
            pc_fit = pc.FitModel(
                np.array(df_obs["phaseAngle"]) * u.deg,
                np.array(df_obs["reduced_mag"]) * u.mag,
                np.array(df_obs["magErr"]) * u.mag,
            )
            pc_fit = pc.InitModelSbpy(pc_fit)

            # Store the fitted values in an AdlerData object
            adler_data.populate_phase_parameters(filt, **pc_fit.__dict__)

            # add to plot
            ax1 = fig.axes[0]
            # TODO: set colours based on filter
            ax1.scatter(df_obs["phaseAngle"], df_obs["reduced_mag"])
            alpha = np.linspace(0, np.amax(obs.phaseAngle)) * u.deg
            ax1.plot(
                alpha.value,
                pc_fit.ReducedMag(alpha).value,
                label="{}, H={:.2f}, G12={:.2f}".format(filt, pc_fit.H.value, pc_fit.phase_parameter_1),
            )

        # TODO: save the figures if an outpath is provided
        ax1.legend()
        if cli_args.outpath:
            fig_file = "{}/phase_curve_{}_{}.png".format(
                cli_args.outpath, cli_args.ssObjectId, int(np.amax(df_obs["midPointMjdTai"]))
            )
            msg = "Save figure: {}".format(fig_file)
            print(msg)
            logger.info(msg)
            fig = plot_errorbar(planetoid, fig=fig, filename=fig_file)  # TODO: add titles with filter name?
            plt.close()
        else:
            plt.show()

        # TODO: output adler values to a database
        print(adler_data.__dict__)

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

            # determine the filt_obs - filt_ref colour
            # generate a plot
            if cli_args.outpath:
                plot_dir = cli_args.outpath
            else:
                plot_dir = None

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
            )

            print(col_dict)


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
        help="Output path location. Default is current working directory.",  # TODO: make adler create the outpath directory on start up if it does not exist?
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
