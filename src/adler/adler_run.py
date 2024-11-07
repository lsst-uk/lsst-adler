import logging
import argparse
import astropy.units as u
from astropy.time import Time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.objectdata.AdlerData import AdlerData
from adler.science.PhaseCurve import PhaseCurve, ADLER_SBPY_DICT
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
    phase_model = cli_args.phase_model  # which phase curve model to fit

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
                logger.info(msg)
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
        if cli_args.plot_show:
            plt.show()
        # Save figures at the outpath location
        else:
            fig_file = "{}/phase_curve_{}_{}_{}.png".format(
                cli_args.outpath, cli_args.ssObjectId, phase_model, int(np.amax(df_obs["midPointMjdTai"]))
            )
            msg = "Save figure: {}".format(fig_file)
            print(msg)
            logger.info(msg)
            fig = plot_errorbar(planetoid, fig=fig, filename=fig_file)  # TODO: add titles with filter name?
            plt.close()

        # Output adler values to a database if a db_name is provided
        print(adler_data.__dict__)
        if cli_args.db_name:
            adler_db = "{}/{}".format(cli_args.outpath, cli_args.db_name)
            msg = "write to {}".format(adler_db)
            print(msg)
            logger.info(msg)
            adler_data.write_row_to_database(adler_db)

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
            if cli_args.plot_show:
                plot_dir = None
            else:
                plot_dir = cli_args.outpath

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
                plot_show=cli_args.plot_show,
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
        # help="Stem filename of output database. If this doesn't exist, it will be created. Default: adler_out.",
        # type=str,
        # default="adler_out",
        help="Optional filename of output database, used to store Adler results in a db if provided.",
        type=str,
        default=None,
    )
    optional_group.add_argument(
        "-i",
        "--sql_filename",
        help="Optional input path location of a sql database file containing observations.",
        type=str,
        default=None,
    )
    optional_group.add_argument(
        "-p",
        "--plot_show",
        help="Optional flag to display plots interactively instead of saving to file.",
        action="store_true",
    )
    # TODO: Add a model_name parameter
    optional_group.add_argument(
        "-m",
        "--phase_model",
        help="Select the phase parameter model_name. LIST OPTIONS AND DEFAULT",
        type=str,
        default="HG12_Pen16",
    )

    args = parser.parse_args()

    cli_args = AdlerCLIArguments(args)

    adler_logger = setup_adler_logging(cli_args.outpath)

    cli_args.logger = adler_logger

    runAdler(cli_args)


if __name__ == "__main__":
    main()
