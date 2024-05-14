import logging
import argparse
import astropy.units as u

from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid
from adler.science.PhaseCurve import PhaseCurve
from adler.utilities.AdlerCLIArguments import AdlerCLIArguments
from adler.utilities.adler_logging import setup_adler_logging
from adler.utilities.readin_utilities import read_in_SSObjectID_file

logger = logging.getLogger(__name__)


def runAdler(cli_args):
    logger.info("Beginning Adler.")

    if cli_args.ssObjectId_list:
        ssObjectId_list = read_in_SSObjectID_file(cli_args.ssObjectId_list)
    else:
        ssObjectId_list = [cli_args.ssObjectId]

    for i, ssObjectId in enumerate(ssObjectId_list):
        logger.info("Processing object {}/{}.".format(i + 1, len(ssObjectId_list)))
        logger.info("Ingesting all data for object {} from RSP...".format(cli_args.ssObjectId))

        planetoid = AdlerPlanetoid.construct_from_RSP(ssObjectId, cli_args.filter_list, cli_args.date_range)

        logger.info("Data successfully ingested.")
        logger.info("Calculating phase curves...")

        # now let's do some phase curves!

        # get the r filter SSObject metadata
        sso_r = planetoid.SSObject_in_filter("r")

        # get the RSP r filter model
        pc = PhaseCurve(
            abs_mag=sso_r.H * u.mag,
            phase_param=sso_r.G12,
            model_name="HG12_Pen16",
        )
        print(pc)
        print(pc.abs_mag, pc.phase_param)

        # get the r filter observations
        obs_r = planetoid.observations_in_filter("r")
        alpha = obs_r.phaseAngle * u.deg
        red_mag = obs_r.reduced_mag * u.mag
        mag_err = obs_r.magErr * u.mag

        # do a simple fit to all data
        pc_fit = pc.FitModel(alpha, red_mag, mag_err)
        print(pc_fit)


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

    args = parser.parse_args()

    cli_args = AdlerCLIArguments(args)

    adler_logger = setup_adler_logging(cli_args.outpath)

    cli_args.logger = adler_logger

    runAdler(cli_args)


if __name__ == "__main__":
    main()
