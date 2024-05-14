import logging
import argparse
import astropy.units as u
import numpy as np

from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid
from adler.science.PhaseCurve import PhaseCurve
from adler.utilities.AdlerCLIArguments import AdlerCLIArguments
from adler.utilities.adler_logging import setup_adler_logging
from adler.utilities.readin_utilities import read_in_SSObjectID_file

logger = logging.getLogger(__name__)


def runAdler(cli_args):
    logger.info("Beginning Adler.")

    N_pc_fit = 10  # minimum number of data points to fit phase curve

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
    
        # # operate on each filter in turn
        for filt in planetoid.filter_list:
    
            print("fit {} filter data".format(filt))
    
            # get the filter SSObject metadata
            sso = planetoid.SSObject_in_filter(filt)
    
            # get the LSST phase curve filter model
            pc = PhaseCurve(
                abs_mag=sso.H * u.mag,
                phase_param=sso.G12,
                model_name="HG12_Pen16",
            )
            print(pc)
            print(pc.abs_mag)
            print(pc.phase_param)
    
            # get the filter observations
            obs = planetoid.observations_in_filter(filt)
            alpha = obs.phaseAngle * u.deg
            red_mag = obs.reduced_mag * u.mag
            mag_err = obs.magErr * u.mag
    
            # when fitting, consider only the previous observations (not tonight's)
            mjd = obs.midPointMjdTai
            obs_mask = mjd < int(np.amax(mjd))
            # also drop any past outlying observations?
    
            print("total N obs = {}".format(len(mjd)))
            print("number of past obs = {}".format(sum(obs_mask)))
            print("number of tonight's obs = {}".format(sum(~obs_mask)))
    
            if sum(obs_mask) < N_pc_fit:
                # use an assumed value of G12 until more data is available
                pc_fit = PhaseCurve(abs_mag=sso.H * u.mag, phase_param=0.62, model_name="HG12_Pen16")
            else:
                # do a simple HG12_Pen16 fit to the past data
                pc_fit = pc.FitModel(alpha[obs_mask], red_mag[obs_mask], mag_err[obs_mask])
                pc_fit = pc.InitModelSbpy(pc_fit)
    
            print(pc_fit)
            print(pc_fit.abs_mag)
            print(pc_fit.phase_param)
    
            # now check if the new observations are outlying
    
            # also check for past observations that are outlying?
    
            # output results
            # flag outlying observations
            # output the phase curve model parameters

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
