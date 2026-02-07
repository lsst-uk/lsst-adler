import argparse
import glob
import math
import os
import sqlite3
import subprocess
import sys
from typing import List
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.stats import sigma_clip as astropy_sigma_clip
import astropy.units as u
from astropy.time import Time
from tqdm import tqdm

from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.objectdata.objectdata_utilities import generate_summary_csvs
from adler.science.AvgMagModel import AvgMagModel
from adler.science.PhaseCurve import PhaseCurve
from adler.objectdata.AdlerData import VALID_AVG_MAG_MODELS, VALID_PHASE_MODELS
import adler.utilities.science_utilities as sci_utils
from adler.objectdata.objectdata_utilities import mpc_file_preprocessing
from adler.utilities.adler_logging import setup_adler_logging

logger = logging.getLogger(__name__)


def parse_args(argv: List[str] = None):
    p = argparse.ArgumentParser(description="Run outlier detection for a single process_mjd")
    p.add_argument("--process-mjd", type=float, help="MJD to process (e.g. 60799.5)")
    p.add_argument(
        "--process-isot",
        type=str,
        help="ISOT date to process (e.g. 2025-05-04T12:00:00.000 (YYYY-MM-DDTHH:MM:SS.sss) or 2025-05-04 (YYYY-MM-DD))",
    )
    p.add_argument(
        "--input-sql-file",
        type=str,
        default="rubin.sqlite",
        help="Path to input SQLite file",
    )
    p.add_argument("--schema", type=str, choices=["MPC", "dp03_catalogs_10yr"], default="MPC")
    p.add_argument(
        "--filter-list",
        nargs="+",
        default=["u", "g", "r", "i", "z", "y"],
        help="Filters to consider",
    )
    p.add_argument("--model-name", type=str, default="median", help="Name of model to use")
    p.add_argument(
        "--data-timespan",
        type=int,
        nargs="+",
        default=[30],
        help="Number of nights of observations to use (one or more values)",
    )
    p.add_argument(
        "--n-new-nights",
        type=int,
        nargs="+",
        default=[3],
        help="Number of nights of obsversations to consider as new nights (one or more values)",
    )
    p.add_argument("--diff-cut", type=float, default=1.5, help="Magnitude difference threshold")
    p.add_argument("--std-cut", type=float, default=5.0, help="Sigma-space threshold")
    p.add_argument("--sig-clip-val", type=float, default=3.0, help="Data sigma-clipping value")
    p.add_argument("--make-plots", action="store_true", help="Generate output plots")
    p.add_argument("--write-model-data", action="store_false", help="Generate output plots")
    p.add_argument("--output-dir", type=str, default="adler_output_files", help="Output directory")
    p.add_argument(
        "--force-overwrite",
        action="store_true",
        help="Flag to force overwriting of any existing output files",
    )
    p.add_argument("--logs-dir", type=str, default="adler_logs", help="Logs directory")
    p.add_argument(
        "--rsp-path",
        type=str,
        help="Path to where you wish to upload the output files to the RSP. This will parse your RSP username from the path specified and use the rsp_sync.sh script. See the SSSC Prompt Products Bandaid for instructions on using this script.",
    )
    return p.parse_args(argv)


def run_outliers(
    process_mjd: float,
    input_sql_file: str = "rubin.sqlite",
    schema: str = "MPC",
    filter_list: List[str] = ["u", "g", "r", "i", "z", "y"],
    model_name: str = "median",
    data_timespan_arr: List[int] = [30],
    n_new_nights_arr: List[int] = [3],
    diff_cut: float = 1.5,
    std_cut: float = 5.0,
    sig_clip_val: float = 3.0,
    make_plots: bool = False,
    write_model_data: bool = False,
    output_dir: str = "demo_output_files",
    force_overwrite: bool = False,
    rsp_path: str = None,
):
    if schema == "MPC":
        logger.info(
            "Ensuring MPC file is compatible with adler, running adler.objectdata.objectdata_utilities.mpc_file_preprocessing"
        )
        mpc_file_preprocessing(input_sql_file)

    os.makedirs(output_dir, exist_ok=True)

    conn = sqlite3.connect(input_sql_file)

    mag_col = "reduced_mag"
    magErr_col = "magErr"

    start_of_night_mjd = process_mjd - 1
    if schema == "MPC":
        q = f"SELECT DISTINCT provid FROM obs_sbn WHERE mjd_tai BETWEEN '{start_of_night_mjd}' AND '{process_mjd}'"
        obj_df = pd.read_sql_query(q, conn)
        unique_obj_ids = obj_df.provid.to_numpy()
    else:
        q = f"SELECT DISTINCT ssObjectId FROM diaSource WHERE midPointMjdTai BETWEEN '{start_of_night_mjd}' AND '{process_mjd}'"
        obj_df = pd.read_sql_query(q, conn)
        unique_obj_ids = obj_df.ssObjectId.to_numpy()

    # TODO implement unique_obj_ids.isin(<user_supplied_list_of_ids>)

    logger.info(f"{len(unique_obj_ids)} objects to analyze for {process_mjd}")
    if len(unique_obj_ids) == 0:
        logger.warning("No objects to process - exiting")
        return

    for data_timespan in data_timespan_arr:
        logger.info(f"Processing data_timespan of {data_timespan} nights")
        for n_new_nights in n_new_nights_arr:
            logger.info(f"Processing n_new_nights of {n_new_nights} nights")

            output_db = f"{output_dir}/adler_output_{model_name}_{process_mjd:.1f}_{data_timespan}n_{n_new_nights}n.sqlite"

            if os.path.isfile(output_db) and not force_overwrite:
                logger.error(
                    f"Output database already exists and force_overwrite not specified. Delete/rename the file at {output_db} or specifiy --force-overwrite in the input arguments"
                )
                raise ValueError(
                    f"Output database already exists and force_overwrite not specified. Delete/rename the file at {output_db} or specifiy --force-overwrite in the input arguments"
                )

            for ssObjectId in tqdm(unique_obj_ids, desc=f"Objects to process for {process_mjd}"):
                logger.info(f"Processing {ssObjectId}")
                if schema == "MPC":
                    planetoid = AdlerPlanetoid.construct_from_mpc_obs_sbn(
                        ssObjectId=ssObjectId,
                        sql_filename=input_sql_file,
                        filter_list=filter_list,
                        date_range=[process_mjd - data_timespan, process_mjd],
                    )
                elif schema == "dp03_catalogs_10yr":
                    planetoid = AdlerPlanetoid.construct_from_SQL(
                        ssObjectId=ssObjectId,
                        sql_filename=input_sql_file,
                        filter_list=filter_list,
                        date_range=[process_mjd - data_timespan, process_mjd],
                    )
                else:
                    logger.error("Schema not recognised")
                    raise ValueError("Schema not recognised")

                for filt in planetoid.filter_list:
                    # Set the modelId here
                    planetoid.AdlerData.set_modelId(model_name, process_mjd, data_timespan, n_new_nights)

                    df_obs = sci_utils.get_df_obs_filt(planetoid, filt=filt)

                    err_flag = df_obs.magErr.isnull().all()
                    if err_flag:
                        logger.info("All magErr values are NaNs, proceed with caution")
                    else:
                        # Remove observations with large errorbars
                        magErr_mask = sci_utils.large_magErr_mask(df_obs)
                        df_obs = df_obs[magErr_mask]

                    # Split into previous observations and observations from the most recent night(s)
                    df_obs_old, df_obs_new, *_ = sci_utils.split_obs(
                        df_obs, process_mjd=process_mjd, n_new_nights=n_new_nights
                    )
                    logger.info(
                        "Previous observations (date < {}): {}".format(
                            process_mjd - n_new_nights, len(df_obs_old)
                        )
                    )
                    logger.info(
                        "New observations ({} <= date < {}): {}".format(
                            process_mjd - n_new_nights, process_mjd, len(df_obs_new)
                        )
                    )
                    if len(df_obs_old) < 2:
                        logger.info(
                            "Insufficient number of previous observations, continuing to next band/object"
                        )
                        continue
                    if len(df_obs_new) == 0:
                        logger.info(f"No new observations in {filt}, continuing to next band/object")
                        continue

                    # Sigma clip old observations
                    sig_clip_mask = astropy_sigma_clip(df_obs_old[mag_col], sigma=sig_clip_val).mask
                    df_obs_old = df_obs_old[~sig_clip_mask].copy()
                    if len(df_obs_old) < 2:
                        logger.info(
                            "Insufficient number of previous observations after sigma clipping, continuing to next band/object"
                        )
                        continue

                    # Populate summary AdlerData params for this filter and particular model
                    ad_params = {}
                    ad_params["phaseAngle_min"] = np.amin(df_obs_old["phaseAngle"])
                    ad_params["phaseAngle_range"] = np.ptp(df_obs_old["phaseAngle"])
                    ad_params["observationTime_max"] = np.amax(df_obs_old["midPointMjdTai"])
                    ad_params["arc"] = np.ptp(df_obs_old["midPointMjdTai"])
                    ad_params["nobs"] = len(df_obs_old)

                    # Fit model
                    if model_name in VALID_AVG_MAG_MODELS:
                        model = AvgMagModel().InitModelObs(
                            mag=df_obs_old[mag_col],
                            magErr=df_obs_old[magErr_col],
                            model_name=model_name,
                        )
                        ad_params.update(model.__dict__)  # store model values in ad_params
                        planetoid.AdlerData.populate_avg_mag_parameters(filt, **ad_params)

                        res = np.array(df_obs_new[mag_col]) - model.avg_mag
                    elif model_name in VALID_PHASE_MODELS:
                        model = PhaseCurve(H=np.amin(df_obs_old["reduced_mag"]), model_name=model_name)
                        pc_fit = model.FitModel(
                            phase_angle=np.radians(df_obs_old["phaseAngle"]),
                            reduced_mag=np.array(df_obs_old["reduced_mag"]),
                            mag_err=np.array(df_obs_old[magErr_col]),
                        )
                        model = model.InitModelSbpy(pc_fit)
                        ad_params.update(model.__dict__)  # store model values in ad_params
                        planetoid.AdlerData.populate_phase_parameters(filt, **ad_params)

                        res = df_obs_new[mag_col].to_numpy() - model.ReducedMag(
                            np.radians(df_obs_new["phaseAngle"])
                        )
                    else:
                        logger.info(f"Model '{model_name}' not recognised")
                        raise ValueError(f"Model '{model_name}' not recognised")

                    ####
                    # Check for individual outlying observations
                    ####

                    # Simple magnitude difference
                    diff_cut_outlier_arr = sci_utils.outlier_diff(res, diff_cut=diff_cut)
                    df_obs_new["mag_diff"] = np.zeros(shape=len(df_obs_new), dtype=float)
                    df_obs_new.loc[diff_cut_outlier_arr, "mag_diff"] = res[diff_cut_outlier_arr]

                    # Outliers in sigma-space
                    magErrs = df_obs_new[magErr_col].to_numpy()
                    std_cut_outlier_arr = sci_utils.outlier_sigma_diff(res, magErrs, std_sigma=std_cut)
                    df_obs_new["std_diff"] = np.zeros(shape=len(df_obs_new), dtype=float)
                    res_sig_diff = res / magErrs
                    df_obs_new.loc[std_cut_outlier_arr, "std_diff"] = res_sig_diff[std_cut_outlier_arr]

                    combined_outlier_arr = diff_cut_outlier_arr | std_cut_outlier_arr

                    planetoid.AdlerData.populate_source_flags(
                        filt,
                        planetoid.AdlerData.modelId,
                        df_obs_new.loc[combined_outlier_arr],
                    )

                    ####
                    # Sustained outburst checks
                    ####

                    # Identify timegaps in case there's only one night of new data (in the case where we are using n_new_nights>1)
                    df_obs_new.sort_values(by="midPointMjdTai", inplace=True)
                    time_gaps = sci_utils.apparition_gap_finder(df_obs_new.midPointMjdTai.to_numpy(), dx=0.5)
                    if len(time_gaps) == 0:
                        # If there is only one night of new data, we continue to the next band/object
                        logger.info(
                            "Insufficient number of nights with new observations, continuing to next band/object"
                        )
                        continue

                    # Sigma clip new observations
                    sig_clip_mask_new = astropy_sigma_clip(df_obs_new[mag_col], sigma=sig_clip_val).mask
                    df_obs_new_sigclip = df_obs_new[~sig_clip_mask_new].copy()

                    if len(df_obs_new_sigclip) == 0:
                        logger.info(
                            f"No new observations in {filt} after sigma clipping, continuing to next band/object"
                        )
                        continue

                    # Calculate average magnitude of the new observations (slightly obtuse to use the model here but future proofing maybe)
                    if model_name in VALID_AVG_MAG_MODELS:
                        new_obs_model = AvgMagModel().InitModelObs(
                            mag=df_obs_new_sigclip[mag_col],
                            magErr=df_obs_new_sigclip[magErr_col],
                            model_name=model_name,
                        )
                        mag_change = np.abs(new_obs_model.avg_mag - model.avg_mag)
                        # Difference in magnitude space
                        sustained_diff_cut_flag = mag_change > diff_cut
                        # Check if sustained outlier detected and record the different between the old and new avg_mag values in AdlerData
                        if sustained_diff_cut_flag:
                            planetoid.AdlerData.populate_filter_dependent_parameters(
                                filter_name=filt, **{"sustained_outliers": mag_change}
                            )
                        else:
                            logger.info(f"No sustained outlier detected")
                    elif model_name in VALID_PHASE_MODELS:
                        logger.info(f"Sustained outlier detection in phase model regime not implemented yet")
                        pass
                        # TODO implement sustained difference in phase model regime?

                    ####
                    # Plotting
                    ####
                    # TODO make optional plotting routine into a function
                    # Make this work from planetoid and put in plotting_utilities as `plot_outliers`, perhaps with model_name as argument to control the different behaviours for AvgMag vs PhaseModel
                    if make_plots:
                        os.makedirs(f"{output_dir}/plots", exist_ok=True)
                        # FIXME sigclip sometimes removes values that shouldn't be plotted/should be marked different for the new obs?
                        fig, ax = plt.subplots()
                        if model_name in VALID_AVG_MAG_MODELS:
                            x_param = "midPointMjdTai"
                            ax.axhline(model.avg_mag, c="k", ls="-")
                            ax.axhline(new_obs_model.avg_mag, c="c", ls="--")
                        elif model_name in VALID_PHASE_MODELS:
                            x_param = "phaseAngle"
                            alpha = np.linspace(0, np.amax(df_obs.phaseAngle)) * u.deg
                            ax.plot(
                                alpha.to_value(),
                                model.ReducedMag(alpha),
                            )
                        else:
                            raise ValueError("Invalid phase model.")

                        ax.errorbar(
                            df_obs_old[x_param],
                            df_obs_old[mag_col],
                            df_obs_old[magErr_col],
                            ls="",
                            marker=".",
                            color="k",
                            label="Previous observations",
                        )
                        ax.errorbar(
                            df_obs_new[x_param],
                            df_obs_new[mag_col],
                            df_obs_new[magErr_col],
                            ls="",
                            marker=".",
                            color="c",
                            label="New observations",
                        )
                        ax.errorbar(
                            df_obs_new[diff_cut_outlier_arr][x_param],
                            df_obs_new[diff_cut_outlier_arr][mag_col],
                            df_obs_new[diff_cut_outlier_arr][magErr_col],
                            ls="",
                            marker="x",
                            color="b",
                            label="Outliers",
                        )
                        ax.invert_yaxis()
                        ax.set_xlabel(x_param)
                        ax.set_ylabel(f"{filt}-band Reduced Magnitude")
                        plot_filepath = f"{output_dir}/plots/{planetoid.ssObjectId}_{planetoid.AdlerData.modelId}_{filt}_outliers.png"
                        fig.savefig(
                            plot_filepath,
                            bbox_inches="tight",
                            pad_inches=0.05,
                        )
                        plt.close(fig)
                        logger.info(f"Plot saved to {plot_filepath}")

                # Write out AdlerData info
                planetoid.AdlerData.write_to_database(output_db, write_model_data=write_model_data)

            generate_summary_csvs(
                output_db,
                f"{output_dir}/{planetoid.AdlerData.modelId}",
                output_cols=["ssObjectId"],
                filter_list=planetoid.AdlerData.filter_list,
            )

    ####
    # Upload results to RSP
    ####
    if rsp_path:
        logger.info("Starting upload of output files to RSP...")

        # Parse RSP username from path: /home/<rsp-username>/path/to/file.sqlite
        # Extract the username from the path
        path_parts = rsp_path.split("/")
        if len(path_parts) >= 3 and path_parts[1] == "home":
            rsp_username = path_parts[2]
        else:
            logger.error(
                f"Could not parse RSP username from path: {rsp_path}. Expected /home/<rsp-username>/..."
            )
            raise ValueError(
                f"Could not parse RSP username from path: {rsp_path}. Expected /home/<rsp-username>/..."
            )

        # Find all .sqlite and .csv files in the output directory
        output_files = glob.glob(
            os.path.join(output_dir, f"*{planetoid.AdlerData.modelId}*.sqlite")
        ) + glob.glob(os.path.join(output_dir, f"*{planetoid.AdlerData.modelId}*.csv"))

        for file_path in output_files:
            file_name = os.path.basename(file_path)
            # Construct the remote path in format: <rsp-username>@rsp:<rsp-path>
            remote_path = f"{rsp_username}@rsp:{rsp_path}/{file_name}".replace(
                "//", "/"
            )  # Replace any double slashes with single slash (guard against filepath ending being directory without end slash)
            logger.info(f"Uploading {file_path} -> {remote_path}")
            try:
                subprocess.run(["./rsp_sync.sh", file_path, remote_path], check=True)
                logger.info(f"Successfully uploaded {file_name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to upload {file_name}: {e}")
                raise
        logger.info("Upload to RSP complete!")
    else:
        logger.info("No RSP path specified, skipping upload")


def main(argv=None):
    args = parse_args(argv)

    if not (args.process_mjd or args.process_isot):
        # Determine process_mjd from database
        logger.info("Processing date not specified, calculating from database max time")

        # Get the maximum obstime from the database
        conn = sqlite3.connect(args.input_sql_file)
        cur = conn.cursor()
        cur.execute("SELECT MAX(obstime) FROM obs_sbn")
        max_time_isot_utc = cur.fetchone()[0]
        conn.close()

        # Convert to TAI MJD
        max_time_mjd_tai = Time(max_time_isot_utc, format="isot", scale="utc").tai.mjd

        # Round to nearest 0.5
        args.process_mjd = math.ceil(max_time_mjd_tai - 0.5) + 0.5

        logger.info(f"Max time (ISOT UTC): {max_time_isot_utc}")
        logger.info(f"Max time (TAI MJD): {max_time_mjd_tai}")
        logger.info(f"Process MJD (rounded to nearest 0.5): {args.process_mjd}")
    elif args.process_isot:
        logger.info(f"Processing date {args.process_isot} specified in ISOT format, converting to MJD...")
        args.process_mjd = Time(args.process_isot, format="isot", scale="utc").tai.mjd

    os.makedirs(args.logs_dir, exist_ok=True)

    setup_adler_logging(
        log_location=args.logs_dir,
        log_file_info=f"adler_{args.model_name}_{args.process_mjd:.1f}.log",
        log_file_error=f"adler_{args.model_name}_{args.process_mjd:.1f}.err",
    )

    run_outliers(
        process_mjd=args.process_mjd,
        input_sql_file=args.input_sql_file,
        schema=args.schema,
        filter_list=args.filter_list,
        model_name=args.model_name,
        data_timespan_arr=args.data_timespan,
        n_new_nights_arr=args.n_new_nights,
        diff_cut=args.diff_cut,
        std_cut=args.std_cut,
        sig_clip_val=args.sig_clip_val,
        make_plots=args.make_plots,
        write_model_data=args.write_model_data,
        force_overwrite=args.force_overwrite,
        output_dir=args.output_dir,
        rsp_path=args.rsp_path,
    )


if __name__ == "__main__":
    main()
