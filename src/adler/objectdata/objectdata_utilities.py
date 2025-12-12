import numpy as np
import pandas as pd
import sqlite3
import warnings
import logging
from astropy.time import Time
import erfa
import astropy.units as u

logger = logging.getLogger(__name__)


def get_data_table(sql_query, service=None, sql_filename=None):
    """Gets a table of data based on a SQL query. Table is pulled from either the RSP or a local SQL database:
    this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

    Parameters
    -----------

    sql_query : str
        The SQL query made to the RSP or SQL database.

    service : pyvo.dal.tap.TAPService object or None
        TAPService object linked to the RSP. Default=None.

    sql_filename : str or None
        Filepath to a SQL database. Default=None.

    Returns
    -----------

    data_table : DALResultsTable or Pandas dataframe
        Data table containing the results of the SQL query.

    """

    if service:  # pragma: no cover
        data_table = service.search(sql_query)
    elif sql_filename:
        cnx = sqlite3.connect(sql_filename)
        data_table = pd.read_sql_query(
            sql_query, cnx
        )  # would really like to move away from Pandas for this...

        # Pandas is triggering a useless FutureWarning here which I am choosing to suppress.
        # The code is already written to account for the FutureWarning, but it triggers anyway. Thanks, Pandas.
        # Note that pd.option_context("future.no_silent_downcasting", True) would also work, but only for Pandas >2.0.3.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            data_table = data_table.fillna(value=np.nan).infer_objects(
                copy=False
            )  # changes Nones to NaNs because None forces dtype=object: bad.

    return data_table


def get_from_table(data_table, column_name, data_type, table_name="default"):
    """Retrieves information from the data_table and forces it to be a specified type.

    Parameters
    -----------
    data_table : DALResultsTable or Pandas dataframe
        Data table containing columns of interest.

    column_name : str
        Column name under which the data of interest is stored.

    data_type : type
        Data type. Should be int, float, str or np.ndarray.

    table_name : str
        Name of the table. This is mostly for more informative error messages. Default="default".

    Returns
    -----------
    data_val : str, float, int or nd.array
        The data requested from the table cast to the type required.

    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=UserWarning
        )  # RSP tables mask unpopulated elements, which get converted to NaN here and trigger a warning we don't care about.
        try:
            if data_type == str:
                data_val = str(data_table[column_name][0])
            elif data_type == float:
                data_val = float(data_table[column_name][0])
            elif data_type == int:
                data_val = int(data_table[column_name][0])
            elif data_type == np.ndarray:
                data_val = np.array(data_table[column_name], ndmin=1)
            else:
                logger.error(
                    "TypeError: Type for argument data_type not recognised for column {} in table {}: must be str, float, int or np.ndarray.".format(
                        column_name, table_name
                    )
                )
                raise TypeError(
                    "Type for argument data_type not recognised for column {} in table {}: must be str, float, int or np.ndarray.".format(
                        column_name, table_name
                    )
                )
        except ValueError:
            logger.error("ValueError: Could not cast column name to type.")
            raise ValueError("Could not cast column name to type.")

    # here we alert the user if one of the values is unpopulated and change it to a NaN
    data_val = check_value_populated(data_val, data_type, column_name, table_name)

    return data_val


def get_from_dictionary(data_dict, key_name, data_type, table_name="default"):
    """Retrieves information from a dictionary and forces it to be a specified type.

    Parameters
    -----------
    data_dict : dict or dict-like object
        Dictionary containing columns of interest.

    key_name : str
        Key name under which the data of interest is stored.

    data_type : type
        Data type. Should be int, float, str or np.ndarray.

    table_name : str
        Name of the table or dictionary. This is mostly for more informative error messages. Default="default".

    Returns
    -----------
    data_val : str, float, int or nd.array
        The data requested from the dictionary cast to the type required.

    """

    try:
        if data_type == str:
            data_val = str(data_dict[key_name])
        elif data_type == float:
            data_val = float(data_dict[key_name])
        elif data_type == int:
            data_val = int(data_dict[key_name])
        elif data_type == np.ndarray:
            data_val = np.array(data_dict[key_name], ndmin=1)
        else:
            print("type not recognised")

    except ValueError:
        print("error message")

    data_val = check_value_populated(data_val, data_type, key_name, "dictionary")

    return data_val


def check_value_populated(data_val, data_type, column_name, table_name):
    """Checks to see if data_val populated properly and prints a helpful warning if it didn't.
    Usually this will trigger because the RSP or Cassandra database hasn't populated that
    field for this particular object.

    Parameters
    -----------
    data_val : str, float, int or nd.array
        The value to check.

    data_type: type
        Data type. Should be int, float, str or np.ndarray.

    column_name: str
        Column name under which the data of interest is stored.

    table_name : str
        Name of the table. This is mostly for more informative error messages. Default="default".

    Returns
    -----------
    data_val : str, float, int, nd.array or np.nan
        Either returns the original data_val or an np.nan if it detected that the value was not populated.

    """

    array_length_zero = data_type == np.ndarray and len(data_val) == 0
    number_is_nan = data_type in [float, int] and np.isnan(data_val)
    str_is_empty = data_type == str and len(data_val) == 0

    if array_length_zero or number_is_nan or str_is_empty:
        logger.warning(
            "{} unpopulated in {} table for this object. Storing NaN instead.".format(column_name, table_name)
        )
        data_val = np.nan

    return data_val


def convertTime(timestamps, input_fmt="mjd", input_scale="utc", output_fmt="mjd", output_scale="tai"):
    """
    Convenience function for converting timestamps between formats and scales.

    Parameters
    -----------
    timestamps : float, int or nd.array
        Timestamp(s) to be converted.

    input_fmt: str
        Input format of timestamps. See https://docs.astropy.org/en/latest/api/astropy.time.Time.html for allowable formats. Default 'mjd'.

    input_scale: str
        Input scale of timestamps. See https://docs.astropy.org/en/latest/api/astropy.time.Time.html for allowable scales. Default 'utc'.

    output_fmt: str
        Desired format of output values. See https://docs.astropy.org/en/latest/api/astropy.time.Time.html for allowable formats. Default 'mjd'.

    output_scale: str
        Desired scale of output values. See https://docs.astropy.org/en/latest/api/astropy.time.Time.html for allowable scales. Default 'tai'.

    Returns
    -----------
    output_timestamps : astropy.time.Time object
        The output timestamps in the desired format and scale as an astropy.time.Time object.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=".*dubious year.*", category=erfa.ErfaWarning, module="erfa.core"
        )  # Converting dates really far in the future throws a warning here that we ignore
        output_timestamps = getattr(
            getattr(Time(timestamps, format=input_fmt, scale=input_scale), output_scale), output_fmt
        )
        return output_timestamps


def sqlite_column_exists(conn, table, column):
    """
    Function for checking if column exists in a given table in a given database connection.

    Parameters
    -----------
    conn : sqlite3.Connection object
        The connection to the database that we wish to check for a given column.

    table : str
        The name of the table we wish to check for the given column in.

    column : str
        The name of the column we wish to check for.

    Returns
    -----------
    result : bool
        Boolean flag for whether the collect exists or not.
    """
    cur = conn.cursor()
    cur.execute(f"PRAGMA table_info({table})")
    result = any(row[1] == column for row in cur.fetchall())
    return result


def add_column_if_not_exists(conn, table, column, coltype):
    """
    Function for adding a column to a given table in a given database if it does not exist.

    Parameters
    -----------
    conn : sqlite3.Connection object
        The connection to the database that we wish to add the column to in the given table.

    table : str
        The name of the table we wish add the column to.

    column : str
        The name of the column we wish to add.

    coltype : str
        SQlite datatype for the column we wish to add. See https://sqlite.org/datatype3.html

    Returns
    -----------
    None
    """
    if not sqlite_column_exists(conn, table, column):
        cur = conn.cursor()
        cur.execute(f"ALTER TABLE {table} ADD COLUMN {column} {coltype};")
        logger.info(f"{column} added to {table}")
    else:
        logger.info(f"{column} already exists in {table}")


def mpc_file_preprocessing(sql_filename, jplhorizons_filename):  # pragma: no cover
    """
    Function for performing pre-processing steps on the obs_sbn table in the MPC file format.
    The function strips the leading 'L' from the band in the obs_sbn file;
    adds a mjd_tai column with the observation time in MJD in the TAI scale;
    and adds the topocentricDist, heliocentricDist and phaseAngle from JPL Horizons.

    Parameters
    -----------
    sql_filename : str
        Filepath to the local SQL database.

    jplhorizons_filename : str
        Filepath to the local csv file containing the topocentricDist, heliocentricDist and phaseAngle from JPL Horizons.

    Returns
    -----------
    None
    """
    conn = sqlite3.connect(sql_filename)
    cursor = conn.cursor()

    # Strip the leading L that is used for some of the observations in the MPC file to ensure compatbility with adler
    # TODO In future we may change this but this is an easy fix for now
    cursor.execute("UPDATE obs_sbn SET band=substr(band, 2) WHERE band LIKE 'L%';")
    conn.commit()

    logger.info(f"Leading 'L' stripped from band names")

    # Add MJD TAI column
    add_column_if_not_exists(conn, "obs_sbn", "mjd_tai", "REAL")
    cursor.execute(f"SELECT rowid, mjd_utc FROM obs_sbn;")
    rows = cursor.fetchall()
    rowids = np.array([r[0] for r in rows])
    mjd_utc_vals = np.array([r[1] for r in rows], dtype=float)
    mjd_tai_vals = convertTime(
        mjd_utc_vals, input_fmt="mjd", input_scale="utc", output_fmt="mjd", output_scale="tai"
    )
    data_to_update = list(zip(mjd_tai_vals.tolist(), rowids.tolist()))
    cursor.executemany(f"UPDATE obs_sbn SET mjd_tai = ? WHERE rowid = ?;", data_to_update)
    conn.commit()

    logger.info(f"Added mjd_tai column to obs_sbn")

    # Load in JPL horizons data
    # Check if columns with alternative names exist:
    # Once we fix the light travel time considerations this won't be technically wrong
    if sqlite_column_exists(
        conn, "obs_sbn", "r"
    ):  # heliocentricDist (without ltt corretion) is called r in some versions of this file
        cursor.execute("ALTER TABLE obs_sbn RENAME COLUMN r TO heliocentricDist")
        conn.commit()
        logger.warning(
            "Column r renamed to heliocentricDist in obs_sbn. Be wary of light travel time as this may not have been accounted for yet"
        )

    if sqlite_column_exists(
        conn, "obs_sbn", "delta"
    ):  # topocentricDist (without ltt corretion) is called delta in some versions of this file
        cursor.execute("ALTER TABLE obs_sbn RENAME COLUMN delta TO topocentricDist")
        conn.commit()
        logger.warning(
            "Column delta renamed to topocentricDist in obs_sbn. Be wary of light travel time as this may not have been accounted for yet"
        )

    if sqlite_column_exists(
        conn, "obs_sbn", "alpha"
    ):  # phaseAngle is called alpha in JPL Horizons and may not have been changed in the file
        cursor.execute("ALTER TABLE obs_sbn RENAME COLUMN alpha TO phaseAngle")
        conn.commit()
        logger.warning(
            "Column alpha renamed to phaseAngle in obs_sbn. Be wary of light travel time as this may not have been accounted for yet"
        )

    if (
        sqlite_column_exists(conn, "obs_sbn", "heliocentricDist")
        and sqlite_column_exists(conn, "obs_sbn", "topocentricDist")
        and sqlite_column_exists(conn, "obs_sbn", "phaseAngle")
    ):
        logger.info(f"heliocentricDist, topocentricDist and phaseAngle information exist in obs_sbn already.")
    else:
        add_column_if_not_exists(conn, "obs_sbn", "heliocentricDist", "REAL")
        add_column_if_not_exists(conn, "obs_sbn", "topocentricDist", "REAL")
        add_column_if_not_exists(conn, "obs_sbn", "phaseAngle", "REAL")

        jplhorizons_df = pd.read_csv(jplhorizons_filename)
        jplhorizons_df.to_sql(
            "temp_updates", conn, if_exists="replace", index=False
        )  # Create a temporary table

        cursor.execute(
            """
            UPDATE obs_sbn
            SET
            heliocentricDist = temp_updates.r,
            topocentricDist = temp_updates.delta,
            phaseAngle = temp_updates.alpha
            FROM temp_updates
            WHERE obs_sbn.obsid = temp_updates.obsid;
        """
        )
        conn.commit()

        # TODO implement a check/warning/error if not all information available in the JPL file

        cursor.execute("DROP TABLE temp_updates;")
        conn.commit()
        conn.close()

        logger.info(
            f"heliocentricDist, topocentricDist and phaseAngle information from JPL Horizons file added to obs_sbn"
        )
def flux_to_magnitude(flux, flux_err=np.nan):
    """Converts a flux measurement (with units of nanoJanskys) and its associated error
    into AB magnitudes. If no flux error is provided, the returned magnitude error
    will be NaN.

    Parameters
    -----------
    flux : float or astropy.units.Quantity
        Flux value in nanoJanskys (can be specified with Astropy units of nanoJanskys (u.nJy).

    flux_err : float or astropy.units.Quantity, optional
        Flux error with units of nanoJanskys (u.nJy). Default is np.nan (dimensionless),
        in which case the magnitude error will be returned as NaN.

    Returns
    -----------
    magnitude : float
        The flux converted into AB magnitude (unitless scalar).

    magnitude_err : float
        The propagated uncertainty in AB magnitude (unitless scalar).
        Returns NaN if flux_err is not provided.

    """
    # TODO Handle the masked arrays better here
    # (ideally I think we want to keep magnitude as a masked array rather than making magErr non-masked)
    magnitude = flux.to(u.ABmag).value
    magnitude_err = ((2.5 / np.log(10)) * (flux_err / flux)).value
    return magnitude, magnitude_err
