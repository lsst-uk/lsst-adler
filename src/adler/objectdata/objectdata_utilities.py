import numpy as np
import pandas as pd
import sqlite3
import warnings
import logging
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


def flux_to_magnitude(flux, flux_err=np.nan):
    """Converts a flux measurement (with units of nanoJanskys) and its associated error
    into AB magnitudes. If no flux error is provided, the returned magnitude error
    will be NaN.

    Parameters
    -----------
    flux : astropy.units.Quantity
        Flux value with units of nanoJanskys (u.nJy).

    flux_err : astropy.units.Quantity, optional
        Flux error with units of nanoJanskys (u.nJy). Default is np.nan (dimensionless),
        in which case the magnitude error will be returned as NaN.

    Returns
    -----------
    magnitude : float
        The flux converted into AB magnitude (unitless scalar).

    magnitude_err : float
        The propagated uncertainty in AB magnitude (unitless scalar).
        Returns NaN if flux_err is not provided.

    Notes
    -----------
    - This function assumes that `flux` and `flux_err` are `astropy.units.Quantity`
      objects in units of `u.nJy`.
    """
    # TODO Handle the masked arrays better here
    # (ideally I think we want to keep magnitude as a masked array rather than making magErr non-masked)
    magnitude = flux.to(u.ABmag).value
    magnitude_err = ((2.5 / np.log(10)) * (flux_err / flux)).value
    return magnitude, magnitude_err
