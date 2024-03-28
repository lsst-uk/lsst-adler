import numpy as np
import pandas as pd
import sqlite3
import warnings


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

    if service:
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


def get_from_table(data_table, column_name, data_type):
    """Retrieves information from the data_table class variable and forces it to be a specified type.

    Parameters
    -----------
    column_name : str
        Column name under which the data of interest is stored.
    type : str
        String delineating data type. Should be "str", "float", "int" or "array".

    Returns
    -----------
    data : any type
        The data requested from the table cast to the type required.

    """
    try:
        if data_type == "str":
            return str(data_table[column_name][0])
        elif data_type == "float":
            return float(data_table[column_name][0])
        elif data_type == "int":
            return int(data_table[column_name][0])
        elif data_type == "array":
            return np.array(data_table[column_name])
        else:
            raise TypeError(
                "Type for argument data_type not recognised: must be one of 'str', 'float', 'int', 'array'."
            )
    except ValueError:
        raise ValueError("Could not cast column name to type.")
