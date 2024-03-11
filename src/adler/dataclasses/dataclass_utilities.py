import numpy as np


def get_data_table(sql_query, service=None, sql_filename=None):
    """Gets a table of data based on a SQL query. Table is pulled from either the RSP or a local SQL database:
    this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

    Parameters
    -----------

    sql_query : str
        The SQL query made to the RSP or SQL database.

    service : pyvo.dal.tap.TAPService object
        TAPService object linked to the RSP. Default=None.

    sql_filename : str
        Filepath to a SQL database. Default=None.

    Returns
    -----------

    data_table : DALResultsTable
        Data table containing the results of the SQL query.

    """

    if service:
        data_table = service.search(sql_query)
    elif sql_filename:
        # to-do
        pass

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
            print("Type not recognised.")
    except ValueError:
        raise ValueError("Could not cast column name to type.")
