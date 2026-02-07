import os
from pathlib import Path


def get_test_data_path():  # pragma: no cover
    """Gets the absolute path of the tests/data directory.

    Returns
    -----------
    path_to_data : str
        The absolute path to the tests/data directory.

    """

    # where is this file?
    path_to_file = os.path.abspath(__file__)

    # the test data folder is thus:
    path_to_data = os.path.join(str(Path(path_to_file).parents[3]), "tests/data/")

    return path_to_data


def get_test_data_filepath(filename):  # pragma: no cover
    """Gets the absolute path of a supplied file in the tests/data directory.


    Parameters
    -----------
    filename : str
        The filename of the desired test file.

    Returns
    -----------
    str
        The absolute path to the specified file in the tests/data directory.

    """

    filepath = get_test_data_path()

    return os.path.join(filepath, filename)


# def fetch_adler_row_from_db(db_filepath, ssObjectId, table_name):
#     """Return the latest row for the supplied ssObjectId as a pandas Series for the given table.

#     This helper can be used in tests to inspect raw DB values before asserting them.

#     Parameters
#     -----------
#     db_filepath : str
#         Filepath to the input database.

#     ssObjectId : str or float
#         ssObjectId of interest.

#     table_name : str
#         Name of table to query.
#     """
#     if not os.path.isfile(db_filepath):
#         raise FileNotFoundError(f"Database not found: {db_filepath}")

#     con = sqlite3.connect(db_filepath)
#     df = pd.read_sql_query(
#         f"SELECT * FROM {table_name} WHERE CAST(ssObjectId AS TEXT) = '{ssObjectId}' ORDER BY updatedMJD DESC",
#         con,
#     )
#     con.close()

#     if df.empty:
#         raise ValueError("No data found for ssObjectId in database")

#     # return the most recent row as a Series
#     return df.iloc[0]


# def test_fetch_adler_row_helper():
#     # Smoke-test the DB helper to ensure it returns a Series
#     # MPC object
#     db_path = get_test_data_filepath("adler_output_median_60799.5_7n_3n.sqlite")
#     ssoid = "2025 MX40Fake"

#     row = fetch_adler_row_from_db(db_path, ssoid, "AdlerData")
#     assert isinstance(row, pd.Series)
#     # must contain ssObjectId column
#     assert "ssObjectId" in row.index

#     # DP0.3 object
#     db_path = get_test_data_filepath("adler_output_HG12_Pen16_63335.5_400n_31n.sqlite")
#     ssoid = "6098332225018000"

#     row = fetch_adler_row_from_db(db_path, ssoid, "AdlerData")
#     assert isinstance(row, pd.Series)
#     # must contain ssObjectId column
#     assert "ssObjectId" in row.index


# def fetch_adler_source_flags_from_db(db_filepath, ssObjectId):
#     """Return the latest AdlerSourceFlags data for the supplied ssObjectId as a pandas DataFrame.

#     This helper can be used in tests to inspect raw DB values before asserting them.

#     Parameters
#     -----------
#     db_filepath : str
#         Filepath to the input database.

#     ssObjectId : str or float
#         ssObjectId of interest.
#     """

#     if not os.path.isfile(db_filepath):
#         raise FileNotFoundError(f"Database not found: {db_filepath}")

#     con = sqlite3.connect(db_filepath)
#     df = pd.read_sql_query(
#         f"SELECT * FROM AdlerSourceFlags WHERE CAST(ssObjectId AS TEXT) = '{ssObjectId}'",
#         con,
#     )
#     con.close()

#     if df.empty:
#         raise ValueError("No data found for ssObjectId in database")

#     # return the most recent row as a Series
#     return df


# def test_fetch_adler_source_flags_helper():
#     # Smoke-test the DB helper to ensure it returns a DataFrame
#     # MPC object
#     db_path = get_test_data_filepath("adler_output_median_60799.5_7n_3n.sqlite")
#     ssoid = "2025 MX40Fake"

#     row = fetch_adler_source_flags_from_db(db_path, ssoid)
#     assert isinstance(row, pd.DataFrame)

#     # DP0.3 object
#     db_path = get_test_data_filepath("adler_output_HG12_Pen16_63335.5_400n_31n.sqlite")
#     ssoid = "6098332225018000"

#     row = fetch_adler_source_flags_from_db(db_path, ssoid)
#     assert isinstance(row, pd.DataFrame)
