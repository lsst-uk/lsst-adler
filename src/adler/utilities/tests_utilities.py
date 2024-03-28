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
