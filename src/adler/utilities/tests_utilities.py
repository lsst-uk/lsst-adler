import os
from pathlib import Path


def get_test_data_path():
    # where is this file?
    path_to_file = os.path.abspath(__file__)

    # the test data folder is thus:
    path_to_data = os.path.join(str(Path(path_to_file).parents[3]), "tests/data/")

    return path_to_data


def get_test_data_filepath(filename):
    filepath = get_test_data_path()

    return os.path.join(filepath, filename)
