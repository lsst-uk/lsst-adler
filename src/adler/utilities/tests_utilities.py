import os
from pathlib import Path


def get_test_data_path():
    # This file's path: `<base_directory>/src/sorcha/utilities/test_data_utilities.py`
    # THIS_DIR = `<base_directory>/`
    THIS_DIR = Path(__file__).parent.parent.parent.parent

    # Returned path: `<base_directory>/tests/data/filename`
    return os.path.join(THIS_DIR, "tests/data")


def get_test_data_filepath(filename):
    filepath = get_test_data_path()

    return os.path.join(filepath, filename)
