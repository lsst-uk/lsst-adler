import glob
import os
import pytest
import tempfile


def test_setup_adler_logging():
    from adler.utilities.adler_logging import setup_adler_logging

    with tempfile.TemporaryDirectory() as dir_name:
        logger = setup_adler_logging(dir_name)

        # Check that the files get created.
        errlog = glob.glob(os.path.join(dir_name, "*-adler.err"))
        datalog = glob.glob(os.path.join(dir_name, "*-adler.log"))

        assert os.path.exists(errlog[0])
        assert os.path.exists(datalog[0])

        # Log some information.
        logger.info("Test1")
        logger.info("Test2")
        logger.error("Error1")
        logger.info("Test3")

        # Check that all five lines exist in the INFO file.
        with open(datalog[0], "r") as f_info:
            log_data = f_info.read()
            assert "Test1" in log_data
            assert "Test2" in log_data
            assert "Error1" in log_data
            assert "Test3" in log_data

        # Check that only error and critical lines exist in the ERROR file.
        with open(errlog[0], "r") as f_err:
            log_data = f_err.read()
            assert "Test1" not in log_data
            assert "Test2" not in log_data
            assert "Error1" in log_data
            assert "Test3" not in log_data
