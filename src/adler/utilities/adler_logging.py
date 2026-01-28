import logging
import os
from datetime import datetime


def setup_adler_logging(
    log_location,
    log_format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s ",
    log_file_info="adler.log",
    log_file_error="adler.err",
):
    """Sets up two logs for Adler: info, which logs everything, and error, which logs warnings and errors.

    Parameters
    -----------
    log_location : string
        Filepath to directory in which to save logs.

    log_format : string, optional
        Format for log filename.
        Default = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s "

    log_file_info : string, optional
        Name with which to save info log. Default = "adler.log"

    log_file_error : string, optional
        Name with which to save error log. Default = "adler.err"

    Returns
    ----------
    log : logging object
        Log object.
    """

    log = logging.getLogger()  # ROOT LOGGER

    # Prevent duplicate handlers if called more than once
    if log.handlers:
        return log

    log.setLevel(logging.INFO)

    log_formatter = logging.Formatter(log_format)

    dstr = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    cpid = os.getpid()

    log_file_info = os.path.join(log_location, f"{dstr}-p{cpid}-{log_file_info}")
    log_file_error = os.path.join(log_location, f"{dstr}-p{cpid}-{log_file_error}")

    file_handler_info = logging.FileHandler(log_file_info, mode="w")
    file_handler_info.setFormatter(log_formatter)
    file_handler_info.setLevel(logging.INFO)

    file_handler_error = logging.FileHandler(log_file_error, mode="w")
    file_handler_error.setFormatter(log_formatter)
    file_handler_error.setLevel(logging.WARNING)

    log.addHandler(file_handler_info)
    log.addHandler(file_handler_error)

    return log
