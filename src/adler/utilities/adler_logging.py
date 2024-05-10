import logging
import os
from datetime import datetime


def setup_adler_logging(
    log_location,
    log_format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s ",
    log_name="",
    log_file_info="adler.log",
    log_file_error="adler.err",
):
    log = logging.getLogger(log_name)
    log_formatter = logging.Formatter(log_format)

    # comment this to suppress console output
    # stream_handler = logging.StreamHandler()
    # stream_handler.setFormatter(log_formatter)
    # log.addHandler(stream_handler)

    dstr = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    cpid = os.getpid()

    log_file_info = os.path.join(log_location, dstr + "-p" + str(cpid) + "-" + log_file_info)
    log_file_error = os.path.join(log_location, dstr + "-p" + str(cpid) + "-" + log_file_error)

    # this log will log pretty much everything: basic info, but also warnings and errors
    file_handler_info = logging.FileHandler(log_file_info, mode="w")
    file_handler_info.setFormatter(log_formatter)
    file_handler_info.setLevel(logging.INFO)
    log.addHandler(file_handler_info)

    # this log only logs warnings and errors, so they can be looked at quickly without a lot of scrolling
    file_handler_error = logging.FileHandler(log_file_error, mode="w")
    file_handler_error.setFormatter(log_formatter)
    file_handler_error.setLevel(logging.WARN)
    log.addHandler(file_handler_error)

    # I don't know why we need this line but info logging doesn't work without it, upsettingly
    log.setLevel(logging.INFO)

    return log
