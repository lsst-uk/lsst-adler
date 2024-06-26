import os
import logging

logger = logging.getLogger(__name__)


class AdlerCLIArguments:
    """
    Class for storing abd validating Adler command-line arguments.

    Attributes
    -----------
    args : argparse.Namespace object
        argparse.Namespace object created by calling parse_args().

    """

    def __init__(self, args):
        self.ssObjectId = args.ssObjectId
        self.ssObjectId_list = args.ssObjectId_list
        self.filter_list = args.filter_list
        self.date_range = args.date_range
        self.outpath = args.outpath
        self.db_name = args.db_name
        self.sql_filename = args.sql_filename

        self.validate_arguments()

    def validate_arguments(self):
        self._validate_filter_list()
        self._validate_date_range()
        self._validate_outpath()

        if self.ssObjectId:
            self._validate_ssObjectId()

        if self.ssObjectId_list:
            self._validate_ssObjectId_list()

        if self.sql_filename:
            self._validate_sql_filename()

    def _validate_filter_list(self):
        expected_filters = ["u", "g", "r", "i", "z", "y"]

        if not set(self.filter_list).issubset(expected_filters):
            logging.error(
                "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
            )
            raise ValueError(
                "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
            )

    def _validate_ssObjectId(self):
        try:
            int(self.ssObjectId)
        except ValueError:
            logging.error("--ssObjectId command-line argument does not appear to be a valid ssObjectId.")
            raise ValueError("--ssObjectId command-line argument does not appear to be a valid ssObjectId.")

    def _validate_date_range(self):
        for d in self.date_range:
            try:
                float(d)
            except ValueError:
                logging.error(
                    "One or both of the values for the --date_range command-line argument do not seem to be valid numbers."
                )
                raise ValueError(
                    "One or both of the values for the --date_range command-line argument do not seem to be valid numbers."
                )

        if any(d > 250000 for d in self.date_range):
            logging.error(
                "Dates for --date_range command-line argument seem rather large. Did you input JD instead of MJD?"
            )
            raise ValueError(
                "Dates for --date_range command-line argument seem rather large. Did you input JD instead of MJD?"
            )

    def _validate_outpath(self):
        # make it an absolute path if it's relative!
        self.outpath = os.path.abspath(self.outpath)

        if not os.path.isdir(self.outpath):
            logging.error("The output path for the command-line argument --outpath cannot be found.")
            raise ValueError("The output path for the command-line argument --outpath cannot be found.")

    def _validate_ssObjectId_list(self):
        self.ssObjectId_list = os.path.abspath(self.ssObjectId_list)

        if not os.path.exists(self.ssObjectId_list):
            logging.error(
                "The file supplied for the command-line argument --ssObjectId_list cannot be found."
            )
            raise ValueError(
                "The file supplied for the command-line argument --ssObjectId_list cannot be found."
            )

    def _validate_sql_filename(self):
        self.sql_filename = os.path.abspath(self.sql_filename)

        if not os.path.exists(self.sql_filename):
            logging.error("The file supplied for the command-line argument --sql_filename cannot be found.")
            raise ValueError(
                "The file supplied for the command-line argument --sql_filename cannot be found."
            )
