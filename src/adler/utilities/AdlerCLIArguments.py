import os
import logging
import numpy as np

logger = logging.getLogger(__name__)


class AdlerCLIArguments:
    """
    Class for storing and validating Adler command-line arguments.

    Attributes
    -----------
    args : argparse.Namespace object
        argparse.Namespace object created by calling parse_args().

    """

    def __init__(self, args):
        self.ssObjectId = args.ssObjectId
        self.ssObjectId_list = args.ssObjectId_list
        self.filter_list = args.filter_list
        self.colour_list = args.colour_list
        self.date_range = args.date_range
        self.outpath = args.outpath
        self.db_name = args.db_name
        self.sql_filename = args.sql_filename
        self.phase_model = args.phase_model
        self.plot_show = args.plot_show
        self.no_plot = args.no_plot

        self.validate_arguments()

    def validate_arguments(self):
        """Checks and validates the command-line arguments."""

        self._validate_filter_list()
        self._validate_date_range()
        self._validate_outpath()

        if self.ssObjectId:
            self._validate_ssObjectId()

        if self.ssObjectId_list:
            self._validate_ssObjectId_list()

        if self.sql_filename:
            self._validate_sql_filename()

        if self.colour_list:
            self._validate_colour_list()

        if self.phase_model:
            self._validate_phase_model()

        if self.plot_show or self.no_plot:
            self._validate_plot_options()

    def _validate_filter_list(self):
        """Validation checks for the filter_list command-line argument."""
        expected_filters = ["u", "g", "r", "i", "z", "y"]

        # TODO: more informative error message, show an example of required filter_list format
        if not set(self.filter_list).issubset(expected_filters):
            logging.error(
                "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
            )
            raise ValueError(
                "Unexpected filters found in --filter_list command-line argument. --filter_list must be a list of LSST filters."
            )

    def _validate_colour_list(self):
        # the expected format is a list of "g-r", "r-i" etc
        expected_filters = ["u", "g", "r", "i", "z", "y"]
        err_msg_filt1 = "Unexpected filters found in --colour_list command-line argument. --colour_list must contain LSST filters in the format 'filter2-filter1'."

        try:
            _colour_filters = np.array([x.split("-") for x in self.colour_list]).flatten()
        except:
            logging.error(err_msg_filt1)
            raise ValueError(err_msg_filt1)

        if not set(_colour_filters).issubset(expected_filters):
            logging.error(err_msg_filt1)
            raise ValueError(err_msg_filt1)

        if not set(_colour_filters).issubset(self.filter_list):
            err_msg = "The filters required to calculate the colours have not been requested in --filter-list"
            logging.error(err_msg)
            raise ValueError(err_msg)

    def _validate_ssObjectId(self):
        """
        Validation checks for the ssObjectId command-line argument.
        """
        try:
            int(self.ssObjectId)
        except ValueError:
            logging.error("--ssObjectId command-line argument does not appear to be a valid ssObjectId.")
            raise ValueError("--ssObjectId command-line argument does not appear to be a valid ssObjectId.")

    def _validate_date_range(self):
        """
        Validation checks for the date_range command-line argument.
        """
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
        """
        Validation checks for the outpath command-line argument.
        """
        # make it an absolute path if it's relative!
        self.outpath = os.path.abspath(self.outpath)

        if not os.path.isdir(self.outpath):
            logging.error("The output path for the command-line argument --outpath cannot be found.")
            raise ValueError("The output path for the command-line argument --outpath cannot be found.")

    def _validate_ssObjectId_list(self):
        """
        Validation checks for the ssObjectId_list command-line argument.
        """
        self.ssObjectId_list = os.path.abspath(self.ssObjectId_list)

        if not os.path.exists(self.ssObjectId_list):
            logging.error(
                "The file supplied for the command-line argument --ssObjectId_list cannot be found."
            )
            raise ValueError(
                "The file supplied for the command-line argument --ssObjectId_list cannot be found."
            )

    def _validate_sql_filename(self):
        """
        Validation checks for the sel_filename command-line argument.
        """
        self.sql_filename = os.path.abspath(self.sql_filename)

        if not os.path.exists(self.sql_filename):
            logging.error("The file supplied for the command-line argument --sql_filename cannot be found.")
            raise ValueError(
                "The file supplied for the command-line argument --sql_filename cannot be found."
            )

    def _validate_phase_model(self):
        """Validation checks for the phase_model command-line argument."""
        expected_models = ["HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc"]
        err_msg_model = (
            "Unexpected model in --phase_model command-line arguments. Please select from {}".format(
                expected_models
            )
        )

        if self.phase_model not in expected_models:
            logging.error(err_msg_model)
            raise ValueError(err_msg_model)

    def _validate_plot_options(self):
        """Validation checks for the plotting options."""

        err_msg_plot = "Both the --plot_show and --no_plot flags have been selected. Please choose only one."

        if self.plot_show and self.no_plot:
            logging.error(err_msg_plot)
            raise ValueError(err_msg_plot)
