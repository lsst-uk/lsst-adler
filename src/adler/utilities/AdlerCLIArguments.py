class AdlerCLIArguments:
    """
    Class for storing abd validating Adler command-line arguments.

    Attributes:
    -----------
    args : argparse.Namespace object
        argparse.Namespace object created by calling parse_args().

    """

    def __init__(self, args):
        self.ssObjectId = args.ssObjectId
        self.filter_list = args.filter_list
        self.date_range = args.date_range

        self.validate_arguments()

    def validate_arguments(self):
        self._validate_filter_list()
        self._validate_ssObjectId()
        self._validate_date_range()

    def _validate_filter_list(self):
        expected_filters = ["u", "g", "r", "i", "z", "y"]

        if not set(self.filter_list).issubset(expected_filters):
            raise ValueError(
                "Unexpected filters found in filter_list command-line argument. filter_list must be a list of LSST filters."
            )

    def _validate_ssObjectId(self):
        try:
            int(self.ssObjectId)
        except ValueError:
            raise ValueError("ssoid command-line argument does not appear to be a valid ssObjectId.")

    def _validate_date_range(self):
        for d in self.date_range:
            try:
                float(d)
            except ValueError:
                raise ValueError(
                    "One or both of the values for the date_range command-line argument do not seem to be valid numbers."
                )

        if any(d > 250000 for d in self.date_range):
            raise ValueError(
                "Dates for date_range command-line argument seem rather large. Did you input JD instead of MJD?"
            )
