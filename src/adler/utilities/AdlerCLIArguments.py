class AdlerCLIArguments:
    def __init__(self, args):
        self.ssObjectId = args.ssoid
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
                "Unexpected filters found in filter_list command-line argument. filter_list must be a comma-separated list of LSST filters."
            )

    def _validate_ssObjectId(self):
        try:
            int(self.ssObjectId)
        except ValueError:
            raise ValueError("ssoid command-line argument does not appear to be a valid ssObjectId.")

    def _validate_date_range(self):
        if len(self.date_range) != 2:
            raise ValueError("date_range command-line argument must be of length 2.")
