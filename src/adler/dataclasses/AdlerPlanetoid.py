from lsst.rsp import get_tap_service

from adler.dataclasses.Observations import Observations
from adler.dataclasses.MPCORB import MPCORB
from adler.dataclasses.SSObject import SSObject
from adler.dataclasses.AdlerData import AdlerData
from adler.dataclasses.dataclass_utilities import get_data_table
from adler.science.DummyScience import DummyScience


class AdlerPlanetoid:
    """AdlerPlanetoid class. Contains the Observations, MPCORB and SSObject dataclass objects."""

    def __init__(
        self, ssObjectId, filter_list, date_range, observations_by_filter, mpcorb, ssobject, adler_data
    ):
        """Initialises the AdlerPlanetoid object.

        Attributes
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of int
            The minimum and maximum dates of the desired observations.

        observations_by_filter : list of Observations objects
            A list of Observations objects holding joined DIASource/SSSource observations of the planetoid specified by ssObjectId. Each item in the list holds observations of a different filter, in the order specified by filter_list.

        mpcorb : MPCORB object
            An MPCORB object, holding the MPCORB database information of the planetoid specified by ssObjectId.

        ssobject : SSObject object
            An SSObject object, holding the SSObject database information of the planetoid specified by ssObjectId.

        adler_data : AdlerData object
            An empty AdlerData object ready to store Adler-calculated values.

        """
        self.ssObjectId = ssObjectId
        self.filter_list = filter_list
        self.date_range = date_range
        self.observations_by_filter = observations_by_filter
        self.MPCORB = mpcorb
        self.SSObject = ssobject
        self.AdlerData = adler_data

    @classmethod
    def construct_from_SQL(cls):
        # to-do
        pass

    @classmethod
    def construct_from_RSP(
        cls, ssObjectId, filter_list=["u", "g", "r", "i", "z", "y"], date_range=[60000.0, 67300.0]
    ):
        """Custom constructor which builds the AdlerPlanetoid object and the associated Observations, MPCORB and SSObject objects
        from the RSP.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of int
            The minimum and maximum dates of the desired observations.

        """

        if len(date_range) != 2:
            raise Exception("date_range argument must be of length 2.")

        service = get_tap_service("ssotap")
        observations_by_filter = cls.populate_observations(
            cls, ssObjectId, filter_list, date_range, service=service
        )
        mpcorb = cls.populate_MPCORB(cls, ssObjectId, service=service)
        ssobject = cls.populate_SSObject(cls, ssObjectId, filter_list, service=service)

        adler_data = AdlerData(filter_list)

        return cls(ssObjectId, filter_list, date_range, observations_by_filter, mpcorb, ssobject, adler_data)

    def populate_observations(self, ssObjectId, filter_list, date_range, service=None, sql_filename=None):
        """Populates the observations_by_filter class attribute. Can populate from either the RSP for a SQL database:
        this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of int
            The minimum and maximum dates of the desired observations.

        service : pyvo.dal.tap.TAPService object
            TAPService object linked to the RSP. Default=None.

        sql_filename : str
            Filepath to a SQL database. Default=None.

        """

        observations_by_filter = []

        for filter_name in filter_list:
            observations_sql_query = f"""
                SELECT
                    ssObject.ssObjectId, mag, magErr, band, midpointMjdTai, ra, dec, phaseAngle,
                    topocentricDist, heliocentricDist
                FROM
                    dp03_catalogs_10yr.ssObject
                    JOIN dp03_catalogs_10yr.diaSource ON dp03_catalogs_10yr.ssObject.ssObjectId   = dp03_catalogs_10yr.diaSource.ssObjectId
                    JOIN dp03_catalogs_10yr.ssSource  ON dp03_catalogs_10yr.diaSource.diaSourceId = dp03_catalogs_10yr.ssSource.diaSourceId
                WHERE
                    ssObject.ssObjectId = {ssObjectId} AND band = '{filter_name} AND midpointMjdTai BETWEEN {date_range[0]} AND {date_range[1]}'
                """

            data_table = get_data_table(observations_sql_query, service=service, sql_filename=sql_filename)

            observations_by_filter.append(
                Observations.construct_from_data_table(ssObjectId, filter_name, data_table)
            )

        return observations_by_filter

    def populate_MPCORB(self, ssObjectId, service=None, sql_filename=None):
        """Populates the MPCORB object class attribute. Can populate from either the RSP for a SQL database:
        this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        service : pyvo.dal.tap.TAPService object
            TAPService object linked to the RSP. Default=None.

        sql_filename : str
            Filepath to a SQL database. Default=None.

        """
        MPCORB_sql_query = f"""
            SELECT
                ssObjectId, mpcDesignation, mpcNumber, mpcH, mpcG, epoch, peri, node, incl, e, n, q, 
                uncertaintyParameter, flags
            FROM
                dp03_catalogs_10yr.MPCORB
            WHERE
                ssObjectId = {ssObjectId}
        """

        data_table = get_data_table(MPCORB_sql_query, service=service, sql_filename=sql_filename)

        return MPCORB.construct_from_data_table(ssObjectId, data_table)

    def populate_SSObject(self, ssObjectId, filter_list, service=None, sql_filename=None):
        """Populates the SSObject class attribute. Can populate from either the RSP for a SQL database:
        this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        service : pyvo.dal.tap.TAPService object
            TAPService object linked to the RSP. Default=None.

        sql_filename : str
            Filepath to a SQL database. Default=None.

        """

        filter_dependent_columns = ""

        for filter_name in filter_list:
            filter_string = "{}_H, {}_G12, {}_Herr, {}_G12err, {}_nData, ".format(
                filter_name, filter_name, filter_name, filter_name, filter_name
            )

            filter_dependent_columns += filter_string

        SSObject_sql_query = f"""
            SELECT
                discoverySubmissionDate, firstObservationDate, arc, numObs, 
                {filter_dependent_columns}
                maxExtendedness, minExtendedness, medianExtendedness
            FROM
                dp03_catalogs_10yr.SSObject
            WHERE
                ssObjectId = {ssObjectId}
        """

        data_table = get_data_table(SSObject_sql_query, service=service, sql_filename=sql_filename)

        return SSObject.construct_from_data_table(ssObjectId, filter_list, data_table)

    def observations_in_filter(self, filter_name):
        """User-friendly helper function. Returns the Observations object for a given filter.

        Parameters
        -----------
        filter_name : str
            The desired filter.

        Returns
        -----------
        Observations object
            The Observations object in self.observations_by_filter corresponding to the desired filter.

        """

        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            raise ValueError("Filter {} is not in AdlerPlanetoid.filter_list.".format(filter_name))

        return self.observations_by_filter[filter_index]

    def do_pretend_science(self):
        self.DummyScienceResult = DummyScience().science_result

        print(self.DummyScienceResult)
