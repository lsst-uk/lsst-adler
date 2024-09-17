from lsst.rsp import get_tap_service
import pandas as pd
import numpy as np
import logging
import json

from adler.dataclasses.Observations import Observations
from adler.dataclasses.MPCORB import MPCORB
from adler.dataclasses.SSObject import SSObject
from adler.dataclasses.AdlerData import AdlerData
from adler.dataclasses.dataclass_utilities import get_data_table

from adler.lasair.cassandra_fetcher import CassandraFetcher

logger = logging.getLogger(__name__)


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

        date_range : list of float
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
    def construct_from_SQL(
        cls,
        ssObjectId,
        sql_filename,
        filter_list=["u", "g", "r", "i", "z", "y"],
        date_range=[60000.0, 67300.0],
        schema=None,
    ):
        """Custom constructor which builds the AdlerPlanetoid object and the associated Observations, MPCORB and SSObject objects from
        a local SQL database. Mostly used for testing.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        sql_filename : str
            Filepath to the local SQL database.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of float
            The minimum and maximum dates of the desired observations.

        schema : str or None
            Schema/database from which to select the data tables. Can be None. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        """

        if len(date_range) != 2:
            logger.error("ValueError: date_range attribute must be of length 2.")
            raise ValueError("date_range attribute must be of length 2.")

        observations_by_filter = cls.populate_observations(
            cls, ssObjectId, filter_list, date_range, sql_filename=sql_filename, schema=schema
        )

        if len(observations_by_filter) == 0:
            logger.error(
                "No observations found for this object in the given filter(s). Check SSOID and try again."
            )
            raise Exception(
                "No observations found for this object in the given filter(s). Check SSOID and try again."
            )

        if len(filter_list) > len(observations_by_filter):
            logger.info(
                "Not all specified filters have observations. Recalculating filter list based on past observations."
            )
            filter_list = [obs_object.filter_name for obs_object in observations_by_filter]
            logger.info("New filter list is: {}".format(filter_list))

        mpcorb = cls.populate_MPCORB(cls, ssObjectId, sql_filename=sql_filename, schema=schema)
        ssobject = cls.populate_SSObject(
            cls, ssObjectId, filter_list, sql_filename=sql_filename, schema=schema
        )

        adler_data = AdlerData(ssObjectId, filter_list)

        return cls(ssObjectId, filter_list, date_range, observations_by_filter, mpcorb, ssobject, adler_data)

    @classmethod
    def construct_from_cassandra(
        cls, ssObjectId, filter_list=["u", "g", "r", "i", "z", "y"], date_range=[60000.0, 67300.0]
    ):  # pragma: no cover

        fetcher = CassandraFetcher(cassandra_hosts=["10.21.3.123"])

        MPCORB_dict = fetcher.fetch_MPCORB(ssObjectId)
        SSObject_dict = fetcher.fetch_SSObject(ssObjectId, filter_list)
        observations_dict = fetcher.fetch_observations(ssObjectId)

        # note that Cassandra doesn't allow filters/joins
        # instead we pull all observations for this ID, then filter with Pandas later
        observations_table = pd.DataFrame(observations_dict)
        observations_table.rename(columns={"decl": "dec"}, inplace=True)

        observations_by_filter = []
        for filter_name in filter_list:
            obs_slice = observations_table[
                (observations_table["band"] == filter_name)
                & (observations_table["midpointmjdtai"].between(date_range[0], date_range[1]))
            ]

            if len(obs_slice) == 0:
                logger.warning(
                    "No observations found in {} filter for this object. Skipping this filter.".format(
                        filter_name
                    )
                )
            else:
                observations = Observations.construct_from_data_table(ssObjectId, filter_name, obs_slice)
                observations_by_filter.append(observations)

        if len(observations_by_filter) == 0:
            logger.error(
                "No observations found for this object in the given filter(s). Check SSOID and try again."
            )
            raise Exception(
                "No observations found for this object in the given filter(s). Check SSOID and try again."
            )

        if len(filter_list) > len(observations_by_filter):
            logger.info(
                "Not all specified filters have observations. Recalculating filter list based on past observations."
            )
            filter_list = [obs_object.filter_name for obs_object in observations_by_filter]
            logger.info("New filter list is: {}".format(filter_list))

        mpcorb = MPCORB.construct_from_dictionary(ssObjectId, MPCORB_dict)
        ssobject = SSObject.construct_from_dictionary(ssObjectId, filter_list, SSObject_dict)

        adler_data = AdlerData(ssObjectId, filter_list)

        return cls(ssObjectId, filter_list, date_range, observations_by_filter, mpcorb, ssobject, adler_data)

    @classmethod
    def construct_from_RSP(
        cls, ssObjectId, filter_list=["u", "g", "r", "i", "z", "y"], date_range=[60000.0, 67300.0]
    ):  # pragma: no cover
        """Custom constructor which builds the AdlerPlanetoid object and the associated Observations, MPCORB and SSObject objects
        from the RSP.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of float
            The minimum and maximum dates of the desired observations.

        """

        if len(date_range) != 2:
            raise Exception("date_range argument must be of length 2.")

        service = get_tap_service("ssotap")
        logger.info("Getting past observations from DIASource/SSSource...")
        observations_by_filter = cls.populate_observations(
            cls, ssObjectId, filter_list, date_range, service=service
        )

        if len(observations_by_filter) == 0:
            logger.error(
                "No observations found for this object in the given filter(s). Check SSOID and try again."
            )
            raise Exception(
                "No observations found for this object in the given filter(s). Check SSOID and try again."
            )

        if len(filter_list) > len(observations_by_filter):
            logger.info(
                "Not all specified filters have observations. Recalculating filter list based on past observations."
            )
            filter_list = [obs_object.filter_name for obs_object in observations_by_filter]
            logger.info("New filter list is: {}".format(filter_list))

        logger.info("Populating MPCORB metadata...")
        mpcorb = cls.populate_MPCORB(cls, ssObjectId, service=service)
        logger.info("Populating SSObject metadata...")
        ssobject = cls.populate_SSObject(cls, ssObjectId, filter_list, service=service)

        adler_data = AdlerData(ssObjectId, filter_list)

        return cls(ssObjectId, filter_list, date_range, observations_by_filter, mpcorb, ssobject, adler_data)

    def populate_observations(
        self,
        ssObjectId,
        filter_list,
        date_range,
        service=None,
        sql_filename=None,
        schema="dp03_catalogs_10yr",
    ):
        """Populates the observations_by_filter class attribute. Can populate from either the RSP for a SQL database:
        this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of float
            The minimum and maximum dates of the desired observations.

        service : pyvo.dal.tap.TAPService object or None
            TAPService object linked to the RSP. Default=None.

        sql_filename : str or None
            Filepath to a SQL database. Default=None.

        schema : str or None
            Schema/database from which to select the data tables. Can be None. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        """

        if schema:  # pragma: no cover
            schema = schema + "."
        else:
            schema = ""

        observations_by_filter = []

        for filter_name in filter_list:
            observations_sql_query = f"""
                SELECT
                    ssObject.ssObjectId, mag, magErr, band, midPointMjdTai, ra, dec, phaseAngle,
                    topocentricDist, heliocentricDist, heliocentricX, heliocentricY, heliocentricZ,
                    topocentricX, topocentricY, topocentricZ,
                    eclipticLambda, eclipticBeta
                FROM
                    {schema}ssObject
                    JOIN {schema}diaSource ON {schema}ssObject.ssObjectId   = {schema}diaSource.ssObjectId
                    JOIN {schema}ssSource  ON {schema}diaSource.diaSourceId = {schema}ssSource.diaSourceId
                WHERE
                    ssObject.ssObjectId = {ssObjectId} AND band = '{filter_name}' AND midPointMjdTai BETWEEN {date_range[0]} AND {date_range[1]}
                """

            data_table = get_data_table(observations_sql_query, service=service, sql_filename=sql_filename)

            if len(data_table) == 0:
                logger.warning(
                    "No observations found in {} filter for this object. Skipping this filter.".format(
                        filter_name
                    )
                )
            else:
                observations_by_filter.append(
                    Observations.construct_from_data_table(ssObjectId, filter_name, data_table)
                )

        return observations_by_filter

    def populate_MPCORB(self, ssObjectId, service=None, sql_filename=None, schema="dp03_catalogs_10yr"):
        """Populates the MPCORB object class attribute. Can populate from either the RSP for a SQL database:
        this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        service : pyvo.dal.tap.TAPService object or None
            TAPService object linked to the RSP. Default=None.

        sql_filename : str or None
            Filepath to a SQL database. Default=None.

        schema : str or None
            Schema/database from which to select the data tables. Can be None. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        """

        if schema:  # pragma: no cover
            schema = schema + "."
        else:
            schema = ""

        MPCORB_sql_query = f"""
            SELECT
                ssObjectId, mpcDesignation, fullDesignation, mpcNumber, mpcH, mpcG, epoch, tperi, peri, node, incl, e, n, q, 
                uncertaintyParameter, flags
            FROM
                {schema}MPCORB
            WHERE
                ssObjectId = {ssObjectId}
        """

        data_table = get_data_table(MPCORB_sql_query, service=service, sql_filename=sql_filename)

        if len(data_table) == 0:
            logger.error("No MPCORB data for this object could be found for this SSObjectId.")
            raise Exception("No MPCORB data for this object could be found for this SSObjectId.")

        return MPCORB.construct_from_data_table(ssObjectId, data_table)

    def populate_SSObject(
        self, ssObjectId, filter_list, service=None, sql_filename=None, schema="dp03_catalogs_10yr"
    ):
        """Populates the SSObject class attribute. Can populate from either the RSP for a SQL database:
        this behaviour is controlled by the service and sql_filename parameters, one of which must be supplied.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        service : pyvo.dal.tap.TAPService object or None
            TAPService object linked to the RSP. Default=None.

        sql_filename : str or None
            Filepath to a SQL database. Default=None.

        schema : str or None
            Schema/database from which to select the data tables. Can be None. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        """

        if schema:  # pragma: no cover
            schema = schema + "."
        else:
            schema = ""

        filter_dependent_columns = ""

        for filter_name in filter_list:
            filter_string = "{}_H, {}_G12, {}_HErr, {}_G12Err, {}_Ndata, ".format(
                filter_name, filter_name, filter_name, filter_name, filter_name
            )

            filter_dependent_columns += filter_string

        SSObject_sql_query = f"""
            SELECT
                discoverySubmissionDate, firstObservationDate, arc, numObs, 
                {filter_dependent_columns}
                maxExtendedness, minExtendedness, medianExtendedness
            FROM
                {schema}SSObject
            WHERE
                ssObjectId = {ssObjectId}
        """

        data_table = get_data_table(SSObject_sql_query, service=service, sql_filename=sql_filename)

        if len(data_table) == 0:
            logger.error("No SSObject data for this object could be found for this SSObjectId.")
            raise Exception("No SSObject data for this object could be found for this SSObjectId.")

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
            logger.error("ValueError: Filter {} is not in AdlerPlanetoid.filter_list.".format(filter_name))
            raise ValueError("Filter {} is not in AdlerPlanetoid.filter_list.".format(filter_name))

        return self.observations_by_filter[filter_index]

    def SSObject_in_filter(self, filter_name):
        """User-friendly helper function. Returns the filter-dependent values from SSObject for a given filter.

        Parameters
        -----------
        filter_name : str
            The desired filter.

        Returns
        -----------
        ssobject_in_filter : SSObject


        """

        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            logger.error("ValueError: Filter {} is not in AdlerPlanetoid.filter_list.".format(filter_name))
            raise ValueError("Filter {} is not in AdlerPlanetoid.filter_list.".format(filter_name))

        return self.SSObject.filter_dependent_values[filter_index]

    def attach_previous_adler_data(self, filepath):
        """Attaches and returns an AdlerData object containing the most recent AdlerData
        for this ssObjectId.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.
        """

        self.PreviousAdlerData = AdlerData(self.ssObjectId, self.filter_list)
        self.PreviousAdlerData.populate_from_database(filepath)

        return self.PreviousAdlerData
