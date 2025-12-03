from lsst.rsp import get_tap_service
import pandas as pd
import numpy as np
import logging
import json
import astropy.units as u

from adler.objectdata.Observations import Observations
from adler.objectdata.MPCORB import MPCORB
from adler.objectdata.SSObject import SSObject
from adler.objectdata.AdlerData import AdlerData
from adler.objectdata.objectdata_utilities import get_data_table, flux_to_magnitude

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
        cls,
        ssObjectId,
        filter_list=["u", "g", "r", "i", "z", "y"],
        date_range=[60000.0, 67300.0],
        cassandra_hosts=["10.21.3.123"],
    ):  # pragma: no cover
        """Custom constructor which builds the AdlerPlanetoid object and the associated Observations, MPCORB and SSObject objects from
        a Cassandra database. Used only for Lasair integration.

        TODO: move method to its own class which inherits from AdlerPlanetoid and move to adler-lasair repo?

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of float
            The minimum and maximum dates of the desired observations.

        cassandra_hosts : list of str
            Location of the Cassandra database - usually an IP address. Default is ["10.21.3.123"].

        """
        # do not move this import! CassandraFetcher requires the non-mandatory
        # cassandra-driver library - if not installed, and this import is at the top,
        # test collection will break.
        from adler.lasair.cassandra_fetcher import CassandraFetcher

        fetcher = CassandraFetcher(cassandra_hosts=cassandra_hosts)

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
        cls,
        ssObjectId,
        filter_list=["u", "g", "r", "i", "z", "y"],
        date_range=[60000.0, 67300.0],
        schema="dp03_catalogs_10yr",
        flux_flag=None,
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

        schema : str or None
            Schema/database from which to select the data tables. Can be None. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        flux_flag : str or None
            Name of the flux column to select from DP1 DiaSource table. Determines FluxErr and ra/dec columns to select also. Default is None (selects mag/magErr/ra/dec for DP0.3)

        """

        if len(date_range) != 2:
            raise Exception("date_range argument must be of length 2.")

        # Select correct TAP service depending on schema chosen
        if schema == "dp03_catalogs_10yr":
            service = get_tap_service("ssotap")
        elif schema == "dp1":
            service = get_tap_service("tap")
        else:
            logger.error(f"Schema {schema} not recognised.")
            raise Exception(f"Schema {schema} not recognised.")
        logger.info("Getting past observations from DIASource/SSSource...")

        observations_by_filter = cls.populate_observations(
            cls, ssObjectId, filter_list, date_range, service=service, schema=schema, flux_flag=flux_flag
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
        mpcorb = cls.populate_MPCORB(cls, ssObjectId, service=service, schema=schema)
        logger.info("Populating SSObject metadata...")
        ssobject = cls.populate_SSObject(cls, ssObjectId, filter_list, service=service, schema=schema)

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
        flux_flag=None,
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

        flux_flag : str or None
            Name of the flux column to select from DP1 DiaSource table. Determines FluxErr and ra/dec columns to select also. Default is None (selects mag/magErr/ra/dec for DP0.3)
        """

        if schema:  # pragma: no cover
            sql_schema = schema + "."
        else:
            sql_schema = ""

        # Convenient dict for setting which columns to include in SQL query given schema and desired flux flag
        # TODO better handling of None case (this is probably bad Python)
        schema_config_dict = {
            None: {
                None: dict(
                    fluxmag_column="mag", fluxmag_err_column="magErr", ra_column="ra", dec_column="dec"
                )
            },
            "dp03_catalogs_10yr": {
                None: dict(
                    fluxmag_column="mag", fluxmag_err_column="magErr", ra_column="ra", dec_column="dec"
                )
            },
            "dp1": {
                "apFlux": dict(
                    fluxmag_column="apFlux", fluxmag_err_column="apFluxErr", ra_column="ra", dec_column="dec"
                ),
                "trailFlux": dict(
                    fluxmag_column="trailFlux",
                    fluxmag_err_column="psfFluxErr",
                    ra_column="trailRa",
                    dec_column="trailDec",
                ),
                "psfFlux": dict(
                    fluxmag_column="psfFlux",
                    fluxmag_err_column="psfFluxErr",
                    ra_column="ra",
                    dec_column="dec",
                ),
            },
        }

        try:
            selected_config = schema_config_dict[schema][flux_flag]
        except KeyError:
            if schema not in schema_config_dict:
                logger.error(f"Schema {schema} not recognised.")
                raise Exception(f"Schema {schema} not recognised.")
            else:
                logger.error(f"Flux column {flux_flag} not recognised for schema {schema}.")
                raise Exception(f"Flux column {flux_flag} not recognised for schema {schema}.")

        fluxmag_column = selected_config["fluxmag_column"]
        fluxmag_err_column = selected_config["fluxmag_err_column"]
        ra_column = selected_config["ra_column"]
        dec_column = selected_config["dec_column"]

        observations_by_filter = []

        for filter_name in filter_list:
            observations_sql_query = f"""
                SELECT
                    SSObject.ssObjectId, SSSource.diaSourceId, {fluxmag_column}, {fluxmag_err_column}, band, midPointMjdTai, {ra_column} AS ra, {dec_column} AS dec, phaseAngle,
                    topocentricDist, heliocentricDist, heliocentricX, heliocentricY, heliocentricZ,
                    topocentricX, topocentricY, topocentricZ,
                    eclipticLambda, eclipticBeta
                FROM
                    {sql_schema}SSObject
                    JOIN {sql_schema}DiaSource ON {sql_schema}SSObject.ssObjectId   = {sql_schema}DiaSource.ssObjectId
                    JOIN {sql_schema}SSSource  ON {sql_schema}DiaSource.diaSourceId = {sql_schema}SSSource.diaSourceId
                WHERE
                    SSObject.ssObjectId = {ssObjectId} AND band = '{filter_name}' AND midPointMjdTai BETWEEN {date_range[0]} AND {date_range[1]}
                """

            # This function submits the query and gets the results (or pulls from the SQL database)
            data_table = get_data_table(observations_sql_query, service=service, sql_filename=sql_filename)

            if len(data_table) == 0:
                logger.warning(
                    "No observations found in {} filter for this object. Skipping this filter.".format(
                        filter_name
                    )
                )
            else:
                if schema in [None, "dp03_catalogs_10yr"]:  # TODO probably better way to do this
                    observations_by_filter.append(
                        Observations.construct_from_data_table(ssObjectId, filter_name, data_table)
                    )
                elif schema == "dp1":
                    # Convert to astropy table so we can operate on it and add mag,magErr columns
                    data_table_astropy = data_table.to_table()

                    # Compute magnitudes
                    mag, mag_err = flux_to_magnitude(
                        data_table_astropy[fluxmag_column], data_table_astropy[fluxmag_err_column]
                    )

                    # Insert the new columns at the same positions
                    data_table_astropy.add_column(
                        mag, name="mag", index=data_table_astropy.colnames.index(fluxmag_column)
                    )
                    data_table_astropy.add_column(
                        mag_err, name="magErr", index=data_table_astropy.colnames.index(fluxmag_err_column)
                    )

                    # Remove the old flux columns
                    data_table_astropy.remove_columns([fluxmag_column, fluxmag_err_column])

                    observations_by_filter.append(
                        Observations.construct_from_data_table(ssObjectId, filter_name, data_table_astropy)
                    )
                else:
                    logger.error(f"Schema {schema} not recognised.")
                    raise Exception(f"Schema {schema} not recognised.")

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
            sql_schema = schema + "."
        else:
            sql_schema = ""

        if schema in [None, "dp03_catalogs_10yr"]:
            # Query for DP0.3. Compatible with subsequent adler code
            MPCORB_sql_query = f"""
                SELECT
                    ssObjectId, mpcDesignation, fullDesignation, mpcNumber, mpcH, mpcG, epoch, tperi, peri, node, incl, e, n, q, uncertaintyParameter, flags
                FROM
                    {sql_schema}MPCORB
                WHERE
                    ssObjectId = {ssObjectId}
            """
        elif schema == "dp1":
            # Query for DP1. Selecting the columns that still exist in the DP1 table
            # We select t_p (MJD of pericentric passage) as tperi for consistency with DP0.3
            MPCORB_sql_query = f"""
                SELECT
                    ssObjectId, mpcDesignation, mpcH, epoch, t_p AS tperi, peri, node, incl, e, q
                FROM
                    {sql_schema}MPCORB
                WHERE
                    ssObjectId = {ssObjectId}
            """
        else:
            logger.error(f"Schema {schema} not recognised.")
            raise Exception(f"Schema {schema} not recognised.")

        data_table = get_data_table(MPCORB_sql_query, service=service, sql_filename=sql_filename)

        if len(data_table) == 0:
            logger.error("No MPCORB data for this object could be found for this SSObjectId.")
            raise Exception("No MPCORB data for this object could be found for this SSObjectId.")

        if schema in [None, "dp03_catalogs_10yr"]:
            return MPCORB.construct_from_data_table(ssObjectId, data_table)
        elif schema == "dp1":
            # TODO get_data_table (above) NaN fills if we, e.g., SELECT NULL AS mpcNumber, which may be fine and remove the need for this
            # Convert to astropy Table and add in NaNs/0/empty strings for the columns that do not appear in DP1
            data_table_astropy = data_table.to_table()
            data_table_astropy.add_columns(
                cols=[
                    np.full(len(data_table_astropy), ""),  # fullDesignation (str)
                    np.full(len(data_table_astropy), 0),  # mpcNumber (int)
                    np.full(len(data_table_astropy), np.nan),  # mpcG (float)
                    np.full(len(data_table_astropy), np.nan),  # n (float)
                    np.full(len(data_table_astropy), ""),  # uncertaintyParameter (str)
                    np.full(
                        len(data_table_astropy), ""
                    ),  # flags (str)
                ],
                names=["fullDesignation", "mpcNumber", "mpcG", "n", "uncertaintyParameter", "flags"],
            )

            # Reorder columns to match DP0.3 expected order
            data_table_astropy = data_table_astropy[
                [
                    "ssObjectId",
                    "mpcDesignation",
                    "fullDesignation",
                    "mpcNumber",
                    "mpcH",
                    "mpcG",
                    "epoch",
                    "tperi",
                    "peri",
                    "node",
                    "incl",
                    "e",
                    "n",
                    "q",
                    "uncertaintyParameter",
                    "flags",
                ]
            ]
            return MPCORB.construct_from_data_table(ssObjectId, data_table_astropy)
        else:
            logger.error(f"Schema {schema} not recognised.")
            raise Exception(f"Schema {schema} not recognised.")

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
            sql_schema = schema + "."
        else:
            sql_schema = ""

        filter_dependent_columns = ""

        for filter_name in filter_list:
            filter_string = "{}_H, {}_G12, {}_HErr, {}_G12Err, {}_Ndata, ".format(
                filter_name, filter_name, filter_name, filter_name, filter_name
            )

            filter_dependent_columns += filter_string

        if schema in [None, "dp03_catalogs_10yr"]:
            # Query for DP0.3. Compatible with subsequent adler code
            SSObject_sql_query = f"""
                SELECT
                    discoverySubmissionDate, firstObservationDate, arc, numObs, 
                    {filter_dependent_columns}
                    maxExtendedness, minExtendedness, medianExtendedness
                FROM
                    {sql_schema}SSObject
                WHERE
                    ssObjectId = {ssObjectId}
            """
        elif schema == "dp1":
            # Query for DP1. Selecting the columns that still exist in the DP1 table
            SSObject_sql_query = f"""
                SELECT
                    discoverySubmissionDate, numObs
                FROM
                    {sql_schema}SSObject
                WHERE
                    ssObjectId = {ssObjectId}
            """
        else:
            logger.error(f"Schema {schema} not recognised.")
            raise Exception(f"Schema {schema} not recognised.")

        data_table = get_data_table(SSObject_sql_query, service=service, sql_filename=sql_filename)

        if len(data_table) == 0:
            logger.error("No SSObject data for this object could be found for this SSObjectId.")
            raise Exception("No SSObject data for this object could be found for this SSObjectId.")

        if schema in [None, "dp03_catalogs_10yr"]:
            return SSObject.construct_from_data_table(ssObjectId, filter_list, data_table)
        elif schema == "dp1":
            # Convert to Table
            data_table_astropy = data_table.to_table()

            # Add non-filter-dependent columns and populate with NaNs
            data_table_astropy.add_columns(
                cols=[
                    np.full(len(data_table_astropy), np.nan),  # firstObservationDate
                    np.full(len(data_table_astropy), np.nan),  # arc
                    np.full(len(data_table_astropy), np.nan),  # maxExtendedness
                    np.full(len(data_table_astropy), np.nan),  # minExtendedness
                    np.full(len(data_table_astropy), np.nan),  # medianExtendedness
                ],
                names=[
                    "firstObservationDate",
                    "arc",
                    "maxExtendedness",
                    "minExtendedness",
                    "medianExtendedness",
                ],
            )

            # Add all filter-dependent columns and populate with NaNs
            for filter_name in filter_list:
                data_table_astropy.add_columns(
                    cols=[
                        np.full(len(data_table_astropy), np.nan),  # f"{filter_name}_H"
                        np.full(len(data_table_astropy), np.nan),  # f"{filter_name}_G12"
                        np.full(len(data_table_astropy), np.nan),  # f"{filter_name}_HErr"
                        np.full(len(data_table_astropy), np.nan),  # f"{filter_name}_G12Err"
                        np.full(len(data_table_astropy), 0),  # Ndata
                    ],
                    names=[
                        f"{filter_name}_H",
                        f"{filter_name}_G12",
                        f"{filter_name}_HErr",
                        f"{filter_name}_G12Err",
                        f"{filter_name}_Ndata",
                    ],
                )

            # Reorder columns to match DP0.3 expected order
            dp03_cols_order = ["discoverySubmissionDate", "firstObservationDate", "arc", "numObs"]
            for filter_name in filter_list:
                dp03_cols_order += [
                    f"{filter_name}_H",
                    f"{filter_name}_G12",
                    f"{filter_name}_HErr",
                    f"{filter_name}_G12Err",
                    f"{filter_name}_Ndata",
                ]
            dp03_cols_order += ["maxExtendedness", "minExtendedness", "medianExtendedness"]

            data_table_astropy = data_table_astropy[dp03_cols_order]

            return SSObject.construct_from_data_table(ssObjectId, filter_list, data_table_astropy)
        else:
            logger.error(f"Schema {schema} not recognised.")
            raise Exception(f"Schema {schema} not recognised.")

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
