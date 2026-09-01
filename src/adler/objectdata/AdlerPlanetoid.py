# from lsst.rsp import get_tap_service
from lsst.rsp import RSPDiscovery
import pandas as pd
import numpy as np
import logging
import json
import astropy.units as u
from astropy.table import Table
import os

from adler.objectdata.Observations import Observations
from adler.objectdata.MPCORB import MPCORB
from adler.objectdata.SSObject import SSObject
from adler.objectdata.AdlerData import AdlerData
from adler.objectdata.objectdata_utilities import get_data_table, flux_to_magnitude, get_tap_service_api

logger = logging.getLogger(__name__)

# Load the adler schema file that maps different input schema onto adler
schema_file = os.path.join(
    os.path.dirname(__file__), "adler_schema_map.csv"
)  # ensure the filepath is always relative to the location of this file, AdlerPlanetoid.py
ADLER_SCHEMA = pd.read_csv(schema_file, index_col=0).to_dict()

# Convenient dict for setting which columns to include in SQL query given schema and desired flux flag
# Also defines the name of the MPCORB table (dp03, dp1: MPCORB; dp2: mpc_orbits) and it id column (dp03, dp1: ssObjectId; dp2: designation)
# TODO: better handling of None case (this is probably bad Python)
SCHEMA_CONFIG_DICT = {
    # None: {None: dict(fluxmag_column="mag", fluxmag_err_column="magErr", ra_column="ra", dec_column="dec")},
    "dp03_catalogs_10yr": {
        None: dict(
            fluxmag_column="mag",
            fluxmag_err_column="magErr",
            fluxunit=u.nJy,
            ra_column="ra",
            dec_column="dec",
        ),
        "mpc_table": "MPCORB",
        "mpc_id": "ssObjectId",
    },
    "dp1": {
        "apFlux": dict(
            fluxmag_column="apFlux",
            fluxmag_err_column="apFluxErr",
            fluxunit=u.nJy,
            ra_column="ra",
            dec_column="dec",
        ),
        "trailFlux": dict(
            fluxmag_column="trailFlux",
            fluxmag_err_column="psfFluxErr",  # TODO: warn that DP1 does not have an uncertainty for trailFlux?
            fluxunit=u.nJy,
            ra_column="trailRa",
            dec_column="trailDec",
        ),
        "psfFlux": dict(
            fluxmag_column="psfFlux",
            fluxmag_err_column="psfFluxErr",
            fluxunit=u.nJy,
            ra_column="ra",
            dec_column="dec",
        ),
        "mpc_table": "MPCORB",
        "mpc_id": "ssObjectId",
    },
    "dp2": {
        "apFlux": dict(
            fluxmag_column="apFlux",
            fluxmag_err_column="apFluxErr",
            fluxunit=u.nJy,
            ra_column="ra",
            dec_column="dec",
        ),
        "trailFlux": dict(
            fluxmag_column="trailFlux",
            fluxmag_err_column="trailFluxErr",
            fluxunit=u.nJy,
            ra_column="trailRa",
            dec_column="trailDec",
        ),
        "psfFlux": dict(
            fluxmag_column="psfFlux",
            fluxmag_err_column="psfFluxErr",
            fluxunit=u.nJy,
            ra_column="ra",
            dec_column="dec",
        ),
        "mpc_table": "mpc_orbits",
        "mpc_id": "designation",
    },
    "MPC": {
        None: dict(
            fluxmag_column="mag",
            fluxmag_err_column="magErr",
            fluxunit=u.nJy,
            ra_column="ra",
            dec_column="dec",
        ),
        "mpc_table": "mpc_orbits",
        "mpc_id": "fullDesignation",
    },
}

# Define the tap service setup for each schema for queries on the RSP
RSP_TAP_CONFIG_DICT = {"dp03_catalogs_10yr": "ssotap", "dp1": "tap", "dp2": "tap"}


class AdlerPlanetoid:
    """AdlerPlanetoid class. Contains the Observations, MPCORB and SSObject dataclass objects."""

    def __init__(
        self,
        ssObjectId,
        filter_list,
        date_range,
        observations_by_filter,
        mpcorb,
        ssobject,
        adler_data,
    ):
        """Initialises the AdlerPlanetoid object.

        Attributes
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of float or None
            Optional. The minimum and maximum dates of the desired observations (MJD), e.g. [60000.0, 67300.0]

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
        date_range=None,
        schema="dp03_catalogs_10yr",
        flux_flag=None,
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

        date_range : list of float or None
            Optional. The minimum and maximum dates of the desired observations (MJD), e.g. [60000.0, 67300.0]

        schema : str
            Schema/database from which to select the data tables. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        flux_flag : str or None
            Name of the flux column to select from DP1 DiaSource table. Determines FluxErr and ra/dec columns to select also. Default is None (selects mag/magErr/ra/dec for DP0.3)

        """

        if date_range is not None:
            if len(date_range) != 2:
                logger.error("ValueError: date_range attribute must be of length 2.")
                raise ValueError("date_range attribute must be of length 2.")

        observations_by_filter = cls.populate_observations(
            cls,
            ssObjectId,
            filter_list,
            date_range,
            sql_filename=sql_filename,
            schema=schema,
            flux_flag=flux_flag,
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

        return cls(
            ssObjectId,
            filter_list,
            date_range,
            observations_by_filter,
            mpcorb,
            ssobject,
            adler_data,
        )

    @classmethod
    def construct_from_cassandra(
        cls,
        ssObjectId,
        filter_list=["u", "g", "r", "i", "z", "y"],
        date_range=None,
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

        date_range : list of float or None
            Optional. The minimum and maximum dates of the desired observations (MJD), e.g. [60000.0, 67300.0]

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
            obs_filt_mask = observations_table["band"] == filter_name
            if date_range is not None:
                obs_filt_mask = obs_filt_mask & (
                    observations_table["midpointmjdtai"].between(date_range[0], date_range[1])
                )
            obs_slice = observations_table[obs_filt_mask]

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

        return cls(
            ssObjectId,
            filter_list,
            date_range,
            observations_by_filter,
            mpcorb,
            ssobject,
            adler_data,
        )

    @classmethod
    def construct_from_RSP(
        cls,
        ssObjectId,
        filter_list=["u", "g", "r", "i", "z", "y"],
        date_range=None,
        schema="dp03_catalogs_10yr",
        api_token_path=None,
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

        date_range : list of float or None
            Optional. The minimum and maximum dates of the desired observations (MJD), e.g. [60000.0, 67300.0]

        schema : str or None
            Schema/database from which to select the data tables. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        api_token_path : str or None
            Path to user RSP API token if running not on RSP. See https://rsp.lsst.io/guides/auth/creating-user-tokens.html and lsst-adler/notebooks/adler_demo/adler_demo_rsp_api.ipynb for guide on setting this up.

        flux_flag : str or None
            Name of the flux column to select from DP1 DiaSource table. Determines FluxErr and ra/dec columns to select also. Default is None (selects mag/magErr/ra/dec for DP0.3)

        """

        if date_range is not None:
            if len(date_range) != 2:
                logger.error("ValueError: date_range attribute must be of length 2.")
                raise ValueError("date_range attribute must be of length 2.")

        rsp_tap_path = RSP_TAP_CONFIG_DICT[schema]  # TODO give better name

        # Select correct TAP service depending on schema chosen
        if api_token_path:
            service = get_tap_service_api(rsp_tap_path, api_token_path=api_token_path)
        else:
            # service = get_tap_service(rsp_tap_path)
            discovery = RSPDiscovery(rsp_tap_path)
            service = discovery.get_tap_client()

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

        return cls(
            ssObjectId,
            filter_list,
            date_range,
            observations_by_filter,
            mpcorb,
            ssobject,
            adler_data,
        )

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

        date_range : list of float or None
            Optional. The minimum and maximum dates of the desired observations (MJD), e.g. [60000.0, 67300.0]

        service : pyvo.dal.tap.TAPService object or None
            TAPService object linked to the RSP. Default=None.

        sql_filename : str or None
            Filepath to a SQL database. Default=None.

        schema : str
            Schema/database from which to select the data tables. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        flux_flag : str or None
            Name of the flux column to select from DP1 DiaSource table. Determines FluxErr and ra/dec columns to select also. Default is None (selects mag/magErr/ra/dec for DP0.3)
        """

        if sql_filename:
            sql_schema = ""
        else:  # pragma: no cover
            sql_schema = schema + "."

        try:
            selected_config = SCHEMA_CONFIG_DICT[schema][flux_flag]
        except KeyError:
            if schema not in SCHEMA_CONFIG_DICT:
                logger.error(f"Schema {schema} not recognised.")
                raise Exception(f"Schema {schema} not recognised.")
            else:
                logger.error(f"Flux column {flux_flag} not recognised for schema {schema}.")
                raise Exception(f"Flux column {flux_flag} not recognised for schema {schema}.")

        fluxmag_column = selected_config["fluxmag_column"]
        fluxmag_err_column = selected_config["fluxmag_err_column"]
        fluxunit = selected_config["fluxunit"]
        ra_column = selected_config["ra_column"]
        dec_column = selected_config["dec_column"]

        observations_by_filter = []

        for filter_name in filter_list:
            observations_sql_query = f"""
                SELECT
                    SSObject.ssObjectId, SSSource.diaSourceId, {fluxmag_column}, {fluxmag_err_column}, band, {ADLER_SCHEMA[schema]['midpointMjdTai']} AS midpointMjdTai, {ra_column} AS ra, {dec_column} AS dec, {ADLER_SCHEMA[schema]['phaseAngle']} AS phaseAngle,
                    {ADLER_SCHEMA[schema]['topocentricDist']} AS topocentricDist, {ADLER_SCHEMA[schema]['heliocentricDist']} AS heliocentricDist, {ADLER_SCHEMA[schema]['heliocentricX']} AS heliocentricX, {ADLER_SCHEMA[schema]['heliocentricY']} AS heliocentricY, {ADLER_SCHEMA[schema]['heliocentricZ']} AS heliocentricZ,
                    {ADLER_SCHEMA[schema]['topocentricX']} AS topocentricX, {ADLER_SCHEMA[schema]['topocentricY']} AS topocentricY, {ADLER_SCHEMA[schema]['topocentricZ']} AS topocentricZ,
                    {ADLER_SCHEMA[schema]['eclipticLambda']} AS eclipticLambda, {ADLER_SCHEMA[schema]['eclipticBeta']} AS eclipticBeta 
                FROM
                    {sql_schema}SSObject
                    JOIN {sql_schema}DiaSource ON {sql_schema}SSObject.ssObjectId   = {sql_schema}DiaSource.ssObjectId
                    JOIN {sql_schema}SSSource  ON {sql_schema}DiaSource.diaSourceId = {sql_schema}SSSource.diaSourceId
                WHERE
                    SSObject.ssObjectId = {ssObjectId} AND band = '{filter_name}'
                """
            # TODO: log the query

            if date_range is not None:
                observations_sql_query += f" AND midPointMjdTai BETWEEN {date_range[0]} AND {date_range[1]}"

            # This function submits the query and gets the results (or pulls from the SQL database)
            data_table = get_data_table(observations_sql_query, service=service, sql_filename=sql_filename)

            # check for any observations, skip if none available
            if len(data_table) == 0:
                logger.warning(
                    "No observations found in {} filter for this object. Skipping this filter.".format(
                        filter_name
                    )
                )
            else:

                # Convert to astropy table so we can operate on it and add mag,magErr columns
                # TODO: temporary fix, get_data_table returns two possible objects (DALResultsTable or Pandas dataframe) that need different handling to convert to astropy tables

                # if isinstance(data_table, pd.DataFrame):
                #     # data_table is DataFrame
                #     data_table_astropy = Table.from_pandas(data_table)
                # else:
                #     # data_table is not Dataframe (DALResultsTable/TAPResults?)
                #     data_table_astropy = data_table.to_table()

                # add the mag and magErr columns if missing
                if ("mag" not in data_table.colnames) & ("magErr" not in data_table.colnames):

                    # ensure flux columns are correct units
                    # Only add units if required, use the fluxunit key to SCHEMA_CONFIG_DICT
                    for x in [fluxmag_column, fluxmag_err_column]:
                        if hasattr(data_table[x], "unit"):
                            print(type(data_table[x]))
                            print(data_table[x].unit, fluxunit)
                            # if the column has units, make sure they are correct as defined in SCHEMA_CONFIG_DICT
                            if (data_table[x].unit != fluxunit) & (data_table[x].unit is not None):
                                data_table[x] = data_table[x].to(fluxunit)
                            else:
                                # catch when the table column is dimensionless
                                data_table[x] = data_table[x].value * fluxunit
                        else:
                            # if the column does not have units, add them
                            data_table[x] *= fluxunit

                    # Compute magnitudes
                    mag, mag_err = flux_to_magnitude(
                        data_table[fluxmag_column], data_table[fluxmag_err_column]
                    )

                    # Insert the new columns at the same positions
                    data_table.add_column(mag, name="mag", index=data_table.colnames.index(fluxmag_column))
                    data_table.add_column(
                        mag_err, name="magErr", index=data_table.colnames.index(fluxmag_err_column)
                    )

                    # Remove the old flux columns
                    data_table.remove_columns([fluxmag_column, fluxmag_err_column])

                # TODO: If mag columns already exist (dp03) do we need to add mag units?

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
            Schema/database from which to select the data tables. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        """

        if sql_filename:
            sql_schema = ""
        else:  # pragma: no cover
            sql_schema = schema + "."

        if schema in SCHEMA_CONFIG_DICT:

            mpc_id = SCHEMA_CONFIG_DICT[schema]["mpc_id"]
            mpc_table = SCHEMA_CONFIG_DICT[schema]["mpc_table"]

            # determine the query_id (ssObjectId or designation) and construct the constraint
            if (mpc_id != "ssObjectId") & (schema!="MPC"):
                constraint = f"JOIN {schema}.SSObject AS sso ON mpc.{mpc_id} = sso.{mpc_id} WHERE sso.ssObjectId={ssObjectId}"
            else:
                constraint = f"WHERE mpc.{mpc_id} = '{ssObjectId}'"

            # get the MPCORB field names from the schema
            query_fields = [
                x
                for x in ADLER_SCHEMA[schema]
                if (ADLER_SCHEMA[schema + "_table"][x] == mpc_table)
                & (not pd.isnull(ADLER_SCHEMA[schema][x]))
                & (x != mpc_id)
            ]
            query_fields = ["mpc.{} AS {}".format(ADLER_SCHEMA[schema][x], x) for x in query_fields]
            query_fields = ["mpc.{}".format(mpc_id)] + query_fields

            MPCORB_sql_query = f"""
                SELECT
                    {",".join(query_fields)}
                FROM
                    {sql_schema}{mpc_table} as mpc
                {constraint}
            """
            # TODO: log the query
            print(MPCORB_sql_query)
        else:
            logger.error(f"Schema {schema} not recognised.")
            raise Exception(f"Schema {schema} not recognised.")

        data_table = get_data_table(MPCORB_sql_query, service=service, sql_filename=sql_filename)
        print(data_table)
        
        if len(data_table) == 0:
            err_message = "No {} data for this object could be found for {}={}.".format(
                mpc_table, mpc_id, query_id
            )
            logger.error(err_message)
            raise Exception(err_message)

        # Any required values are added below if missing
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

        schema : str
            Schema/database from which to select the data tables. Default is currently "dp03_catalogs_10yr" for testing using DP0.3.

        """

        if sql_filename:
            sql_schema = ""
        else:  # pragma: no cover
            sql_schema = schema + "."

        if schema in SCHEMA_CONFIG_DICT:

            # if schema == "MPC":
            #     ssobject_id_col = "provid"
            #     ssobject_table = "obs_sbn"
            #     constraint = f"{ssobject_id_col} = '{ssObjectId}' LIMIT 1" # TODO: a bit unnecessary but we only need to get one row to check the object is indeed in the table - try replace MPC case with just a blank SSObject!
            # else:
            #     ssobject_id_col = "ssObjectId"
            #     ssobject_table = "SSObject"
            #     constraint = f"{ssobject_id_col} = '{ssObjectId}'

            # get the SSObject field names from the schema
            query_fields = [
                x
                for x in ADLER_SCHEMA[schema]
                if (ADLER_SCHEMA[schema + "_table"][x] == "SSObject")
                & (not pd.isnull(ADLER_SCHEMA[schema][x]))
            ]
            query_fields = ["{} AS {}".format(ADLER_SCHEMA[schema][x], x) for x in query_fields]

            # seperate out the filter dependent columns
            filter_fields = [x for x in query_fields if x.startswith("{filt}")]
            query_fields = [x for x in query_fields if not x.startswith("{filt}")]

            # add columns for required for each filter
            filter_dependent_columns = []
            for filter_name in filter_list:
                filter_fields_list = [x.replace("{filt}", filter_name) for x in filter_fields]
                filter_dependent_columns += filter_fields_list
            query_fields += filter_dependent_columns

            query_fields = ["ssObjectId"] + query_fields

            SSObject_sql_query = f"""
                SELECT
                    {",".join(query_fields)}
                FROM
                    {sql_schema}SSObject
                WHERE
                    ssObjectId = '{ssObjectId}'
            """
            # TODO: log the query
            print(SSObject_sql_query)
        else:
            logger.error(f"Schema {schema} not recognised.")
            raise Exception(f"Schema {schema} not recognised.")

        data_table = get_data_table(SSObject_sql_query, service=service, sql_filename=sql_filename)
        print(data_table)
        print(len(data_table))
        
        if len(data_table) == 0:
            err_message = "No SSObject data for this object could be found for ssObjectId={}.".format(
                ssObjectId
            )
            logger.error(err_message)
            raise Exception(err_message)

        return SSObject.construct_from_data_table(ssObjectId, filter_list, data_table)

    @classmethod
    def construct_from_mpc_obs_sbn(
        cls,
        ssObjectId,
        sql_filename,
        filter_list=["u", "g", "r", "i", "z", "y"],
        date_range=None,
    ):
        """Custom constructor which builds the AdlerPlanetoid object and the associated Observations, MPCORB and SSObject objects
        from the MPC obs_sbn database. This is designed specifically for the SSSC Prompt Products Database Bandaid.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        sql_filename : str
            Filepath to the local SQL database.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of float or None
            Optional. The minimum and maximum dates of the desired observations (MJD), e.g. [60000.0, 67300.0]

        """

        if date_range is not None:
            if len(date_range) != 2:
                logger.error("ValueError: date_range attribute must be of length 2.")
                raise ValueError("date_range attribute must be of length 2.")

        observations_by_filter = cls.populate_observations_from_mpc_obs_sbn(
            cls, ssObjectId, filter_list, date_range, sql_filename=sql_filename
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

        # mpcorb = cls.populate_MPCORB_from_mpc_obs_sbn(cls, ssObjectId, sql_filename=sql_filename)
        mpcorb = cls.populate_MPCORB(cls, ssObjectId, sql_filename=sql_filename, schema = "MPC")
        # ssobject = cls.populate_SSObject_from_mpc_obs_sbn(
        #     cls, ssObjectId, filter_list, sql_filename=sql_filename
        # )
        # TODO: make a blank (default values) SSObject and add the numObs from the planetoid Observations
        ssobject = SSObject(ssObjectId,filter_list = filter_list, numObs = sum([len(obs.__dict__) for obs in observations_by_filter]))
        
        adler_data = AdlerData(ssObjectId, filter_list)

        return cls(
            ssObjectId,
            filter_list,
            date_range,
            observations_by_filter,
            mpcorb,
            ssobject,
            adler_data,
        )

    def populate_observations_from_mpc_obs_sbn(self, ssObjectId, filter_list, date_range, sql_filename):
        """Populates the observations_by_filter class attribute. This version is specific to the construct_from_mpc_obs_sbn function.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        date_range : list of float or None
            Optional. The minimum and maximum dates of the desired observations (MJD), e.g. [60000.0, 67300.0]

        sql_filename : str
            Filepath to an SQL database.

        """

        logger.warning(
            f"Constructing from the MPC obs_sbn table populates the following LSST schema columns as their best case obs_sbn analogs (LSST column name = obs_sbn column name):"
        )
        logger.warning(f"SSObjectId = provid; diaSourceId = obsid; magErr = rmsmag")
        logger.warning(f"mjd_utc is converted to mjd_tai and presented as midpointMjdTai")
        logger.warning(
            f"phaseAngle, topocentricDist and heliocentricDist are not currently corrected for light travel time effects"
        )
        logger.warning(
            f"heliocentricX, heliocentricY, heliocentricZ, topocentricX, topocentricY, topocentricZ, eclipticLambda, eclipticBeta are unpopulated and selected as NULLs."
        )

        observations_by_filter = []

        for filter_name in filter_list:
            # TODO: update this query using ADLER_SCHEMA!
            observations_sql_query = f"""
                SELECT
                    provid AS SSObjectId, obsid as diaSourceId, mag, rmsmag AS magErr, band, mjd_tai AS midpointMjdTai, ra, dec,
                    phaseAngle, topocentricDist, heliocentricDist,
                    NULL AS heliocentricX, NULL AS heliocentricY, NULL AS heliocentricZ,
                    NULL AS topocentricX, NULL AS topocentricY, NULL AS topocentricZ,
                    NULL AS eclipticLambda, NULL AS eclipticBeta
                FROM
                    obs_sbn
                WHERE
                    provid='{ssObjectId}' AND band = '{filter_name}'
                """
            if date_range is not None:
                observations_sql_query += f" AND mjd_tai BETWEEN '{date_range[0]}' AND '{date_range[1]}'"

            # This function submits the query and gets the results from the SQL database supplied
            # Explicitly setting service=None here for clarity as this version does not query from non-local databases
            data_table = get_data_table(observations_sql_query, service=None, sql_filename=sql_filename)

            if len(data_table) == 0:
                logger.warning(
                    "No observations found in {} filter for this object. Skipping this filter.".format(
                        filter_name
                    )
                )
            else:
                # DP1 discoveries have no magErr values so fill with NaNs
                # for some reason (TODO: check why) this means that this column has dtype object so we force it to be float64 here
                data_table["magErr"] = data_table["magErr"].astype(
                    np.float64
                )  # TODO: add dtype checking somewhere, populate_observations etc?

                observations_by_filter.append(
                    Observations.construct_from_data_table(ssObjectId, filter_name, data_table)
                )

        return observations_by_filter

    # def populate_MPCORB_from_mpc_obs_sbn(self, ssObjectId, sql_filename):
    #     """Populates the MPCORB object class attribute. This version is specific to the construct_from_mpc_obs_sbn function.

    #     Parameters
    #     -----------
    #     ssObjectId : str
    #         ssObjectId of the object of interest.

    #     sql_filename : str or None
    #         Filepath to an SQL database.

    #     """

    #     logger.warning(
    #         f"Constructing from the MPC obs_sbn table populates the following LSST schema columns as their best case obs_sbn analogs (LSST column name = obs_sbn column name):"
    #     )
    #     logger.warning(f"ssObjectId = fullDesignation; fullDesignation = fullDesignation; tperi = t_p")
    #     logger.warning(
    #         f"mpcDesignation, mpcNumber, mpcG, n, uncertaintyParameter, flags are unpopulated and selected as NULL/0."
    #     )
    #     mpc_orbits_sql_query = f"""
    #         SELECT
    #             fullDesignation AS ssObjectId, NULL AS mpcDesignation, fullDesignation AS fullDesignation, 0 AS mpcNumber,
    #             mpcH, NULL AS mpcG, epoch, t_p AS tperi, peri, node, incl, e, NULL AS n, q, NULL AS uncertaintyParameter, NULL AS flags
    #         FROM
    #             mpc_orbits
    #         WHERE
    #             fullDesignation = '{ssObjectId}'
    #     """

    #     # Explicitly setting service=None here for clarity as this version does not query from non-local databases
    #     data_table = get_data_table(mpc_orbits_sql_query, service=None, sql_filename=sql_filename)
    #     print(data_table)
    #     print(len(data_table))
        
    #     if len(data_table) == 0:
    #         logger.error("No mpc_orbits data for this object could be found for this SSObjectId.")
    #         raise Exception("No mpc_orbits data for this object could be found for this SSObjectId.")

    #     return MPCORB.construct_from_data_table(ssObjectId, data_table)

    # def populate_SSObject_from_mpc_obs_sbn(self, ssObjectId, filter_list, sql_filename):
    #     """Populates the SSObject class attribute. This version is specific to the construct_from_mpc_obs_sbn function.

    #     Parameters
    #     -----------
    #     ssObjectId : str
    #         ssObjectId of the object of interest.

    #     filter_list : list of str
    #         A comma-separated list of the filters of interest.

    #     sql_filename : str or None
    #         Filepath to an SQL database.

    #     """

    #     filter_dependent_columns = ""

    #     for filter_name in filter_list:
    #         # Counting number of observations in given filter in the query here
    #         filter_string = "NULL AS {}_H, NULL AS {}_G12, NULL AS {}_HErr, NULL AS {}_G12Err, (SELECT COUNT(*) FROM obs_sbn WHERE band='{}' and provid='{}') AS {}_Ndata, ".format(
    #             filter_name, filter_name, filter_name, filter_name, filter_name, ssObjectId, filter_name
    #         )

    #         filter_dependent_columns += filter_string

    #     logger.warning(
    #         f"Constructing from the MPC obs_sbn table populates the following LSST schema columns as their best case obs_sbn analogs (LSST column name = obs_sbn column name):"
    #     )
    #     logger.warning(f"All columns other than numObs/'band'_Ndata are selected as NULL/0.")
    #     # TODO: update this query with ADLER_SCHEMA!
    #     SSObject_sql_query = f"""
    #         SELECT
    #             NULL AS discoverySubmissionDate, NULL AS firstObservationDate, NULL AS arc, count(*) AS numObs, 
    #             {filter_dependent_columns}
    #             NULL AS maxExtendedness, NULL AS minExtendedness, NULL AS medianExtendedness
    #         FROM
    #             obs_sbn
    #         WHERE
    #             provid = '{ssObjectId}'
    #     """
    #     print(SSObject_sql_query)

    #     # Explicitly setting service=None here for clarity as this version does not query from non-local databases
    #     data_table = get_data_table(SSObject_sql_query, service=None, sql_filename=sql_filename)
    #     print(data_table)
    #     print(len(data_table))
        
    #     # TODO probably add some warnings for these as there isn't actually any SSObject data for these things in MPC file
    #     if (len(data_table) == 0) or (data_table["numObs"].values == 0):
    #         logger.error("No SSObject data for this object could be found for this SSObjectId.")
    #         raise Exception("No SSObject data for this object could be found for this SSObjectId.")

    #     return SSObject.construct_from_data_table(ssObjectId, filter_list, data_table)

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

    def observations_within_time(self, start=None, stop=None):
        """Get a dataframe of all observations (across all filters) taken within a given time interval.


        Parameters
        ----------
        start, stop : float
            The time limits as modified Julian dates. Optional, if not declared then get all data


        Returns
        -------
        observations : pandas.DataFrame

        """

        result = pd.DataFrame()
        for obs in self.observations_by_filter:
            df = pd.DataFrame(obs.__dict__ | {"filter_name": [obs.filter_name] * obs.num_obs})
            result = pd.concat([result, df]).reset_index(drop=True)

        if start is None:
            start = np.amin(result["midpointMjdTai"])
        if stop is None:
            stop = np.amax(result["midpointMjdTai"])
        i = (result.midpointMjdTai >= start) * (result.midpointMjdTai <= stop)
        result = result[i]
        result = result.sort_values("midpointMjdTai")
        return result

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

    def attach_previous_adler_data(self, filepath, modelId=None):
        """Attaches and returns an AdlerData object containing the most recent AdlerData
        for this ssObjectId.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.

        modelId : str, optional
            modelId for the model of interest that should be recovered. Default: None.

        """

        self.PreviousAdlerData = AdlerData(self.ssObjectId, self.filter_list)
        self.PreviousAdlerData.populate_from_database(filepath, modelId=modelId)

        return self.PreviousAdlerData
