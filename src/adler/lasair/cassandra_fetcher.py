import json
import sys
from cassandra.cluster import Cluster, ConsistencyLevel
from cassandra.query import dict_factory, SimpleStatement


class CassandraFetcher:  # pragma: no cover
    """Class to fetch data from a Cassandra database, used for Lasair integration.

    TODO: move to the lasair-adler repo.

    Attributes
    -----------
    cassandra_hosts : list of str
        Location of the Cassandra database - usually an IP address. Default is ["10.21.3.123"].

    """

    def __init__(self, cassandra_hosts):
        self.cluster = Cluster(cassandra_hosts)
        self.session = self.cluster.connect()
        # Set the row_factory to dict_factory, otherwise
        # the data returned will be in the form of object properties.
        self.session.row_factory = dict_factory
        self.session.set_keyspace("adler")

    def fetch_SSObject(self, ssObjectId, filter_list):
        """Fetches the metadata from the SSObject table of a Cassandra database as a dictionary.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        Returns
        -----------
        dict
            A dictionary of metadata for the object of interest in the filters
            of interest.

        """

        filter_dependent_columns = ""
        for filter_name in filter_list:
            filter_string = "{}_H, {}_G12, {}_HErr, {}_G12Err, {}_Ndata, ".format(
                filter_name, filter_name, filter_name, filter_name, filter_name
            )

            filter_dependent_columns += filter_string

        obj = {}

        SSObject_sql_query = f"""
            SELECT
                discoverySubmissionDate, firstObservationDate, arc, numObs, 
                {filter_dependent_columns}
                maxExtendedness, minExtendedness, medianExtendedness
            FROM
                ssobjects
            WHERE
                ssObjectId = {ssObjectId}
        """

        ret = self.session.execute(SSObject_sql_query)

        for ssObject in ret:
            obj = ssObject

        return obj

    def fetch_MPCORB(self, ssObjectId):
        """Fetches the metadata from the MPCORB table of a Cassandra database as a dictionary.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        Returns
        -----------
        dict
            A dictionary of metadata for the object of interest.

        """

        obj = {}

        MPCORB_sql_query = f"""
            SELECT
                ssObjectId, mpcDesignation, fullDesignation, mpcNumber, mpcH, mpcG, epoch, tperi, peri, node, incl, e, n, q, 
                uncertaintyParameter, flags
            FROM
                mpcorbs
            WHERE
                ssObjectId = {ssObjectId}
        """

        ret = self.session.execute(MPCORB_sql_query)

        for MPCORB in ret:
            obj = MPCORB

        return obj

    def fetch_observations(self, ssObjectId):
        """Fetches the source observations from the DIASource and SSSource tables as a dictionary.
        Note that it will retrieve ALL observations for the object regardless of filter and data range,
        so any filtering must be performed later. This is due to restrictions on queries to Cassandra.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        Returns
        -----------
        dict
            A dictionary of metadata for the object of interest in the filters
            of interest.

        """

        sourceDict = {}

        dia_query = f"""
                    SELECT 
                        diasourceid, band, mag, magErr, midPointMjdTai, ra, decl
                    FROM 
                        diasources
                    WHERE
                        ssObjectId = {ssObjectId}
                    """
        ret = self.session.execute(dia_query)

        n = 0
        for diaSource in ret:
            sourceDict[diaSource["diasourceid"]] = diaSource
            n += 1

        ss_query = f"""SELECT diasourceid, phaseAngle, topocentricDist, heliocentricDist, heliocentricX, heliocentricY, heliocentricZ,
                            topocentricX, topocentricY, topocentricZ, eclipticLambda, eclipticBeta
                        FROM sssources
                        WHERE
                        ssObjectId = {ssObjectId}
                        """
        ret = self.session.execute(ss_query)

        n = 0
        for ssSource in ret:
            n += 1
            sourceDict[ssSource["diasourceid"]].update(ssSource)

        sources = []
        for k, v in sourceDict.items():
            sources.append(v)

        return sources
