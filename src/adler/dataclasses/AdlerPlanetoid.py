from adler.dataclasses.DataSchema import Observations, MPCORB, SSObject


class AdlerPlanetoid:
    """AdlerPlanetoid class. Contains the Observations, MPCORB and SSObject objects."""

    def __init__(self, ssObjectId, population_location="RSP", sql_filename=None):
        """Initialises the AdlerPlanetoid object.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.
        population_location : str
            String delineating source of data. Should be "RSP" for Rubin Science Platform or "SQL" for a SQL table.
        sql_filename: str, optional
            Location of local SQL database, if using.

        """
        self.ssObjectId = ssObjectId
        self.population_location = population_location
        self.sql_filename = sql_filename
        # can also include date ranges at some point

        # this creates the AdlerPlanetoid.Observations, AdlerPlanetoid.MPCORB and
        # AdlerPlanetoid.SSObject objects.
        self.populate_observations()
        self.populate_MPCORB()
        self.populate_SSObject()

    def populate_observations(self):
        """Populates the Observations object class attribute."""
        observations_sql_query = f"""
            SELECT
                ssObject.ssObjectId, mag, magErr, band, midpointMjdTai as mjd, ra, dec, phaseAngle,
                topocentricDist, heliocentricDist
            FROM
                dp03_catalogs_10yr.ssObject
                JOIN dp03_catalogs_10yr.diaSource ON dp03_catalogs_10yr.ssObject.ssObjectId   = dp03_catalogs_10yr.diaSource.ssObjectId
                JOIN dp03_catalogs_10yr.ssSource  ON dp03_catalogs_10yr.diaSource.diaSourceId = dp03_catalogs_10yr.ssSource.diaSourceId
            WHERE
                ssObject.ssObjectId = {self.ssObjectId} and band='r'
            """

        self.Observations = Observations(
            self.ssObjectId, self.population_location, observations_sql_query, sql_filename=self.sql_filename
        )

    def populate_MPCORB(self):
        """Populates the MPCORB object class attribute."""
        MPCORB_sql_query = f"""
            SELECT
                ssObjectId, mpcDesignation, mpcNumber, mpcH, mpcG, epoch, peri, node, incl, e, n, q, 
                uncertaintyParameter, flags
            FROM
                dp03_catalogs_10yr.MPCORB
            WHERE
                ssObjectId = {self.ssObjectId}
        """

        self.MPCORB = MPCORB(
            self.ssObjectId, self.population_location, MPCORB_sql_query, sql_filename=self.sql_filename
        )

    def populate_SSObject(self):
        """Populates the SSObject class attribute."""
        SSObject_sql_query = f"""
            SELECT
                discoverySubmissionDate, firstObservationDate, arc, numObs, 
                r_H, r_G12, r_Herr, r_G12err, r_nData,
                maxExtendedness, minExtendedness, medianExtendedness
            FROM
                dp03_catalogs_10yr.SSObject
            WHERE
                ssObjectId = {self.ssObjectId}
        """

        self.SSObject = SSObject(
            self.ssObjectId,
            self.population_location,
            sql_query=SSObject_sql_query,
            sql_filename=self.sql_filename,
        )
